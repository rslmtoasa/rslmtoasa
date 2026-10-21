! Opt-in diagnostics, including an explicit saved-input counterfactual.
! No exchange or Hamiltonian formula is implemented here.
module gf_audit_io_mod
   use precision_mod, only: rp
   implicit none
   private
   public :: gf_audit_selected, gf_audit_stride, gf_audit_matrix, gf_audit_open, gf_audit_h_enabled
   public :: gf_audit_load_sbar
   public :: gf_audit_graph
contains
   subroutine gf_audit_graph(nn, types, h, ho, overlap, enu, soc)
      ! Exhaustive bare-h/weighted-h graph check; no Hamiltonian is modified.
      ! h=h^dagger, O=O^dagger and ho=h*O imply Hermiticity of h-h*O*h+enu+soc.
      integer, intent(in) :: nn(:, :), types(:)
      complex(rp), intent(in) :: h(:,:,:,:), ho(:,:,:,:), overlap(:,:,:), enu(:,:,:), soc(:,:,:)
      complex(rp) :: defect(18,18)
      real(rp) :: hmax, hfro2, homax, ofail, efail, sfail
      integer :: i,j,nb,back,k,ti,tj,hits,edges,u
      integer :: override_status
      character(4096) :: override_path
      if (.not. gf_audit_h_enabled()) return
      hmax=0; hfro2=0; homax=0; ofail=0; efail=0; sfail=0; edges=0
      do ti=1,size(overlap,3)
         ofail=max(ofail,maxval(abs(overlap(:,:,ti)-conjg(transpose(overlap(:,:,ti))))))
         efail=max(efail,maxval(abs(enu(:,:,ti)-conjg(transpose(enu(:,:,ti))))))
         sfail=max(sfail,maxval(abs(soc(:,:,ti)-conjg(transpose(soc(:,:,ti))))))
      end do
      do i=1,size(types)
         ti=types(i)
         do nb=1,nn(i,1)
            j=i
            if (nb > 1) j=nn(i,nb)
            if (j == 0) cycle
            tj=types(j);back=1;hits=1
            if (j /= i) then
               hits=0
               do k=2,nn(j,1)
                  if (nn(j,k) /= i) cycle
                  hits=hits+1;back=k
               end do
            end if
            ! nn can retain unscreened geometric slots with identically zero h.
            ! Such a slot is not a Hamiltonian edge; any nonzero reverse is still checked.
            if (hits == 0 .and. maxval(abs(h(:,:,nb,ti))) == 0.0_rp) cycle
            if (hits /= 1) then
               write(*,*) 'Audit graph edge:',i,j,nb,' reverse hits:',hits,' h scale:',maxval(abs(h(:,:,nb,ti)))
               error stop 'Audit graph missing or ambiguous nonzero reverse edge'
            end if
            defect=h(:,:,nb,ti)-conjg(transpose(h(:,:,back,tj)))
            hmax=max(hmax,maxval(abs(defect)))
            hfro2=hfro2+sum(abs(defect)**2)
            homax=max(homax,maxval(abs(ho(:,:,nb,ti)-matmul(h(:,:,nb,ti),overlap(:,:,tj)))))
            edges=edges+1
         end do
      end do
      open(newunit=u,file='gf-audit-graph.dat',status='replace',action='write')
      write(u,'(a,1x,i0)') 'directed_blocks',edges
      write(u,'(a,1x,es25.16e3)') 'bare_h_max_Ry',hmax
      write(u,'(a,1x,es25.16e3)') 'bare_h_fro_Ry',sqrt(hfro2)
      write(u,'(a,1x,es25.16e3)') 'hO_identity_max',homax
      write(u,'(a,1x,es25.16e3)') 'O_hermiticity_max_per_Ry',ofail
      write(u,'(a,1x,es25.16e3)') 'enu_hermiticity_max_Ry',efail
      write(u,'(a,1x,es25.16e3)') 'soc_hermiticity_max_Ry',sfail
      close(u)
      call get_environment_variable('RSLMTO_AUDIT_SBAR',override_path,status=override_status)
      if (override_status == 0) then
         if (max(hmax,homax,ofail,efail,sfail) > 1.0e-12_rp) &
            error stop 'Reciprocal Sbar failed full graph Hermiticity prerequisites before recursion'
      end if
   end subroutine gf_audit_graph

   subroutine gf_audit_load_sbar(sbar)
      ! Called only by frozen-potential exchange, after the normal structb call.
      ! The normal calculation regenerates Sbar; changing sbar. alone has no effect.
      ! Validate saved versus freshly generated vector/type ORDER before loading.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      complex(rp), intent(inout) :: sbar(:, :, :, :)
      complex(rp), allocatable :: replacement(:, :, :, :)
      character(4096) :: filename, saved_view, line, saved_line, audit
      real(rp) :: row(9)
      integer :: status, unit, generated, labels, ios, saved_ios, typ, pos, k, blocks, nbytes
      integer :: counts(size(sbar, 4))
      call get_environment_variable('RSLMTO_AUDIT_SBAR', filename, status=status)
      if (status == 1) return
      if (status /= 0 .or. len_trim(filename) == 0) error stop 'Invalid audit Sbar path'
      call get_environment_variable('RSLMTO_GF_AUDIT', audit, status=status)
      if (status /= 0 .or. len_trim(audit) == 0) error stop 'Sbar override requires GF audit'
      call get_environment_variable('RSLMTO_AUDIT_SBAR_VIEW', saved_view, status=status)
      if (status /= 0 .or. len_trim(saved_view) == 0) error stop 'Sbar override requires saved vector labels'
      if (size(sbar,1) /= 9 .or. size(sbar,2) /= 9) error stop 'Audit requires nine orbitals'
      flush(16)
      open(newunit=generated, file='view.sbar', status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'Cannot read freshly generated Sbar vector labels'
      open(newunit=labels, file=trim(saved_view), status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'Cannot read saved Sbar vector labels'
      open(newunit=unit, file=trim(filename), form='unformatted', status='old', action='read', iostat=ios)
      if (ios /= 0) error stop 'Cannot read audit Sbar binary'
      allocate(replacement, source=sbar)
      counts=0
      blocks=0
      do
         read(generated, '(a)', iostat=ios) line
         read(labels, '(a)', iostat=saved_ios) saved_line
         if (ios < 0 .and. saved_ios < 0) exit
         if (ios /= 0 .or. saved_ios /= 0) error stop 'Sbar vector files have different lengths'
         if (index(line,'VECTOR=') == 0) cycle
         if (trim(line) /= trim(saved_line)) error stop 'Sbar vector/type/order mismatch'
         pos=index(line,'ICLUS=')
         if (pos == 0) error stop 'Missing Sbar source type'
         read(line(pos+6:), *, iostat=ios) typ
         if (ios /= 0) error stop 'Invalid Sbar source type'
         if (typ < 1 .or. typ > size(sbar,4)) error stop 'Audit source is not a type index'
         counts(typ)=counts(typ)+1
         if (counts(typ) > size(sbar,3)) error stop 'Too many Sbar blocks'
         do k=1,9
            read(unit, iostat=ios) row
            if (ios /= 0) error stop 'Invalid audit Sbar record'
            if (.not. all(ieee_is_finite(row))) error stop 'Nonfinite audit Sbar'
            replacement(k,:,counts(typ),typ)=cmplx(row,0.0_rp,rp)
         end do
         blocks=blocks+1
      end do
      inquire(unit=unit, size=nbytes)
      if (blocks == 0 .or. any(counts == 0)) error stop 'Missing audit Sbar type blocks'
      if (nbytes /= blocks*9*80) error stop 'Unsupported audit Sbar record size/count'
      read(unit, iostat=ios) row
      if (ios >= 0) error stop 'Trailing audit Sbar record'
      close(unit)
      close(labels)
      close(generated)
      sbar=replacement
      write(*,'(a,i0,a,a)') 'DIAGNOSTIC frozen Sbar override: ',blocks,' blocks from ',trim(filename)
   end subroutine gf_audit_load_sbar

   logical function gf_audit_selected(pair)
      integer, intent(in) :: pair
      character(32) :: value
      integer :: status, selected, ios
      gf_audit_selected = .false.
      call get_environment_variable('RSLMTO_GF_AUDIT', value, status=status)
      if (status /= 0) return
      if (trim(value) == 'all') then
         gf_audit_selected = pair >= 1 .and. pair <= 24
      else
         read(value, *, iostat=ios) selected
         if (ios == 0) gf_audit_selected = pair == selected .and. pair >= 1
      end if
   end function

   logical function gf_audit_h_enabled()
      character(8) :: value
      integer :: status
      call get_environment_variable('RSLMTO_GF_AUDIT_H', value, status=status)
      gf_audit_h_enabled = status == 0 .and. trim(value) == '1'
   end function

   integer function gf_audit_stride(ne)
      integer, intent(in) :: ne
      integer :: status, step, ios
      character(32) :: value
      gf_audit_stride = max(1, (ne + 62)/63)
      call get_environment_variable('RSLMTO_GF_AUDIT_STRIDE', value, status=status)
      if (status /= 0) return
      read(value, *, iostat=ios) step
      if (ios == 0) gf_audit_stride = max(gf_audit_stride, step)
   end function

   subroutine gf_audit_open(unit, kind, owner, pair, append)
      integer, intent(out) :: unit
      character(*), intent(in) :: kind
      integer, intent(in) :: owner, pair
      logical, optional, intent(in) :: append
      integer :: ios
      character(160) :: filename
      logical :: add
      add = .false.
      if (present(append)) add = append
      write(filename, '(a,a,a,i0,a,i0,a)') 'gf-audit-',kind,'-r',owner,'-p',pair,'.dat'
      if (add) then
         open(newunit=unit, file=trim(filename), status='old', position='append', action='write', iostat=ios)
      else
         open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
      end if
      if (ios /= 0) then
         unit = 0
         write(*, '(a,a)') 'GF audit could not open ', trim(filename)
      end if
   end subroutine

   subroutine gf_audit_matrix(unit, index, energy, label, matrix)
      integer, intent(in) :: unit, index
      real(rp), intent(in) :: energy
      character(*), intent(in) :: label
      complex(rp), intent(in) :: matrix(:, :)
      integer :: row, col
      write(unit, '(a,1x,i0,1x,es25.16e3,1x,a,1x,i0)') 'M',index,energy,label,size(matrix,1)
      do row = 1, size(matrix, 1)
         write(unit, '(*(es25.16e3,1x))') (real(matrix(row,col),rp), aimag(matrix(row,col)), col=1,size(matrix,2))
      end do
   end subroutine
end module gf_audit_io_mod
