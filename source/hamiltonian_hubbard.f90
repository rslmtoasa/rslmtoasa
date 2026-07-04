submodule(hamiltonian_mod) hamiltonian_hubbard

contains

   module subroutine calculate_hubbard_u_potential_general(this)
      class(hamiltonian), intent(inout) :: this

      integer :: na, l, ispin, i, j, m1, m2, m3, m4, l_index, m_max
      integer :: m1_val, m2_val, m3_val, m4_val
      real(rp) :: f0, f2, f4, f6
      real(rp) :: ubar, jbar, ueff, d1, d2, eps_den, tr_n1mn
      real(rp) :: num_u, num_j, sum_occ_opposite, sum_occ_same_excl
      real(rp) :: common_pref, sum_u_aux, sum_j_aux, dUdn, dJdn, dUeff_dn
      real(rp) :: vdiag_up_avg, vdiag_dn_avg, vdiag_split
      real(rp), dimension(2, 7) :: occ_m
      real(rp), dimension(2, 7, 7) :: pbar
      logical :: use_acbn0
      real(rp), dimension(this%lattice%ntype, 4) :: hub_u, hub_j
      real(rp), dimension(this%lattice%ntype, 4, 4) :: f
      real(rp), dimension(this%lattice%ntype, 4, 2, 7, 7) :: ldm, hub_pot
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin
      real(rp), dimension(this%lattice%ntype, 4) :: n_tot
      type :: int_array
         integer, allocatable :: val(:)
      end type int_array
      type(int_array), dimension(this%lattice%ntype) :: l_arr
      type(int_array), dimension(4) :: ms
      integer :: cntr
      real(rp) :: ldm_trace_up, ldm_trace_dn, ldm_max_abs, hup_max_abs
      real(rp) :: ql_occ_up, ql_occ_dn

      ms(1)%val = [0]
      ms(2)%val = [-1, 0, 1]
      ms(3)%val = [-2, -1, 0, 1, 2]
      ms(4)%val = [-3, -2, -1, 0, 1, 2, 3]

      this%hubbard_u_pot(:, :, :) = 0.0_rp
      hub_u(:, :) = 0.0_rp
      hub_j(:, :) = 0.0_rp
      f(:, :, :) = 0.0_rp
      ldm(:, :, :, :, :) = 0.0_rp
      hub_pot(:, :, :, :, :) = 0.0_rp
      n_spin(:, :, :) = 0.0_rp
      n_tot(:, :) = 0.0_rp

      do na = 1, this%lattice%ntype
         do l = 1, min(4, size(this%lattice%symbolic_atoms(na)%potential%hubbard_u))
            hub_u(na, l) = this%lattice%symbolic_atoms(na)%potential%hubbard_u(l)
         end do
         do l = 1, min(4, size(this%lattice%symbolic_atoms(na)%potential%hubbard_j))
            hub_j(na, l) = this%lattice%symbolic_atoms(na)%potential%hubbard_j(l)
         end do
      end do

      do na = 1, this%lattice%ntype
         f(na, 1, 1) = hub_u(na, 1)
         f(na, 2, 1) = hub_u(na, 2)
         f(na, 2, 2) = hub_j(na, 2)*5.0_rp
         f(na, 3, 1) = hub_u(na, 3)
         f(na, 3, 2) = 14.0_rp*hub_j(na, 3)/1.625_rp
         f(na, 3, 3) = 0.625_rp*f(na, 3, 2)
         f(na, 4, 1) = hub_u(na, 4)
         f(na, 4, 2) = 6435.0_rp*hub_j(na, 4)/(286.0_rp + 195.0_rp*0.67_rp + 250.0_rp*0.49_rp)
         f(na, 4, 3) = 0.67_rp*f(na, 4, 2)
         f(na, 4, 4) = 0.49_rp*f(na, 4, 2)
      end do

      do na = 1, this%lattice%ntype
         cntr = count((abs(hub_u(na, :)) > 1.0e-10_rp) .or. (abs(hub_j(na, :)) > 1.0e-10_rp))
         allocate(l_arr(na)%val(cntr))
         cntr = 0
         do l = 1, 4
            if ((abs(hub_u(na, l)) > 1.0e-10_rp) .or. (abs(hub_j(na, l)) > 1.0e-10_rp)) then
               cntr = cntr + 1
               l_arr(na)%val(cntr) = l
            end if
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 0, min(3, this%control%lmax)
            do ispin = 1, 2
               do i = 1, 2*l + 1
                  do j = 1, 2*l + 1
                     ldm(na, l + 1, ispin, i, j) = this%lattice%symbolic_atoms(na)%potential%ldm(l + 1, ispin, i, j)
                  end do
               end do
            end do
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 1, min(4, this%control%lmax + 1)
            if ((abs(hub_u(na, l)) <= 1.0e-10_rp) .and. (abs(hub_j(na, l)) <= 1.0e-10_rp)) cycle
            ldm_trace_up = 0.0_rp
            ldm_trace_dn = 0.0_rp
            ldm_max_abs = 0.0_rp
            do i = 1, 2*l - 1
               ldm_trace_up = ldm_trace_up + ldm(na, l, 1, i, i)
               ldm_trace_dn = ldm_trace_dn + ldm(na, l, 2, i, i)
               do j = 1, 2*l - 1
                  ldm_max_abs = max(ldm_max_abs, abs(ldm(na, l, 1, i, j)))
                  ldm_max_abs = max(ldm_max_abs, abs(ldm(na, l, 2, i, j)))
               end do
            end do
            if (rank == 0) then
               call g_logger%info('HUBBARD_LDM type='//fmt('i4', na)//' l='//fmt('i2', l - 1)// &
                  ' tr_up='//fmt('f10.6', ldm_trace_up)//' tr_dn='//fmt('f10.6', ldm_trace_dn)// &
                  ' maxabs='//fmt('es12.4', ldm_max_abs), __FILE__, __LINE__)
               call g_logger%info('HUBBARD_OCC type='//fmt('i4', na)//' l='//fmt('i2', l - 1)// &
                  ' U='//fmt('f10.6', hub_u(na, l))//' J='//fmt('f10.6', hub_j(na, l))// &
                  ' nup='//fmt('f10.6', ldm_trace_up)//' ndn='//fmt('f10.6', ldm_trace_dn)// &
                  ' ntot='//fmt('f10.6', ldm_trace_up + ldm_trace_dn), __FILE__, __LINE__)
            end if
            ql_occ_up = 0.0_rp
            ql_occ_dn = 0.0_rp
            if ((l - 1) <= this%lattice%symbolic_atoms(na)%potential%lmax) then
               ql_occ_up = this%lattice%symbolic_atoms(na)%potential%ql(1, l - 1, 1)
               ql_occ_dn = this%lattice%symbolic_atoms(na)%potential%ql(1, l - 1, 2)
            end if
            if (rank == 0) then
               call g_logger%info('HUBBARD_QLCMP type='//fmt('i4', na)//' l='//fmt('i2', l - 1)// &
                  ' ql_up='//fmt('f10.6', ql_occ_up)//' ql_dn='//fmt('f10.6', ql_occ_dn)// &
                  ' ldm_up='//fmt('f10.6', ldm_trace_up)//' ldm_dn='//fmt('f10.6', ldm_trace_dn)// &
                  ' d_up(ql-ldm)='//fmt('f10.6', ql_occ_up - ldm_trace_up)// &
                  ' d_dn(ql-ldm)='//fmt('f10.6', ql_occ_dn - ldm_trace_dn), __FILE__, __LINE__)
            end if
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 1, 4
            do ispin = 1, 2
               do i = 1, 2*l - 1
                  n_spin(na, l, ispin) = n_spin(na, l, ispin) + ldm(na, l, ispin, i, i)
               end do
            end do
            ! Total occupancy for this channel is n_up + n_dn.
            n_tot(na, l) = n_spin(na, l, 1) + n_spin(na, l, 2)
         end do
      end do

      use_acbn0 = (trim(lower(this%hubbard_u_potential_form)) == 'acbn0')
      eps_den = 1.0e-12_rp

      do na = 1, this%lattice%ntype
         do l = 1, size(l_arr(na)%val)
            l_index = l_arr(na)%val(l)
            m_max = 2*l_index - 1

            if (use_acbn0) then
               occ_m(:, :) = 0.0_rp
               pbar(:, :, :) = 0.0_rp
               do ispin = 1, 2
                  do m1 = 1, m_max
                     occ_m(ispin, m1) = ldm(na, l_index, ispin, m1, m1)
                  end do
               end do
               do ispin = 1, 2
                  do m1 = 1, m_max
                     do m2 = 1, m_max
                        pbar(ispin, m1, m2) = (occ_m(ispin, m1) + occ_m(ispin, m2))*ldm(na, l_index, ispin, m1, m2)
                     end do
                  end do
               end do

               d1 = 0.0_rp
               d2 = 0.0_rp
               do m1 = 1, m_max
                  do m2 = 1, m_max
                     d1 = d1 + occ_m(1, m1)*occ_m(2, m2) + occ_m(2, m1)*occ_m(1, m2)
                     if (m1 /= m2) then
                        d1 = d1 + occ_m(1, m1)*occ_m(1, m2) + occ_m(2, m1)*occ_m(2, m2)
                        d2 = d2 + occ_m(1, m1)*occ_m(1, m2) + occ_m(2, m1)*occ_m(2, m2)
                     end if
                  end do
               end do

               f0 = f(na, l_index, 1)
               f2 = f(na, l_index, 2)
               f4 = f(na, l_index, 3)
               f6 = f(na, l_index, 4)
               num_u = 0.0_rp
               num_j = 0.0_rp
               do ispin = 1, 2
                  do i = 1, 2
                     do m1 = 1, m_max
                        do m2 = 1, m_max
                           do m3 = 1, m_max
                              do m4 = 1, m_max
                                 m1_val = ms(l_index)%val(m1)
                                 m2_val = ms(l_index)%val(m2)
                                 m3_val = ms(l_index)%val(m3)
                                 m4_val = ms(l_index)%val(m4)
                                 num_u = num_u + Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6)* &
                                         pbar(ispin, m1, m2)*pbar(i, m3, m4)
                              end do
                           end do
                        end do
                     end do
                  end do
                  do m1 = 1, m_max
                     do m2 = 1, m_max
                        if (m1 == m2) cycle
                        do m3 = 1, m_max
                           do m4 = 1, m_max
                              m1_val = ms(l_index)%val(m1)
                              m2_val = ms(l_index)%val(m2)
                              m3_val = ms(l_index)%val(m3)
                              m4_val = ms(l_index)%val(m4)
                              num_j = num_j + Coulomb_mat(l_index - 1, m1_val, m3_val, m4_val, m2_val, f0, f2, f4, f6)* &
                                      pbar(ispin, m1, m2)*pbar(ispin, m3, m4)
                           end do
                        end do
                     end do
                  end do
               end do
               ubar = num_u/max(d1, eps_den)
               jbar = num_j/max(d2, eps_den)
               ueff = ubar - jbar
               tr_n1mn = 0.0_rp
               do ispin = 1, 2
                  do m1 = 1, m_max
                     tr_n1mn = tr_n1mn + occ_m(ispin, m1)*(1.0_rp - occ_m(ispin, m1))
                  end do
               end do
            end if

            do ispin = 1, 2
               do m1 = 1, m_max
                  do m2 = 1, m_max
                     do m3 = 1, m_max
                        do m4 = 1, m_max
                           m1_val = ms(l_index)%val(m1)
                           m2_val = ms(l_index)%val(m2)
                           m3_val = ms(l_index)%val(m3)
                           m4_val = ms(l_index)%val(m4)
                           f0 = f(na, l_index, 1)
                           f2 = f(na, l_index, 2)
                           f4 = f(na, l_index, 3)
                           f6 = f(na, l_index, 4)
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                                + Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6)*ldm(na, l_index, 3 - ispin, m3, m4) &
                                + (Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6) &
                                -  Coulomb_mat(l_index - 1, m1_val, m3_val, m4_val, m2_val, f0, f2, f4, f6))*ldm(na, l_index, ispin, m3, m4)
                        end do
                     end do
                     if (use_acbn0) then
                        sum_u_aux = 0.0_rp
                        sum_j_aux = 0.0_rp
                        m1_val = ms(l_index)%val(m1)
                        m2_val = ms(l_index)%val(m2)
                        do j = 1, 2
                           do m3 = 1, m_max
                              do m4 = 1, m_max
                                 m3_val = ms(l_index)%val(m3)
                                 m4_val = ms(l_index)%val(m4)
                                 sum_u_aux = sum_u_aux + Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6)*pbar(j, m3, m4)
                              end do
                           end do
                        end do
                        do m3 = 1, m_max
                           do m4 = 1, m_max
                              m3_val = ms(l_index)%val(m3)
                              m4_val = ms(l_index)%val(m4)
                              sum_j_aux = sum_j_aux + Coulomb_mat(l_index - 1, m1_val, m3_val, m4_val, m2_val, f0, f2, f4, f6)*pbar(ispin, m3, m4)
                           end do
                        end do

                        common_pref = occ_m(ispin, m1) + occ_m(ispin, m2)
                        if (m1 == m2) common_pref = common_pref + 2.0_rp*occ_m(ispin, m1)
                        sum_occ_opposite = sum(occ_m(3 - ispin, 1:m_max))
                        sum_occ_same_excl = sum(occ_m(ispin, 1:m_max)) - occ_m(ispin, m1)
                        if (m1 == m2) then
                           dUdn = common_pref*sum_u_aux/max(d1, eps_den) - &
                                  (2.0_rp*ubar/max(d1, eps_den))*(sum_occ_opposite - sum_occ_same_excl)
                           dJdn = common_pref*sum_j_aux/max(d2, eps_den) - &
                                  (2.0_rp*jbar/max(d2, eps_den))*sum_occ_same_excl
                        else
                           dUdn = common_pref*sum_u_aux/max(d1, eps_den)
                           dJdn = common_pref*sum_j_aux/max(d2, eps_den)
                        end if
                        dUeff_dn = dUdn - dJdn
                     end if
                     if (m1 == m2) then
                        if (use_acbn0) then
                           hub_pot(na, l_index, ispin, m1, m2) = 0.0_rp
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) + &
                              ueff*(0.5_rp - ldm(na, l_index, ispin, m1, m2)) + 0.5_rp*dUeff_dn*tr_n1mn
                        else
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                              - hub_u(na, l_index)*(n_tot(na, l_index) - 0.5_rp) &
                              + hub_j(na, l_index)*(n_spin(na, l_index, ispin) - 0.5_rp)
                        end if
                     else if (use_acbn0) then
                        hub_pot(na, l_index, ispin, m1, m2) = -ueff*ldm(na, l_index, ispin, m1, m2) + 0.5_rp*dUeff_dn*tr_n1mn
                     end if
                  end do
               end do
            end do

            ! Diagnostic for expected splitting: average diagonal on-site +U potential per spin.
            vdiag_up_avg = 0.0_rp
            vdiag_dn_avg = 0.0_rp
            do m1 = 1, m_max
               vdiag_up_avg = vdiag_up_avg + hub_pot(na, l_index, 1, m1, m1)
               vdiag_dn_avg = vdiag_dn_avg + hub_pot(na, l_index, 2, m1, m1)
            end do
            vdiag_up_avg = vdiag_up_avg/real(max(1, m_max), rp)
            vdiag_dn_avg = vdiag_dn_avg/real(max(1, m_max), rp)
            vdiag_split = vdiag_up_avg - vdiag_dn_avg
            if (rank == 0) then
               call g_logger%info('HUBBARD_VDIAG type='//fmt('i4', na)//' l='//fmt('i2', l_index - 1)// &
                  ' avg_up='//fmt('f11.6', vdiag_up_avg)//' avg_dn='//fmt('f11.6', vdiag_dn_avg)// &
                  ' split(up-dn)='//fmt('f11.6', vdiag_split), __FILE__, __LINE__)
            end if
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 0, min(3, this%control%lmax)
            do i = 1, 2*l + 1
               do j = 1, 2*l + 1
                  this%hubbard_u_pot(l**2 + i, l**2 + j, na) = hub_pot(na, l + 1, 1, i, j)
                  this%hubbard_u_pot(l**2 + i + spin_off, l**2 + j + spin_off, na) = hub_pot(na, l + 1, 2, i, j)
               end do
            end do
         end do
      end do

      do na = 1, this%lattice%ntype
         hup_max_abs = 0.0_rp
         do i = 1, nb
            do j = 1, nb
               hup_max_abs = max(hup_max_abs, abs(this%hubbard_u_pot(i, j, na)))
            end do
         end do
         if (rank == 0) call g_logger%info('HUBBARD_POT type='//fmt('i4', na)//' maxabs='//fmt('es12.4', hup_max_abs), __FILE__, __LINE__)
      end do

      do na = 1, this%lattice%ntype
         if (allocated(l_arr(na)%val)) deallocate(l_arr(na)%val)
      end do
   end subroutine calculate_hubbard_u_potential_general

   module subroutine calculate_hubbard_v_potential(this)
      class(hamiltonian), intent(inout) :: this

      integer :: itype, ia, nr, m, ja, jt, ispin, li, lj, ii, i0, i1, jdim
      real(rp) :: rmin, rcur, tol, occ_up, occ_dn, nn_shell_tol
      real(rp), dimension(3) :: dr
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin
      logical, save :: warned_proxy_once = .false.

      this%hubbard_v_pot(:, :, :, :) = 0.0_rp
      n_spin(:, :, :) = 0.0_rp
      tol = 1.0e-4_rp

      if (rank == 0 .and. .not. warned_proxy_once) then
         call g_logger%warning('HUBBARD +V currently uses a diagonal local-occupation proxy. Full inter-site n^{JI}_{mm''} from inter-site Green functions is not yet wired into calculate_hubbard_v_potential.', __FILE__, __LINE__)
         warned_proxy_once = .true.
      end if

      ! Spin- and l-resolved occupations from local density matrices.
      do itype = 1, this%lattice%ntype
         do li = 1, min(4, this%control%lmax + 1)
            jdim = 2*li - 1
            do ispin = 1, 2
               do ii = 1, jdim
                  n_spin(itype, li, ispin) = n_spin(itype, li, ispin) + &
                     this%lattice%symbolic_atoms(itype)%potential%ldm(li, ispin, ii, ii)
               end do
            end do
         end do
      end do

      do itype = 1, this%lattice%ntype
         ia = this%lattice%atlist(itype)
         nr = this%lattice%nn(ia, 1)

         ! Find nearest-neighbour shell radius for this atom type.
         rmin = huge(1.0_rp)
         do m = 2, nr
            ja = this%lattice%nn(ia, m)
            if (ja == 0) cycle
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk, this%lattice%cr*this%lattice%alat, ia, ja, dr)
            else
               dr(:) = (this%lattice%cr(:, ja) - this%lattice%cr(:, ia))*this%lattice%alat
            end if
            rcur = norm2(dr)
            if (rcur > tol .and. rcur < rmin) rmin = rcur
         end do
         if (rmin >= huge(1.0_rp)*0.5_rp) cycle

         do m = 2, nr
            ja = this%lattice%nn(ia, m)
            if (ja == 0) cycle
            jt = this%lattice%iz(ja)
            if (jt <= 0 .or. jt > this%lattice%ntype) cycle

            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk, this%lattice%cr*this%lattice%alat, ia, ja, dr)
            else
               dr(:) = (this%lattice%cr(:, ja) - this%lattice%cr(:, ia))*this%lattice%alat
            end if
            rcur = norm2(dr)
            nn_shell_tol = max(5.0e-4_rp, 2.5e-3_rp*rmin)
            if (abs(rcur - rmin) > nn_shell_tol) cycle

            ! Approximate intersite +V contribution (current proxy implementation):
            ! V^I,sigma_mm' += - sum_{J in NN} V^{IJ}_{li,lj} * n^{J,sigma}_{lj} / (2*li-1) * delta_mm'
            !
            ! IMPORTANT LIMITATION:
            ! The full +V expression requires inter-site density matrices n^{JI,sigma}_{mprime,m}
            ! from inter-site Green functions. That machinery is not yet wired here, so only
            ! an orbital-diagonal proxy based on local occupations is applied.
            do li = 1, min(4, this%control%lmax + 1)
               if (li > size(this%hubbard_v, 3)) cycle
               if (li > size(this%hubbard_v, 4)) cycle
               jdim = 2*li - 1
               i0 = (li - 1)*(li - 1) + 1
               i1 = i0 + jdim - 1
               do lj = 1, min(4, size(this%hubbard_v, 4), this%control%lmax + 1)
                  if (lj > size(this%hubbard_v, 3)) cycle
                  if (abs(this%hubbard_v(itype, jt, li, lj)) <= 1.0e-12_rp) cycle
                  occ_up = n_spin(jt, lj, 1)
                  occ_dn = n_spin(jt, lj, 2)
                  do ii = i0, i1
                     this%hubbard_v_pot(ii, ii, m, itype) = this%hubbard_v_pot(ii, ii, m, itype) - &
                        this%hubbard_v(itype, jt, li, lj)*occ_up/max(1, jdim)
                     this%hubbard_v_pot(ii + spin_off, ii + spin_off, m, itype) = &
                        this%hubbard_v_pot(ii + spin_off, ii + spin_off, m, itype) - &
                        this%hubbard_v(itype, jt, li, lj)*occ_dn/max(1, jdim)
                  end do
               end do
            end do
         end do
      end do
   end subroutine calculate_hubbard_v_potential

end submodule hamiltonian_hubbard
