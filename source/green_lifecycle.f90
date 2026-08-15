submodule (green_mod) green_lifecycle
   implicit none

contains

   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(green)
   !---------------------------------------------------------------------------
   module function constructor(dos_obj) result(obj)
      type(green) :: obj
      type(dos), target, intent(in) :: dos_obj

      obj%dos => dos_obj
      obj%recursion => dos_obj%recursion
      obj%en => dos_obj%en
      obj%symbolic_atom => dos_obj%recursion%hamiltonian%charge%lattice%symbolic_atoms
      obj%lattice => dos_obj%recursion%lattice
      obj%control => dos_obj%recursion%lattice%control

      call obj%restore_to_default()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   module subroutine destructor(this)
      type(green) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%g0)) call g_safe_alloc%deallocate('green.g0', this%g0)
      if (allocated(this%gij)) call g_safe_alloc%deallocate('green.gij', this%gij)
      if (allocated(this%gji)) call g_safe_alloc%deallocate('green.gji', this%gji)
      if (allocated(this%ginmag)) call g_safe_alloc%deallocate('green.ginmag', this%ginmag)
      if (allocated(this%gjnmag)) call g_safe_alloc%deallocate('green.gjnmag', this%gjnmag)
      if (allocated(this%gix)) call g_safe_alloc%deallocate('green.gix', this%gix)
      if (allocated(this%giy)) call g_safe_alloc%deallocate('green.giy', this%giy)
      if (allocated(this%giz)) call g_safe_alloc%deallocate('green.giz', this%giz)
      if (allocated(this%gjx)) call g_safe_alloc%deallocate('green.gjx', this%gjx)
      if (allocated(this%gjy)) call g_safe_alloc%deallocate('green.gjy', this%gjy)
      if (allocated(this%gjz)) call g_safe_alloc%deallocate('green.gjz', this%gjz)
      if (allocated(this%g00ij)) call g_safe_alloc%deallocate('green.g00ij', this%g00ij)
      if (allocated(this%g01ij)) call g_safe_alloc%deallocate('green.g01ij', this%g01ij)
      if (allocated(this%g00ji)) call g_safe_alloc%deallocate('green.g00ji', this%g00ji)
      if (allocated(this%g01ji)) call g_safe_alloc%deallocate('green.g01ji', this%g01ji)
      if (allocated(this%gx0ij)) call g_safe_alloc%deallocate('green.gx0ij', this%gx0ij)
      if (allocated(this%gy0ij)) call g_safe_alloc%deallocate('green.gy0ij', this%gy0ij)
      if (allocated(this%gz0ij)) call g_safe_alloc%deallocate('green.gz0ij', this%gz0ij)
      if (allocated(this%gx1ij)) call g_safe_alloc%deallocate('green.gx1ij', this%gx1ij)
      if (allocated(this%gy1ij)) call g_safe_alloc%deallocate('green.gy1ij', this%gy1ij)
      if (allocated(this%gz1ij)) call g_safe_alloc%deallocate('green.gz1ij', this%gz1ij)
      if (allocated(this%gx0ji)) call g_safe_alloc%deallocate('green.gx0ji', this%gx0ji)
      if (allocated(this%gy0ji)) call g_safe_alloc%deallocate('green.gy0ji', this%gy0ji)
      if (allocated(this%gz0ji)) call g_safe_alloc%deallocate('green.gz0ji', this%gz0ji)
      if (allocated(this%gx1ji)) call g_safe_alloc%deallocate('green.gx1ji', this%gx1ji)
      if (allocated(this%gy1ji)) call g_safe_alloc%deallocate('green.gy1ji', this%gy1ji)
      if (allocated(this%gz1ji)) call g_safe_alloc%deallocate('green.gz1ji', this%gz1ji)
#else
      if (allocated(this%g0)) deallocate (this%g0)
      if (allocated(this%gij)) deallocate (this%gij)
      if (allocated(this%gji)) deallocate (this%gji)
      if (allocated(this%ginmag)) deallocate (this%ginmag)
      if (allocated(this%gjnmag)) deallocate (this%gjnmag)
      if (allocated(this%gix)) deallocate (this%gix)
      if (allocated(this%giz)) deallocate (this%giz)
      if (allocated(this%gjx)) deallocate (this%gjx)
      if (allocated(this%gjy)) deallocate (this%gjy)
      if (allocated(this%gjz)) deallocate (this%gjz)
      if (allocated(this%g00ij)) deallocate (this%g00ij)
      if (allocated(this%g01ij)) deallocate (this%g01ij)
      if (allocated(this%g00ji)) deallocate (this%g00ji)
      if (allocated(this%g01ji)) deallocate (this%g01ji)
      if (allocated(this%gx0ij)) deallocate (this%gx0ij)
      if (allocated(this%gy0ij)) deallocate (this%gy0ij)
      if (allocated(this%gz0ij)) deallocate (this%gz0ij)
      if (allocated(this%gx1ij)) deallocate (this%gx1ij)
      if (allocated(this%gy1ij)) deallocate (this%gy1ij)
      if (allocated(this%gz1ij)) deallocate (this%gz1ij)
      if (allocated(this%gx0ji)) deallocate (this%gx0ji)
      if (allocated(this%gy0ji)) deallocate (this%gy0ji)
      if (allocated(this%gz0ji)) deallocate (this%gz0ji)
      if (allocated(this%gx1ji)) deallocate (this%gx1ji)
      if (allocated(this%gy1ji)) deallocate (this%gy1ji)
      if (allocated(this%gz1ji)) deallocate (this%gz1ji)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   module subroutine restore_to_default(this)
      use mpi_mod
      class(green) :: this

#ifdef USE_SAFE_ALLOC
      if (this%lattice%njij == 0) then
         call g_safe_alloc%allocate('green.g0', this%g0, (/nb, nb, this%en%channels_ldos + 10, atoms_per_process/))
      else
         call g_safe_alloc%allocate('green.g0', this%g0, (/nb, nb, this%en%channels_ldos + 10, 4/))
      end if
      call g_safe_alloc%allocate('green.gij', this%gij, (/nb, nb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gji', this%gji, (/nb, nb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.ginmag', this%ginmag, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjnmag', this%gjnmag, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gix', this%gix, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.giy', this%giy, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.giz', this%giz, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjx', this%gjx, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjy', this%gjy, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjz', this%gjz, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g00ij', this%g00ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g00ji', this%g00ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g01ij', this%g01ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g01ji', this%g01ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx0ij', this%gx0ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy0ij', this%gy0ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz0ij', this%gz0ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx1ij', this%gx1ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy1ij', this%gy1ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz1ij', this%gz1ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx0ji', this%gx0ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy0ji', this%gy0ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz0ji', this%gz0ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx1ji', this%gx1ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy1ji', this%gy1ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz1ji', this%gz1ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gij_eta', this%gij_eta, (/64, nb, nb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gji_eta', this%gji_eta, (/64, nb, nb, atoms_per_process/))
      call g_safe_alloc%allocate('green.ginmag_eta', this%ginmag_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjnmag_eta', this%gjnmag_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gix_eta', this%gix_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.giy_eta', this%giy_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.giz_eta', this%giz_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjx_eta', this%gjx_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjy_eta', this%gjy_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjz_eta', this%gjz_eta, (/64, norb, norb, atoms_per_process/))
#else
      if (this%lattice%njij == 0) then
         allocate (this%g0(nb, nb, this%en%channels_ldos + 10, this%lattice%nrec))
      else
         allocate (this%g0(nb, nb, this%en%channels_ldos + 10, 4))
      end if
      allocate (this%gij(nb, nb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gji(nb, nb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%ginmag(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjnmag(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gix(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%giy(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%giz(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjx(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjy(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjz(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g00ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g00ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g01ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g01ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx0ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy0ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz0ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx1ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy1ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz1ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx0ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy0ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz0ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx1ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy1ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz1ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gij_eta(64, nb, nb, atoms_per_process))
      allocate (this%gji_eta(64, nb, nb, atoms_per_process))
      allocate (this%ginmag_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjnmag_eta(64, norb, norb, atoms_per_process))
      allocate (this%gix_eta(64, norb, norb, atoms_per_process))
      allocate (this%giy_eta(64, norb, norb, atoms_per_process))
      allocate (this%giz_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjx_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjy_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjz_eta(64, norb, norb, atoms_per_process))
#endif

      this%g0(:, :, :, :) = (0.0d0, 0.0d0)
      this%gij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gji(:, :, :, :) = (0.0d0, 0.0d0)
      this%ginmag(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjnmag(:, :, :, :) = (0.0d0, 0.0d0)
      this%gix(:, :, :, :) = (0.0d0, 0.0d0)
      this%giy(:, :, :, :) = (0.0d0, 0.0d0)
      this%giz(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjx(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjy(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjz(:, :, :, :) = (0.0d0, 0.0d0)
      this%g00ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%g00ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%g01ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%g01ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx0ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy0ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz0ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx1ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy1ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz1ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx0ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy0ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz0ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx1ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy1ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz1ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gij_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gji_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%ginmag_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjnmag_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gix_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%giy_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%giz_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjx_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjy_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjz_eta(:, :, :, :) = (0.0d0, 0.0d0)
   end subroutine restore_to_default

end submodule green_lifecycle
