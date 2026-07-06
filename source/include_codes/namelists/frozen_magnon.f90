! frozen_magnon: sweeps &hamiltonian's q_ss over a list of points and evaluates
! the energy at each. The first q-point is always the reference point (full SCF
! is run there unconditionally; every other point's E(q)/omega(q) is measured
! relative to it).
!
! mode = 'mft' (default): converge SCF fully at q_ss_list(:,1), then reuse
!   that converged potential for a single-iteration band-energy pass at every
!   other q (magnetic force theorem).
! mode = 'scf': converge SCF fully, independently, at every q.
!
! Prefer q_file for long paths. Its format is:
!   n_q_points
!   qx qy qz
!   ...
! q_coordinates = 'cartesian' (default): q rows are Cartesian units of 2*pi/alat,
! matching &hamiltonian q_ss. q_coordinates = 'direct': q rows are reciprocal
! lattice coordinates and are converted after the lattice is built.
! If q_file is blank, q_ss_list has a fixed compile-time column bound (500)
! because Fortran namelists cannot size an allocatable array from the input
! file alone in a way this codebase already exercises elsewhere; n_q_points
! selects how many of those columns are actually used.
integer, parameter :: frozen_magnon_max_q = 500
character(len=10) :: mode
character(len=sl) :: output_file
character(len=sl) :: q_file
character(len=16) :: q_coordinates
integer :: n_q_points
real(rp), dimension(3, frozen_magnon_max_q) :: q_ss_list

namelist /frozen_magnon/ mode, output_file, q_file, q_coordinates, n_q_points, q_ss_list
