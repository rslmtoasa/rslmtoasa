# Add library
set(cmake_archive_output_directory lib)
message(STATUS "The project name is: ${PROJECT_NAME}")
add_library(${PROJECT_NAME} STATIC
    array.f90
    bands.f90
    calculation.f90
    charge.f90
    conductivity.f90
    control.f90
    density_of_states.f90
    element.f90
    energy.f90
    exchange.f90
    face.F90
    globals.f90
    green.f90
    hamiltonian.f90
    lattice.f90
    # lists.f90
    logger.f90
    math.f90
    mix.f90
    mpi.f90
    namelist_generator.f90
    os.f90
    potential.f90
    precision.f90
    reciprocal.f90
    recursion.f90
    report.f90
    safe_alloc.f90
    self.f90
    bounds.f90
    spglib_interface.f90
    spin_dynamics.f90
    string.f90
    symbolic_atom.f90
    symmetry.f90
    timer.f90
    xc.f90
)
add_subdirectory(./include_codes/abspinlib)
target_link_libraries(${PROJECT_NAME} abspinlib)
