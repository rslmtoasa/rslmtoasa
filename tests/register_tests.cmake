# Test registration. Included from the top-level CMakeLists.txt.
#
# Every test is registered through rslmto_add_test and carries exactly one tier label:
#   quick       each test <= 10 s, whole tier <= 2 min (-j4); guards core behaviour; runs on every push
#   standard    whole tier <= 20 min (-j4); runs nightly
#   validation  long studies; run manually
# A test registered without a valid tier stops the configure step.

include(CTest)
set(TestBinary ${CMAKE_BINARY_DIR}/bin/${PROJECT_BINARY})
set(TestDIR ${CMAKE_SOURCE_DIR}/tests)

# rslmto_add_test(NAME <name> TIER <quick|standard|validation>
#                 COMMAND <cmd...> [LABELS <extra...>] [TIMEOUT <s>] [WORKING_DIRECTORY <dir>]
#                 [PROPERTIES <prop value ...>])
# Default timeouts: quick 30 s, standard 600 s, validation none.
function(rslmto_add_test)
  cmake_parse_arguments(T "" "NAME;TIER;TIMEOUT;WORKING_DIRECTORY" "COMMAND;LABELS;PROPERTIES" ${ARGN})
  if(NOT T_NAME)
    message(FATAL_ERROR "rslmto_add_test: NAME is required")
  endif()
  if(NOT T_TIER)
    message(FATAL_ERROR "rslmto_add_test: test '${T_NAME}' is registered without a TIER (quick|standard|validation)")
  endif()
  if(NOT T_TIER MATCHES "^(quick|standard|validation)$")
    message(FATAL_ERROR "rslmto_add_test: test '${T_NAME}' has invalid TIER '${T_TIER}' (quick|standard|validation)")
  endif()
  if(NOT T_COMMAND)
    message(FATAL_ERROR "rslmto_add_test: test '${T_NAME}' has no COMMAND")
  endif()
  set(_labels "${T_TIER}")
  list(APPEND _labels ${T_LABELS})
  list(REMOVE_DUPLICATES _labels)
  add_test(NAME "${T_NAME}" COMMAND ${T_COMMAND})
  set_tests_properties("${T_NAME}" PROPERTIES LABELS "${_labels}")
  if(T_TIMEOUT)
    set_tests_properties("${T_NAME}" PROPERTIES TIMEOUT "${T_TIMEOUT}")
  elseif(T_TIER STREQUAL "quick")
    set_tests_properties("${T_NAME}" PROPERTIES TIMEOUT 30)
  elseif(T_TIER STREQUAL "standard")
    set_tests_properties("${T_NAME}" PROPERTIES TIMEOUT 600)
  endif()
  if(T_WORKING_DIRECTORY)
    set_tests_properties("${T_NAME}" PROPERTIES WORKING_DIRECTORY "${T_WORKING_DIRECTORY}")
  endif()
  if(T_PROPERTIES)
    set_tests_properties("${T_NAME}" PROPERTIES ${T_PROPERTIES})
  endif()
endfunction()

# Optional semantic labels of one manifest case.
function(read_manifest_case_labels MANIFEST_JSON CASE_INDEX OUT_VAR)
  set(_labels "")
  string(JSON _label_count ERROR_VARIABLE _label_error LENGTH "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "labels")
  if(NOT _label_error AND _label_count GREATER 0)
    math(EXPR _label_last "${_label_count} - 1")
    foreach(_label_index RANGE ${_label_last})
      string(JSON _label GET "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "labels" ${_label_index})
      list(APPEND _labels "${_label}")
    endforeach()
  endif()
  set(${OUT_VAR} "${_labels}" PARENT_SCOPE)
endfunction()

# Effective structure-constant backend of one manifest case (case value, else manifest default, else legacy).
function(read_manifest_case_strux_backend MANIFEST_JSON CASE_INDEX OUT_VAR)
  set(_backend "legacy")
  string(JSON _default_backend ERROR_VARIABLE _default_error GET "${MANIFEST_JSON}" "defaults" "namelists" "lattice" "strux_backend")
  if(NOT _default_error)
    set(_backend "${_default_backend}")
  endif()
  string(JSON _case_backend ERROR_VARIABLE _case_error GET "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "namelists" "lattice" "strux_backend")
  if(NOT _case_error)
    set(_backend "${_case_backend}")
  endif()
  set(${OUT_VAR} "${_backend}" PARENT_SCOPE)
endfunction()

# Tier and timeout of one manifest case; a case without a tier stops the configure step.
function(read_manifest_case_tier MANIFEST_JSON CASE_INDEX OUT_TIER OUT_TIMEOUT)
  string(JSON _tier ERROR_VARIABLE _tier_error GET "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "tier")
  if(_tier_error)
    string(JSON _name GET "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "name")
    message(FATAL_ERROR "manifest case '${_name}' has no \"tier\"")
  endif()
  string(JSON _timeout ERROR_VARIABLE _timeout_error GET "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "timeout")
  if(_timeout_error)
    set(_timeout "")
  endif()
  set(${OUT_TIER} "${_tier}" PARENT_SCOPE)
  set(${OUT_TIMEOUT} "${_timeout}" PARENT_SCOPE)
endfunction()

# Cases of the group "ldau_guards" are run together by tests/ldau/test_guards.py.
function(manifest_case_is_grouped MANIFEST_JSON CASE_INDEX OUT_VAR)
  string(JSON _group ERROR_VARIABLE _group_error GET "${MANIFEST_JSON}" "cases" ${CASE_INDEX} "group")
  if(_group_error)
    set(${OUT_VAR} FALSE PARENT_SCOPE)
  else()
    set(${OUT_VAR} TRUE PARENT_SCOPE)
  endif()
endfunction()

set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(EXAMPLE_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE FILEPATH
    "Python interpreter for example tests (must have f90nml installed)")
set(PY "${EXAMPLE_PYTHON_EXECUTABLE}")

# Every path that a registered test reads must exist (tests/test_inventory.csv).
rslmto_add_test(NAME CheckTestData TIER quick LABELS tooling
  COMMAND "${PY}" "${TestDIR}/check_test_data.py")

##############################################################################################
# Regression (backend variants) and route-agnostic estimator triads
##############################################################################################
if(RUN_REG_TESTS)
  set(REGRESSION_CASES_JSON "${TestDIR}/regression/cases.json")
  file(READ "${REGRESSION_CASES_JSON}" _regression_json)
  string(JSON _regression_n LENGTH "${_regression_json}" "cases")
  math(EXPR _regression_last "${_regression_n} - 1")
  foreach(_i RANGE ${_regression_last})
    string(JSON _name GET "${_regression_json}" "cases" ${_i} "name")
    read_manifest_case_tier("${_regression_json}" ${_i} _tier _timeout)
    read_manifest_case_labels("${_regression_json}" ${_i} _case_labels)
    read_manifest_case_strux_backend("${_regression_json}" ${_i} _strux_backend)
    set(_labels "functional;regression;backend")
    list(APPEND _labels ${_case_labels})
    if(_strux_backend STREQUAL "strux_lib")
      list(APPEND _labels "strux")
    elseif(_strux_backend STREQUAL "legacy")
      list(APPEND _labels "legacy_strux")
    endif()
    rslmto_add_test(NAME "Regression_${_name}" TIER ${_tier} TIMEOUT "${_timeout}" LABELS ${_labels}
      COMMAND "${PY}" "${TestDIR}/regression/run_matrix.py" "--binary" "${TestBinary}"
        "--cases-json" "${REGRESSION_CASES_JSON}" "--case-name" "${_name}"
        "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/regression" "--references" "${TestDIR}/regression/references"
      PROPERTIES SKIP_RETURN_CODE 77)
  endforeach()

  # Each triad runs one observable three ways (recursion, lehmann, dyson) and checks the agreement envelopes.
  set(TRIAD_CASES_JSON "${TestDIR}/regression/triad_cases.json")
  file(READ "${TRIAD_CASES_JSON}" _triad_json)
  string(JSON _triad_n LENGTH "${_triad_json}" "cases")
  math(EXPR _triad_last "${_triad_n} - 1")
  foreach(_i RANGE ${_triad_last})
    string(JSON _name GET "${_triad_json}" "cases" ${_i} "name")
    read_manifest_case_tier("${_triad_json}" ${_i} _tier _timeout)
    rslmto_add_test(NAME "Triad_${_name}" TIER ${_tier} TIMEOUT "${_timeout}" LABELS functional regression triad
      COMMAND "${PY}" "${TestDIR}/regression/run_triad.py" "--binary" "${TestBinary}"
        "--cases-json" "${TRIAD_CASES_JSON}" "--case-name" "${_name}"
        "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/triad" "--references" "${TestDIR}/regression/references_triad")
  endforeach()

  # Validation studies (manual tier).
  foreach(_v IN ITEMS
      "Val05GreenConvergence|val05_green_convergence.py|val05_green_convergence|green;convergence"
      "Val07ExchangeFormulations|val07_exchange_formulations.py|val07_exchange_formulations|exchange;magnetic"
      "Val08DampingInertia|val08_damping_inertia.py|val08_damping_inertia|damping;inertia"
      "Val09KuboBastinTransport|val09_kubo_bastin_transport.py|val09_kubo_bastin_transport|conductivity;transport"
      "Val12LmtoFieldsTorques|val12_lmto_fields_torques.py|val12_lmto_fields_torques|torque;constraints;soc"
      "Val13AbInitioSpinDynamics|val13_abinitio_spin_dynamics.py|val13_abinitio_spin_dynamics|spin_dynamics"
      "Val15MultilayerVacuumElectrostatics|val15_multilayer_vacuum.py|val15_multilayer_vacuum|electrostatics;vacuum")
    string(REPLACE "|" ";" _f "${_v}")
    list(GET _f 0 _vname)
    list(GET _f 1 _vscript)
    list(GET _f 2 _vscratch)
    list(GET _f 3 _vlabels)
    rslmto_add_test(NAME "${_vname}" TIER validation LABELS functional ${_vlabels}
      COMMAND "${PY}" "${TestDIR}/validation/${_vscript}" "--binary" "${TestBinary}"
        "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/validation/${_vscratch}")
  endforeach()

  # Per-atom forced nonmagnetic SCF (Fe magnetic, Co constrained; default flags equal all-false flags).
  foreach(_nsp IN ITEMS 1 4)
    rslmto_add_test(NAME Nonmagnetic_SCF_nsp${_nsp} TIER standard TIMEOUT 720 LABELS functional force_nonmagnetic
      COMMAND "${PY}" "${TestDIR}/nonmagnetic/run_scf.py" "--binary" "${TestBinary}"
        "--workdir" "${CMAKE_BINARY_DIR}/Testing/nonmagnetic/nsp${_nsp}" "--nsp" "${_nsp}")
  endforeach()
  if(ENABLE_MPI)
    rslmto_add_test(NAME Nonmagnetic_SCF_MPI TIER standard TIMEOUT 720 LABELS functional force_nonmagnetic mpi
      COMMAND "${PY}" "${TestDIR}/nonmagnetic/run_scf.py" "--binary" "${TestBinary}"
        "--workdir" "${CMAKE_BINARY_DIR}/Testing/nonmagnetic/mpi" "--nsp" "4"
        "--mpi-exec" "${MPIEXEC_EXECUTABLE}" "--ranks" "2")
  endif()

  # --- Kubo transport operators and backends (pytest) -----------------------------------------
  # The compile-based suites build Debug dense-reference drivers themselves (about 75 s for both);
  # test_transport_backend.py runs the executable of this build (RSLMTO_EXE).
  set(_pytest_env "${CMAKE_COMMAND}" -E env "FC=${CMAKE_Fortran_COMPILER}" "OPENBLAS_NUM_THREADS=1"
      "RS_LMTO_TEST_CMAKE_ARGS=-DENABLE_MARCH_NATIVE=OFF")
  rslmto_add_test(NAME TransportOperatorSource TIER quick LABELS transport regression
    COMMAND ${_pytest_env} "${PY}" -m pytest -q -p no:cacheprovider
      "${TestDIR}/regression/test_octupole_operators.py" "${TestDIR}/regression/test_quadrupole_operators.py")
  rslmto_add_test(NAME TransportRegression TIER standard TIMEOUT 900 LABELS transport regression
    COMMAND ${_pytest_env} "${PY}" -m pytest -q -p no:cacheprovider
      "${TestDIR}/regression/test_hoh_transport.py" "${TestDIR}/regression/test_response_units.py")
  rslmto_add_test(NAME TransportBackend TIER standard TIMEOUT 900 LABELS transport regression
    COMMAND ${_pytest_env} "RSLMTO_EXE=${TestBinary}" "${PY}" -m pytest -q -p no:cacheprovider
      "${TestDIR}/regression/test_transport_backend.py")
endif()

##############################################################################################
# Standalone Fortran unit tests
##############################################################################################
if(RUN_UNIT_TESTS)
  set(_standalone_unit_targets "")

  function(_rslmto_link_unit_target _name)
    if(ENABLE_OPENMP AND ENABLE_MPI AND OpenMP_Fortran_FOUND AND MPI_Fortran_FOUND)
      target_link_libraries(${_name} ${COMMON_LIBRARIES} MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
    elseif(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
      target_link_libraries(${_name} ${COMMON_LIBRARIES} OpenMP::OpenMP_Fortran)
    elseif(ENABLE_MPI AND MPI_Fortran_FOUND)
      target_link_libraries(${_name} ${COMMON_LIBRARIES} MPI::MPI_Fortran)
    else()
      target_link_libraries(${_name} ${COMMON_LIBRARIES})
    endif()
    set_target_properties(${_name} PROPERTIES LINKER_LANGUAGE "Fortran")
    list(APPEND _standalone_unit_targets ${_name})
    set(_standalone_unit_targets "${_standalone_unit_targets}" PARENT_SCOPE)
  endfunction()

  # add_fortran_unit_test(<name> <tier> <source...> [ARGS <args...>] [LABELS ...] [PROPERTIES ...] [WORKING_DIRECTORY d])
  function(add_fortran_unit_test _name _tier)
    cmake_parse_arguments(U "" "WORKING_DIRECTORY" "SOURCES;ARGS;LABELS;PROPERTIES" ${ARGN})
    add_executable(${_name} ${U_SOURCES})
    _rslmto_link_unit_target(${_name})
    set(_standalone_unit_targets "${_standalone_unit_targets}" PARENT_SCOPE)
    set(_wd "")
    if(U_WORKING_DIRECTORY)
      set(_wd WORKING_DIRECTORY "${U_WORKING_DIRECTORY}")
    endif()
    rslmto_add_test(NAME ${_name} TIER ${_tier} LABELS unit ${U_LABELS} ${_wd}
      COMMAND ${_name} ${U_ARGS} PROPERTIES ${U_PROPERTIES})
  endfunction()

  # Driver executables of consolidated tests share the sub-case runner.
  add_library(SubcaseRunner STATIC ${TestDIR}/unit/subcase_runner.f90)
  target_link_libraries(SubcaseRunner ${COMMON_LIBRARIES})

  # --- core units --------------------------------------------------------------------------
  add_fortran_unit_test(UnitBasisDimensions quick SOURCES ${TestDIR}/unit/test_basis_dimensions.f90)
  foreach(_b1 IN ITEMS
      "UnitScfB1GuardCapacityRequiresFatal|capacity|B1 occupation-capacity fatal.*site=3.*channel=d.*iteration=17.*q0="
      "UnitScfB1GuardWarningBoundary|warning|B1 occupation-capacity warning.*site=3.*channel=d.*iteration=17.*q0="
      "UnitScfB1GuardNonfiniteRequiresFatal|nonfinite|B1 non-finite guard.*site=2.*channel=p.*iteration=9.*field=new q1.*value="
      "UnitScfB1GuardBoundary|boundary|")
    string(REPLACE "|" ";" _f "${_b1}")
    list(GET _f 0 _bname)
    list(GET _f 1 _barg)
    list(GET _f 2 _bregex)
    if(_bregex)
      add_fortran_unit_test(${_bname} quick SOURCES ${TestDIR}/unit/test_scf_b1_guard.f90 ARGS ${_barg}
        PROPERTIES PASS_REGULAR_EXPRESSION "${_bregex}")
    else()
      add_fortran_unit_test(${_bname} quick SOURCES ${TestDIR}/unit/test_scf_b1_guard.f90 ARGS ${_barg})
    endif()
  endforeach()
  add_fortran_unit_test(UnitGreenLifecycle quick SOURCES ${TestDIR}/unit/test_green_lifecycle.f90 LABELS green)
  add_fortran_unit_test(UnitNonmagneticProjection quick SOURCES ${TestDIR}/unit/test_nonmagnetic_projection.f90 LABELS force_nonmagnetic)
  add_fortran_unit_test(UnitControlSpinSemantics quick SOURCES ${TestDIR}/unit/test_control_spin_semantics.f90 LABELS control spin)
  add_fortran_unit_test(UnitControlConstraintsRead quick SOURCES ${TestDIR}/unit/test_control_constraints_read.f90 LABELS control constraints)
  add_fortran_unit_test(UnitLehmannChain quick SOURCES ${TestDIR}/unit/test_lehmann_chain.f90)
  add_fortran_unit_test(UnitSimpsonFTemperature quick SOURCES ${TestDIR}/unit/test_simpson_f_temperature.f90)
  add_fortran_unit_test(UnitKspaceGFValidation quick SOURCES ${TestDIR}/unit/test_kspace_gf_validation.f90 LABELS green kspace)
  add_fortran_unit_test(UnitFermiDiracEntropy quick SOURCES ${TestDIR}/unit/test_fermi_dirac_entropy.f90 LABELS thermodynamics)
  add_fortran_unit_test(UnitChebyshevGFOracle quick SOURCES ${TestDIR}/unit/test_chebyshev_gf_oracle.f90 LABELS green chebyshev)
  add_fortran_unit_test(UnitChebyshevGreenPrecision quick SOURCES ${TestDIR}/unit/test_chebyshev_green_precision.f90 LABELS green chebyshev)
  add_fortran_unit_test(UnitChebyshevSpectralPositivity quick SOURCES ${TestDIR}/unit/test_chebyshev_spectral_positivity.f90 LABELS green chebyshev)
  add_fortran_unit_test(UnitDysonEquivalence quick SOURCES ${TestDIR}/unit/test_dyson_equivalence.f90 LABELS green)
  add_fortran_unit_test(UnitGammaSupercell quick SOURCES ${TestDIR}/unit/test_gamma_supercell.f90)
  add_fortran_unit_test(UnitBsfSumRule quick SOURCES ${TestDIR}/unit/test_bsf_sumrule.f90)
  add_fortran_unit_test(UnitMomentKernel quick SOURCES ${TestDIR}/unit/test_moment_kernel.f90)
  add_fortran_unit_test(UnitKpmTransport quick SOURCES ${TestDIR}/unit/test_kpm_transport.f90 LABELS transport kpm)
  add_fortran_unit_test(UnitMadl2dExhGuard quick SOURCES ${TestDIR}/unit/test_madl2d_exh_guard.f90 LABELS electrostatics)
  # This test also inspects the madl2d source, so that the corrected overflow guard cannot silently regress.
  target_compile_definitions(UnitMadl2dExhGuard PRIVATE
    MADL2D_SOURCE="${CMAKE_CURRENT_SOURCE_DIR}/source/charge_madelung_2d.f90")
  add_fortran_unit_test(UnitMadelungConventions quick SOURCES ${TestDIR}/unit/test_madelung_conventions.f90 LABELS electrostatics)
  add_fortran_unit_test(UnitRegionRegistry quick SOURCES ${TestDIR}/unit/test_region_registry.f90 LABELS electrostatics)
  add_fortran_unit_test(UnitInterfaceAlignmentOracle quick SOURCES ${TestDIR}/unit/test_interface_alignment_oracle.f90 LABELS electrostatics)
  add_fortran_unit_test(UnitVacuumLead quick SOURCES ${TestDIR}/unit/test_vacuum_lead.f90 LABELS electrostatics)
  add_fortran_unit_test(UnitVacuumRegionWiring quick SOURCES ${TestDIR}/unit/test_vacuum_region_wiring.f90 LABELS electrostatics)
  add_fortran_unit_test(UnitAlignmentSolver quick SOURCES ${TestDIR}/unit/test_alignment_solver.f90 LABELS electrostatics)
  add_fortran_unit_test(UnitParallelContext quick SOURCES ${TestDIR}/unit/test_parallel_context.f90 LABELS mpi parallel_context)
  rslmto_add_test(NAME UnitParallelContext_cuda_override TIER quick LABELS unit mpi parallel_context
    COMMAND $<TARGET_FILE:UnitParallelContext> PROPERTIES ENVIRONMENT "RSLMTO_CUDA_DEVICE=1")
  add_fortran_unit_test(UnitKpointWorkset quick SOURCES ${TestDIR}/unit/test_kpoint_workset.f90 LABELS reciprocal)
  add_fortran_unit_test(UnitKpointWorksetReducedReject quick SOURCES ${TestDIR}/unit/test_kpoint_workset_reduced_rejection.f90
    LABELS reciprocal PROPERTIES WILL_FAIL TRUE)
  add_fortran_unit_test(UnitArbitraryKEigenpairs quick SOURCES ${TestDIR}/unit/test_arbitrary_k_eigenpairs.f90 LABELS reciprocal)
  add_fortran_unit_test(UnitStructureConstantsBackends quick SOURCES ${TestDIR}/unit/test_structure_constants_backends.f90
    LABELS strux WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/structure_constants")
  target_compile_definitions(UnitStructureConstantsBackends PRIVATE
    TEST_STRUCTURE_CONSTANT_DATABASE="${CMAKE_SOURCE_DIR}/tests/scf/cases/bulk/bccFe")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/structure_constants")
  add_fortran_unit_test(UnitAbspinlibIntegrators quick SOURCES ${TestDIR}/validation/test_abspinlib_integrators.f90 LABELS spin_integrators)
  add_fortran_unit_test(UnitConstrainingField quick SOURCES ${TestDIR}/unit/test_constraining_field.f90 LABELS constraints)

  # --- XC: three consolidated executables ---------------------------------------------------
  add_fortran_unit_test(UnitXcKernels quick LABELS xc
    SOURCES ${TestDIR}/unit/test_xc_kernels.f90
      ${TestDIR}/unit/test_xc_selector_semantics.f90 ${TestDIR}/unit/test_xc_lda_legacy_kernel.f90
      ${TestDIR}/unit/test_xc_legacy_unpolarized_density.f90 ${TestDIR}/unit/test_radial_gga.f90
      ${TestDIR}/unit/test_legacy_pbe_gga.f90 ${TestDIR}/unit/test_vxc0sp_gga_origin.f90
      ${TestDIR}/unit/test_atomsc_residual_controls.f90 ${TestDIR}/unit/test_xc_legacy_polarization_tolerance.f90)
  target_link_libraries(UnitXcKernels SubcaseRunner)
  add_fortran_unit_test(UnitXcGuards quick LABELS xc negative
    SOURCES ${TestDIR}/unit/test_xc_guards.f90
      ${TestDIR}/unit/test_xc_unknown_selector_requires_fatal.f90 ${TestDIR}/unit/test_xc_legacy_polarized_requires_fatal.f90
      ${TestDIR}/unit/test_xc_legacy_polarization_tolerance.f90)
  target_link_libraries(UnitXcGuards SubcaseRunner)
  if(ENABLE_LIBXC)
    add_fortran_unit_test(UnitXcLibxc quick LABELS xc libxc
      SOURCES ${TestDIR}/unit/test_xc_libxc.f90
        ${TestDIR}/unit/test_libxc_xc_baseline.f90 ${TestDIR}/unit/test_libxc_gga_radial.f90
        ${TestDIR}/unit/test_libxc_production_contract.f90 ${TestDIR}/unit/test_xc_composition_semantics.f90
        ${TestDIR}/unit/test_pbe_gga_comparison.f90 ${TestDIR}/unit/test_xc_lda_reconciliation.f90
        ${TestDIR}/unit/test_xc_unsupported_family_requires_fatal.f90 ${TestDIR}/unit/test_xc_capability_requires_fatal.f90
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/xc_lda_reconciliation")
    target_link_libraries(UnitXcLibxc SubcaseRunner)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/xc_lda_reconciliation")
  else()
    add_fortran_unit_test(UnitXcSelectorRequiresLibxc quick LABELS xc negative
      SOURCES ${TestDIR}/unit/test_xc_selector_requires_libxc.f90 PROPERTIES WILL_FAIL TRUE)
  endif()

  # --- LDA+U -------------------------------------------------------------------------------
  add_fortran_unit_test(UnitLdaUTensor quick SOURCES ${TestDIR}/unit/test_lda_u_tensor.f90 LABELS lda_u)
  add_fortran_unit_test(UnitHubbardSpinor quick SOURCES ${TestDIR}/unit/test_hubbard_spinor.f90 LABELS lda_u
    ARGS "${TestDIR}/ldau/data/spinor_oracle_reference.txt")
  add_fortran_unit_test(UnitLdaULocalMoments quick LABELS lda_u
    SOURCES ${TestDIR}/unit/test_ldau_local_moments.f90
      ${TestDIR}/unit/test_local_moments_kernel.f90 ${TestDIR}/unit/test_hubbard_spinor_moments.f90
      ${TestDIR}/unit/test_local_moments_site_swap.f90
    ARGS all "${TestDIR}/ldau/data/a3_m14_reference_fixture.txt" "${TestDIR}/ldau/data/a3_m14_fresh_site_swap_reference.txt")
  target_link_libraries(UnitLdaULocalMoments SubcaseRunner)
  add_fortran_unit_test(UnitLdaUProjectorFrozenH quick LABELS lda_u projector
    SOURCES ${TestDIR}/unit/test_ldau_projector_frozen_h.f90
      ${TestDIR}/unit/test_ldau_projector_common_frozen_h.f90 ${TestDIR}/unit/test_ldau_projector_multisite_frozen_h.f90)
  target_link_libraries(UnitLdaUProjectorFrozenH SubcaseRunner)
  rslmto_add_test(NAME LdaUGuards TIER standard LABELS lda_u guards
    COMMAND "${PY}" "${TestDIR}/ldau/test_guards.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/guards")
  rslmto_add_test(NAME LdaUHohFermiShift TIER standard LABELS lda_u hoh
    COMMAND "${PY}" "${TestDIR}/ldau/test_hoh_u.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/hoh_u")
  rslmto_add_test(NAME LdaUNsp34 TIER standard LABELS lda_u nsp34
    COMMAND "${PY}" "${TestDIR}/ldau/test_nsp34_u.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/nsp34")
  rslmto_add_test(NAME LdaUNsp34Noncollinear TIER standard LABELS lda_u nsp34 TIMEOUT 1200
    COMMAND "${PY}" "${TestDIR}/ldau/test_nsp34_noncollinear.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/nsp34_noncollinear")
  rslmto_add_test(NAME LdaUGbtHoh TIER standard LABELS lda_u hoh gbt TIMEOUT 1800
    COMMAND "${PY}" "${TestDIR}/ldau/test_gbt_hoh_u.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/gbt_hoh")
  rslmto_add_test(NAME LdaUKspaceGbt TIER standard LABELS lda_u gbt
    COMMAND "${PY}" "${TestDIR}/ldau/test_kspace_gbt_u.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/kspace_gbt")
  rslmto_add_test(NAME LdaUKspaceLocalNsp3 TIER standard LABELS lda_u nsp34
    COMMAND "${PY}" "${TestDIR}/ldau/test_kspace_local_nsp3.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/kspace_local_nsp3")
  rslmto_add_test(NAME LdaUImpurityHoh TIER standard LABELS lda_u hoh TIMEOUT 1200
    COMMAND "${PY}" "${TestDIR}/ldau/test_impurity_hoh_u.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/impurity_hoh")
  rslmto_add_test(NAME LdaULocalNoncollinear TIER standard LABELS lda_u nsp34 TIMEOUT 1200
    COMMAND "${PY}" "${TestDIR}/ldau/test_local_noncollinear_u.py" "--binary" "${TestBinary}"
      "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/local_noncollinear")

  # Stored-reference rows (tests/ldau/reference_rows.py; tests/ldau/references/*.json).
  # LdaURefRsImpurityNsp2 is tier validation: on Linux/OpenBLAS the nsp=2 (SOC = 0) run differs from nsp=1 by 1e-6
  # over 4 iterations against the tolerance of 1e-7 set on the developer platform (KNOWN_ISSUES.md).
  foreach(_row IN ITEMS
      "LdaURefRsBulkNsp1U|U_rs|standard"
      "LdaURefKsBulkNsp2SocUTilted|ks_soc_u2_x|standard"
      "LdaURefRsImpurityNsp2|imp_rs_soc0|validation"
      "LdaURefSurfaceNsp2Soc0|surf_soc0|standard"
      "LdaURefNoUKeyOutputs|noU_rs|standard")
    string(REPLACE "|" ";" _f "${_row}")
    list(GET _f 0 _rname)
    list(GET _f 1 _rrow)
    list(GET _f 2 _rtier)
    rslmto_add_test(NAME ${_rname} TIER ${_rtier} LABELS lda_u regression
      COMMAND "${PY}" "${TestDIR}/ldau/reference_rows.py" ${_rrow} "--binary" "${TestBinary}"
        "--scratch-root" "${CMAKE_BINARY_DIR}/Testing/ldau/${_rrow}")
  endforeach()

  # --- spin response -----------------------------------------------------------------------
  # The oracle .nml files are regenerated at test time by tests/spin_response/oracles/make_references*.py
  # into a scratch copy of the oracle directory (tests/spin_response/with_oracles.py).
  set(_sr_oracles "${CMAKE_BINARY_DIR}/Testing/spin_response/oracles")
  add_library(SpinResponseTolerances STATIC ${TestDIR}/unit/spin_response_tolerances.f90)
  target_link_libraries(SpinResponseTolerances ${COMMON_LIBRARIES})
  add_library(SpinResponseFeState STATIC ${TestDIR}/unit/spin_response_fe_state.f90)
  target_link_libraries(SpinResponseFeState ${COMMON_LIBRARIES} SpinResponseTolerances)

  add_executable(UnitSpinResponseKernels ${TestDIR}/unit/test_spin_response_kernel_cases.f90
    ${TestDIR}/unit/test_spin_response_kernels.f90 ${TestDIR}/unit/test_spin_response_channels_kernel.f90
    ${TestDIR}/unit/test_spin_response_a4.f90 ${TestDIR}/unit/test_spin_response_a5.f90 ${TestDIR}/unit/test_spin_response_r1.f90)
  _rslmto_link_unit_target(UnitSpinResponseKernels)
  target_link_libraries(UnitSpinResponseKernels SpinResponseTolerances SubcaseRunner)
  rslmto_add_test(NAME UnitSpinResponseKernels TIER quick LABELS unit spin-response
    COMMAND "${PY}" "${TestDIR}/spin_response/with_oracles.py" "${_sr_oracles}_kernels" --
      $<TARGET_FILE:UnitSpinResponseKernels> all "{oracles}")

  # bcc Fe scratch run, then the unit checks on its output (and on the AF scratch for s2af).
  rslmto_add_test(NAME SpinResponseFeDriver TIER standard LABELS spin-response
    COMMAND "${PY}" "${TestDIR}/spin_response/fe_scratch.py" "${TestBinary}" "${TestDIR}/spin_response/fe_bcc"
      "${CMAKE_BINARY_DIR}/Testing/spin_response/fe_bcc"
    PROPERTIES FIXTURES_SETUP SpinResponseScratch RESOURCE_LOCK SpinResponseScratch)
  rslmto_add_test(NAME SpinResponseAfDriver TIER standard LABELS spin-response
    COMMAND "${PY}" "${TestDIR}/spin_response/af_scratch.py" "${TestBinary}" "${TestDIR}/spin_response/fe_bcc"
      "${CMAKE_BINARY_DIR}/Testing/spin_response/af"
    PROPERTIES FIXTURES_SETUP SpinResponseScratch RESOURCE_LOCK SpinResponseScratch)
  add_executable(UnitSpinResponseDriven ${TestDIR}/unit/test_spin_response_driver_cases.f90
    ${TestDIR}/unit/test_spin_response_c1.f90 ${TestDIR}/unit/test_spin_response_c2.f90
    ${TestDIR}/unit/test_spin_response_c3.f90 ${TestDIR}/unit/test_spin_response_s2.f90
    ${TestDIR}/unit/test_spin_response_channels.f90 ${TestDIR}/unit/test_spin_response_radial.f90
    ${TestDIR}/unit/test_spin_response_s2af.f90)
  _rslmto_link_unit_target(UnitSpinResponseDriven)
  target_link_libraries(UnitSpinResponseDriven SpinResponseFeState SubcaseRunner)
  rslmto_add_test(NAME UnitSpinResponseDriven TIER standard LABELS unit spin-response
    COMMAND "${PY}" "${TestDIR}/spin_response/with_oracles.py" "${_sr_oracles}_driven" --
      $<TARGET_FILE:UnitSpinResponseDriven> all "{oracles}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/spin_response/fe_bcc"
    PROPERTIES FIXTURES_REQUIRED SpinResponseScratch RESOURCE_LOCK SpinResponseScratch)
  set_tests_properties(UnitSpinResponseDriven PROPERTIES DEPENDS "SpinResponseFeDriver;SpinResponseAfDriver")

  rslmto_add_test(NAME SpinResponseSocZeroIdentity TIER standard LABELS spin-response
    COMMAND "${PY}" "${TestDIR}/spin_response/soc_zero_identity.py" "${TestBinary}" "${TestDIR}/spin_response/fe_bcc"
      "${CMAKE_BINARY_DIR}/Testing/spin_response/soc0")
  rslmto_add_test(NAME SpinResponseInversionOdd TIER standard LABELS spin-response
    COMMAND "${PY}" "${TestDIR}/spin_response/inversion_odd.py" "${TestBinary}" "${TestDIR}/spin_response/fe_bcc"
      "${CMAKE_BINARY_DIR}/Testing/spin_response/inversion_odd" "${TestDIR}/spin_response/oracles/tolerances.nml")
  rslmto_add_test(NAME SpinResponseBaselineRepro TIER standard LABELS spin-response regression
    COMMAND "${PY}" "${TestDIR}/spin_response/baseline_repro.py" "${TestBinary}" "${TestDIR}/spin_response/fe_bcc"
      "${CMAKE_BINARY_DIR}/Testing/spin_response/baseline" "${TestDIR}/spin_response/baseline"
      "${TestDIR}/spin_response/oracles/tolerances.nml")
  rslmto_add_test(NAME SpinResponseIsolation TIER quick LABELS spin-response tooling
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR} -P ${TestDIR}/spin_response_isolation.cmake)

  # --- GBT ---------------------------------------------------------------------------------
  add_fortran_unit_test(UnitGbtKspaceForceTheorem quick SOURCES ${TestDIR}/unit/test_gbt_kspace_force_theorem.f90
    ARGS ${TestDIR}/gbt/oracles LABELS gbt)
  add_fortran_unit_test(UnitGbtKspaceHorth quick SOURCES ${TestDIR}/unit/test_gbt_kspace_horth.f90
    ARGS ${TestDIR}/gbt/oracles LABELS gbt)
  add_fortran_unit_test(UnitGbtWp02GaugeShiftedK quick SOURCES ${TestDIR}/unit/test_gbt_wp02_gauge_shifted_k.f90 LABELS gbt)

  # --- tooling -----------------------------------------------------------------------------
  rslmto_add_test(NAME UnitManifestDefaults TIER quick LABELS tooling manifest
    COMMAND ${Python3_EXECUTABLE} ${TestDIR}/unit/test_manifest_defaults.py)

  if(ENABLE_CUDA_PLUGIN)
    add_fortran_unit_test(CudaMpiDeviceMapping standard SOURCES ${TestDIR}/unit/test_cuda_mpi_device_mapping.f90
      LABELS cuda mpi PROPERTIES SKIP_RETURN_CODE 77)
    target_link_libraries(CudaMpiDeviceMapping rsrec_cuda)
  endif()
  if(ENABLE_CUDA_RECIPROCAL)
    add_executable(ReciprocalCudaArbitraryK ${TestDIR}/unit/test_reciprocal_cuda_arbitrary_k.f90)
    target_link_libraries(ReciprocalCudaArbitraryK PRIVATE ${COMMON_LIBRARIES} rslmto_reciprocal_cuda)
    set_target_properties(ReciprocalCudaArbitraryK PROPERTIES LINKER_LANGUAGE "Fortran")
    rslmto_add_test(NAME ReciprocalCudaArbitraryK TIER standard LABELS unit cuda reciprocal
      COMMAND ReciprocalCudaArbitraryK PROPERTIES SKIP_RETURN_CODE 77)
  endif()

  if(ENABLE_MPI AND MPI_Fortran_FOUND)
    # Homebrew Open MPI on macOS may report one slot; the parallel-context unit uses four ranks.
    set(_rslmto_mpi_unit_preflags ${MPIEXEC_PREFLAGS})
    if(APPLE)
      list(APPEND _rslmto_mpi_unit_preflags --map-by :OVERSUBSCRIBE)
    endif()
    foreach(_parallel_context_ranks IN ITEMS 1 2 4)
      rslmto_add_test(NAME "UnitParallelContext_mpi_${_parallel_context_ranks}" TIER quick LABELS unit mpi parallel_context
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_parallel_context_ranks} ${_rslmto_mpi_unit_preflags}
                $<TARGET_FILE:UnitParallelContext> ${MPIEXEC_POSTFLAGS}
        PROPERTIES PROCESSORS ${_parallel_context_ranks})
    endforeach()
  endif()

  add_custom_target(unit_tests DEPENDS ${_standalone_unit_targets})
endif()

##############################################################################################
# Example (application) tests
##############################################################################################
if(RUN_EXAMPLE_TESTS)
  # Register one CTest entry per case of a cases.json manifest.
  # Args: CASES_JSON  SUITE_LABEL  SCRATCH_ROOT  REFS_DIR
  function(register_example_tests CASES_JSON SUITE SCRATCH_ROOT REFS_DIR)
    file(READ "${CASES_JSON}" _json)
    string(JSON _n LENGTH "${_json}" "cases")
    math(EXPR _last "${_n} - 1")
    foreach(_i RANGE ${_last})
      string(JSON _name GET "${_json}" "cases" ${_i} "name")
      string(JSON _recur GET "${_json}" "cases" ${_i} "namelists" "control" "recur")
      manifest_case_is_grouped("${_json}" ${_i} _grouped)
      if(_grouped)
        continue()
      endif()
      read_manifest_case_tier("${_json}" ${_i} _tier _timeout)
      # Cases requiring an optional CMake feature are skipped when that feature is off.
      string(JSON _req_opt ERROR_VARIABLE _req_opt_err GET "${_json}" "cases" ${_i} "requires_cmake_option")
      if(NOT _req_opt_err AND NOT ${_req_opt})
        continue()
      endif()
      read_manifest_case_labels("${_json}" ${_i} _case_labels)
      read_manifest_case_strux_backend("${_json}" ${_i} _case_strux_backend)
      if(_case_strux_backend STREQUAL "strux_lib")
        list(APPEND _case_labels "strux")
      elseif(_case_strux_backend STREQUAL "legacy")
        list(APPEND _case_labels "legacy_strux")
      endif()
      set(_cmd "${PY}" "${TestDIR}/run_test.py" "--binary" "${TestBinary}" "--cases-json" "${CASES_JSON}"
        "--case-name" "${_name}" "--scratch-root" "${SCRATCH_ROOT}" "--compare-ref" "${REFS_DIR}"
        "--abs-tol" "${EXAMPLE_REF_ABS_TOL}" "--rel-tol" "${EXAMPLE_REF_REL_TOL}")
      # Cases listing "launch_modes" register a serial/MPI pair that shares one reference.
      string(JSON _n_modes ERROR_VARIABLE _modes_err LENGTH "${_json}" "cases" ${_i} "launch_modes")
      if(_modes_err)
        set(_default_cmd ${_cmd})
        if(ENABLE_MPI)
          list(APPEND _default_cmd "--mpi-enabled")
        endif()
        rslmto_add_test(NAME "${_name}" TIER ${_tier} TIMEOUT "${_timeout}" LABELS functional example ${SUITE} ${_recur} ${_case_labels}
          COMMAND ${_default_cmd})
      else()
        string(JSON _mpi_procs GET "${_json}" "cases" ${_i} "mpi_procs")
        rslmto_add_test(NAME "${_name}" TIER ${_tier} TIMEOUT "${_timeout}"
          LABELS functional example ${SUITE} ${_recur} launch_serial ${_case_labels}
          COMMAND ${_cmd} "--force-serial")
        if(ENABLE_MPI)
          string(JSON _n_mpi_ranks ERROR_VARIABLE _mpi_ranks_err LENGTH "${_json}" "cases" ${_i} "mpi_rank_counts")
          if(_mpi_ranks_err)
            set(_n_mpi_ranks 1)
            set(_mpi_rank_0 "${_mpi_procs}")
          endif()
          math(EXPR _last_mpi_rank "${_n_mpi_ranks} - 1")
          foreach(_mpi_rank_index RANGE ${_last_mpi_rank})
            if(NOT _mpi_ranks_err)
              string(JSON _mpi_rank GET "${_json}" "cases" ${_i} "mpi_rank_counts" ${_mpi_rank_index})
            else()
              set(_mpi_rank "${_mpi_rank_0}")
            endif()
            if(_mpi_rank EQUAL _mpi_procs)
              set(_mpi_suffix "_mpi")
            else()
              set(_mpi_suffix "_mpi_${_mpi_rank}")
            endif()
            set(_mpi_test_name "${_name}${_mpi_suffix}")
            rslmto_add_test(NAME "${_mpi_test_name}" TIER ${_tier} TIMEOUT "${_timeout}"
              LABELS functional example ${SUITE} ${_recur} launch_mpi ${_case_labels}
              COMMAND ${_cmd} "--mpi-enabled" "--mpi-procs" "${_mpi_rank}" "--scratch-name" "${_mpi_test_name}"
              PROPERTIES PROCESSORS "${_mpi_rank}")
          endforeach()
        endif()
      endif()
    endforeach()
  endfunction()
  register_example_tests("${TestDIR}/scf/cases.json" scf "${CMAKE_BINARY_DIR}/Testing/scf" "${TestDIR}/scf/references")
  register_example_tests("${TestDIR}/postproc/cases.json" postproc "${CMAKE_BINARY_DIR}/Testing/postproc" "${TestDIR}/postproc/references")

  # --- GBT (k-space) application tests ---------------------------------------------------------
  rslmto_add_test(NAME AppGbtKspaceAntiparallel TIER standard LABELS functional gbt kspace
    COMMAND "${PY}" "${TestDIR}/gbt/test_kspace_antiparallel.py" "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/gbt/antiparallel")
  rslmto_add_test(NAME AppGbtKspaceCommensurate TIER standard LABELS functional gbt kspace
    COMMAND "${PY}" "${TestDIR}/gbt/test_kspace_commensurate.py" "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/gbt/commensurate" first)
  rslmto_add_test(NAME AppGbtKspaceSymmetry TIER standard LABELS functional gbt kspace
    COMMAND "${PY}" "${TestDIR}/gbt/test_kspace_symmetry.py" "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/gbt/symmetry")
  rslmto_add_test(NAME AppGbtKspaceScSpiral_converge TIER standard LABELS functional gbt kspace
    COMMAND "${PY}" "${TestDIR}/gbt/test_kspace_scspiral.py" converge "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/gbt/scspiral_converge")
  rslmto_add_test(NAME AppGbtKspaceScSpiral_guards TIER quick LABELS functional gbt kspace
    COMMAND "${PY}" "${TestDIR}/gbt/test_kspace_scspiral.py" guards "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/gbt/scspiral_guards")

  # --- other application tests -----------------------------------------------------------------
  rslmto_add_test(NAME AppConstraintsNamelistMom TIER quick LABELS functional constraints
    COMMAND "${PY}" "${TestDIR}/constraints/test_constraints_namelist.py" mom "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/constraints/mom")
  rslmto_add_test(NAME AppConstraintsNamelistMalformed TIER quick LABELS functional constraints
    COMMAND "${PY}" "${TestDIR}/constraints/test_constraints_namelist.py" malformed "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/constraints/malformed")
  rslmto_add_test(NAME AppPeriodicNcRotation TIER standard LABELS functional noncollinear
    COMMAND "${PY}" "${TestDIR}/noncollinear/test_periodic_nc_rotation.py" "${TestBinary}" "${CMAKE_BINARY_DIR}/Testing/noncollinear/rotation")
  rslmto_add_test(NAME Example_si_chebyshev_kspace_dos_equivalence TIER quick
    LABELS functional example scf equivalence rs kspace chebyshev
    COMMAND "${PY}" "${TestDIR}/scf/test_si_dos_equivalence.py" --binary "${TestBinary}"
      --case-root "${TestDIR}/scf/cases/bulk/diamondSi" --scratch-root "${CMAKE_BINARY_DIR}/Testing/scf/si_dos_equivalence")
  if(ENABLE_LIBXC)
    rslmto_add_test(NAME Val03LibxcScfEquivalence TIER validation LABELS validation functional scf libxc equivalence
      COMMAND "${PY}" "${TestDIR}/validation/val03_libxc_scf.py" --binary "${TestBinary}"
        --scratch-root "${CMAKE_BINARY_DIR}/Testing/validation/val03_libxc_scf" --kspace-only)
  endif()
endif()
