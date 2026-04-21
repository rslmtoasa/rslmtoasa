add_library(strux_lib STATIC
    strux_errors.f90
    strux_utils.f90
    strux_angular.f90
    strux_bessel.f90
    strux_geometry.f90
    strux_linalg.f90
    strux_pair_expansion.f90
    strux_solver_assembly.f90
    strux_support.f90
    strux_lmto47_mode.f90
    strux_api.f90
    strux_tb.f90
)

target_link_libraries(strux_lib PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
