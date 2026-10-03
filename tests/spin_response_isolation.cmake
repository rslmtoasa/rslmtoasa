# Fails if any source/spin_response* file `use`s a module of the archived
# linear-response code (tag lr-campaign-archive-2026-10).
# Run as: cmake -DSOURCE_DIR=<repo root> -P spin_response_isolation.cmake

set(_archived "(linear_response|lr_|lmto_path_operator)[a-z0-9_]*|(lmto_pair_potential|lmto_magnetic_tangent|radial_ground_state|pauli_ground_state_projection|exchange_q)(_mod)?")
set(_use_regex "(^|\n)[ \t]*use[ \t]*(,[^:\n]*::)?[ \t]*(${_archived})([^a-z0-9_]|$)")

file(GLOB _files "${SOURCE_DIR}/source/spin_response*")
set(_violations "")
foreach(_file IN LISTS _files)
  file(READ "${_file}" _text)
  string(TOLOWER "${_text}" _text)
  string(REGEX MATCHALL "${_use_regex}" _hits "${_text}")
  foreach(_hit IN LISTS _hits)
    string(STRIP "${_hit}" _hit)
    list(APPEND _violations "${_file}: ${_hit}")
  endforeach()
endforeach()

list(LENGTH _files _nfiles)
if(_violations)
  string(REPLACE ";" "\n  " _report "${_violations}")
  message(FATAL_ERROR "spin_response* uses archived modules:\n  ${_report}")
endif()
message(STATUS "checked ${_nfiles} spin_response* files: no archived `use`")
