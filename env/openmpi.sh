#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Source this script instead of executing it: source env/openmpi.sh" >&2
  exit 1
fi

: "${RSLMTO_OPENMPI_BIN_DIR:=/usr/bin}"
: "${RSLMTO_OPENMPI_LIB_DIR:=/usr/lib/x86_64-linux-gnu/openmpi/lib}"
: "${RSLMTO_MKL_LIB_DIR:=/opt/intel/oneapi/mkl/2026.1/lib}"

_rslmto_filter_path_var() {
  local var_name="$1"
  shift

  local current_value filtered_value entry pattern skip_entry
  current_value="${!var_name-}"
  filtered_value=""

  IFS=':' read -r -a _rslmto_entries <<< "$current_value"
  for entry in "${_rslmto_entries[@]}"; do
    [[ -z "$entry" ]] && continue
    skip_entry=0
    for pattern in "$@"; do
      if [[ "$entry" =~ $pattern ]]; then
        skip_entry=1
        break
      fi
    done
    [[ $skip_entry -eq 1 ]] && continue
    if [[ -z "$filtered_value" ]]; then
      filtered_value="$entry"
    else
      filtered_value+="::$entry"
      filtered_value="${filtered_value/::/:}"
    fi
  done

  printf -v "$var_name" '%s' "$filtered_value"
  export "$var_name"
}

_rslmto_prepend_path_var() {
  local var_name="$1"
  local new_entry="$2"
  local current_value

  [[ -z "$new_entry" ]] && return
  [[ ! -d "$new_entry" ]] && return

  current_value="${!var_name-}"
  case ":$current_value:" in
    *":$new_entry:"*) ;;
    *)
      if [[ -n "$current_value" ]]; then
        printf -v "$var_name" '%s' "$new_entry:$current_value"
      else
        printf -v "$var_name" '%s' "$new_entry"
      fi
      export "$var_name"
      ;;
  esac
}

unset I_MPI_ROOT MPI_ROOT FI_PROVIDER MPICH_ROOT

_rslmto_filter_path_var PATH \
  '/opt/intel/oneapi/.*/bin$' \
  '/opt/nvidia/hpc_sdk/.*/comm_libs/mpi/bin$'
_rslmto_filter_path_var LIBRARY_PATH \
  '/opt/intel/oneapi/.*/lib($|/)' \
  '/opt/intel/oneapi/.*/libfabric($|/)' \
  '/opt/nvidia/hpc_sdk/.*/comm_libs/mpi/.*/lib($|/)'
_rslmto_filter_path_var LD_LIBRARY_PATH \
  '/opt/intel/oneapi/.*/lib($|/)' \
  '/opt/intel/oneapi/.*/libfabric($|/)' \
  '/opt/nvidia/hpc_sdk/.*/comm_libs/mpi/.*/lib($|/)'
_rslmto_filter_path_var PKG_CONFIG_PATH \
  '/opt/intel/oneapi/.*/lib/pkgconfig$' \
  '/opt/nvidia/hpc_sdk/.*/comm_libs/mpi/.*/pkgconfig$'

_rslmto_prepend_path_var PATH /bin
_rslmto_prepend_path_var LIBRARY_PATH "$RSLMTO_MKL_LIB_DIR"
_rslmto_prepend_path_var LD_LIBRARY_PATH "$RSLMTO_MKL_LIB_DIR"
_rslmto_prepend_path_var PATH "$RSLMTO_OPENMPI_BIN_DIR"
_rslmto_prepend_path_var LIBRARY_PATH "$RSLMTO_OPENMPI_LIB_DIR"
_rslmto_prepend_path_var LD_LIBRARY_PATH "$RSLMTO_OPENMPI_LIB_DIR"

hash -r

export FC="${RSLMTO_FC:-/usr/bin/gfortran}"
export F77="$FC"
export F90="$FC"
export MPI_Fortran_COMPILER="${RSLMTO_MPI_FORTRAN_COMPILER:-/usr/bin/mpifort}"
export MPIEXEC_EXECUTABLE="${RSLMTO_MPIEXEC_EXECUTABLE:-/usr/bin/mpirun}"

echo "Configured RS-LMTO Open MPI environment"
echo "  FC=$FC"
echo "  MPI_Fortran_COMPILER=$MPI_Fortran_COMPILER"
echo "  MPIEXEC_EXECUTABLE=$MPIEXEC_EXECUTABLE"
echo "  PATH head=${PATH%%:*}"
echo "  LIBRARY_PATH=${LIBRARY_PATH}"
echo "  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

unset -f _rslmto_filter_path_var
unset -f _rslmto_prepend_path_var
unset _rslmto_entries