# librsrec -- GPU/CPU recursion backends for RS-LMTO-ASA
#
# Targets:
#   make cpu        librsrec_cpu.so   (OpenMP + LAPACK reference, validated)
#   make gpu        librsrec_gpu.so   (CUDA + cuBLAS, host LAPACK for zheev)
#   make test       run Python validation against the NumPy Fortran mirror
#   make ftest      Fortran wrapper end-to-end test (links the CPU backend)
#
# Both libraries export the same C API (rsrec.h); link whichever one you
# want against recursion_gpu_mod.f90 -- no code changes needed to switch.

CXX      ?= g++
NVCC     ?= nvcc
FC       ?= gfortran
CXXFLAGS ?= -O3 -fopenmp -fPIC -march=native
NVFLAGS  ?= -O3 -arch=native -Xcompiler -fPIC
LAPACK   ?= -llapack -lblas
# In sandboxes without the -dev symlinks, point LAPACK at the versioned
# libraries instead, e.g.:
# LAPACK = /lib/x86_64-linux-gnu/liblapack.so.3 /lib/x86_64-linux-gnu/libblas.so.3

cpu: librsrec_cpu.so
gpu: librsrec_gpu.so

librsrec_cpu.so: rsrec_cpu.cpp rsrec.h
	$(CXX) $(CXXFLAGS) -shared $< -o $@ $(LAPACK)

librsrec_gpu.so: rsrec_gpu.cu rsrec.h
	$(NVCC) $(NVFLAGS) -shared $< -o $@ -lcublas -lcufft $(LAPACK)

recursion_gpu_mod.o: recursion_gpu_mod.f90
	$(FC) -O2 -c $<

test: librsrec_cpu.so
	python3 test_rsrec.py

ftest: librsrec_cpu.so recursion_gpu_mod.o
	$(FC) -O2 test_recursion_gpu.f90 recursion_gpu_mod.o -o test_fortran \
	      -L. -lrsrec_cpu -Wl,-rpath,'$$ORIGIN'
	./test_fortran

clean:
	rm -f *.so *.o *.mod test_fortran

.PHONY: cpu gpu test ftest clean
