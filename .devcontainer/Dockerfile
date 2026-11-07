FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential cmake gfortran libopenblas-dev liblapack-dev libsymspg-dev pkg-config \
        python3 python3-pip python3-venv python3-numpy python3-matplotlib python3-f90nml \
        gnuplot-nox fonts-dejavu-core \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
