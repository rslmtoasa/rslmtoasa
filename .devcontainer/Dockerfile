FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential cmake gfortran libopenblas-dev liblapack-dev \
        python3 python3-pip python3-venv python3-numpy python3-matplotlib \
        gnuplot-nox fonts-dejavu-core \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
