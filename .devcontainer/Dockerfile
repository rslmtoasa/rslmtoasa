FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential cmake gfortran libopenblas-dev liblapack-dev libsymspg-dev pkg-config python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
