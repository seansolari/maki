FROM ubuntu:24.04

ENV PATH="/opt/venv/bin:${PATH}"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        python3 \
        python3-venv \
        python3-dev \
        python3-pip \
        build-essential \
        cmake \
        gcc \
        g++ \
        git \
        libtbb-dev \
        patchelf \
        && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/maki

COPY . /opt/maki

RUN python3 -m venv /opt/venv
RUN . /opt/venv/bin/activate

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install .

RUN apt-get remove -y cmake git && \
    apt-get autoremove -y && \
    apt-get clean

RUN python3 -m pip install nuitka[onefile]

RUN python3 -m nuitka \
    --standalone \
    --onefile \
    --lto=yes \
    --output-filename=maki-ubuntu24 \
    python/maki/cli.py
