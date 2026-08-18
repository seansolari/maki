FROM rockylinux:9

ENV PATH="/opt/venv/bin:${PATH}"

RUN dnf update -y && \
    dnf install -y \
    python3 \
    python3-devel \
    python3-pip \
    cmake \
    gcc \
    gcc-c++ \
    git \
    tbb-devel

WORKDIR /opt/maki

COPY . /opt/maki

RUN python3 -m venv /opt/venv
RUN . /opt/venv/bin/activate

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install .

RUN python3 -m pip install nuitka[onefile]

RUN python3 -m nuitka \
    --standalone \
    --onefile \
    --lto=yes \
    --output-filename=maki-ubuntu24 \
    python/maki/cli.py
