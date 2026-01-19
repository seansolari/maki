# STAGE ONE: BUILD
FROM ubuntu:latest AS builder

ENV MKISRCDIR=/opt/maki/tmp \
    APPDIR=/opt

RUN apt-get update && \
    apt-get install -y \
        cmake \
        build-essential \
        git \
        libatomic1 \
        zlib1g-dev \
        libhdf5-dev \
        libunwind-dev

RUN mkdir -p $MKISRCDIR

COPY maki-HEAD.tar.gz $MKISRCDIR

RUN mkdir -p $APPDIR && \
    tar -C $APPDIR -xvzf $MKISRCDIR/maki-HEAD.tar.gz && \
    rm -rf $MKISRCDIR

WORKDIR $APPDIR/maki-latest

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build --config Release --target maki_build_targets -j $(nproc --all) --

# STAGE TWO: RUN
# TODO: consider a smaller image in the future
FROM ubuntu:latest AS deployer

RUN apt-get update && \
    apt-get install -y \
        libatomic1 \
        zlib1g-dev \
        libhdf5-dev \
        libunwind-dev

RUN mkdir -p \
    /opt/maki-latest/build \
    /opt/maki-latest/data \
    /mnt/INPUT \
    /mnt/OUTPUT \
    /mnt/DATABASE \
    /mnt/LOGS \
    /mnt/TEMP

COPY --from=builder /opt/maki-latest/build /opt/maki-latest/build
COPY --from=builder /opt/maki-latest/data /opt/maki-latest/data

ENTRYPOINT ["/opt/maki-latest/build/bin/maki"]
