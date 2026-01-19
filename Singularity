Bootstrap: docker
From: ubuntu:latest

%setup

    git archive --format=tar.gz --output=maki-HEAD.tar.gz --prefix maki-latest/ HEAD

%files

    maki-HEAD.tar.gz /opt

%labels

    Maintainer ssolari

    Version v2

%environment
    
    export BINDIR="/opt/maki/bin"

%post

    export APPDIR="/opt/maki/tmp"

    export BINDIR="/opt/maki/bin"


    apt-get update

    apt-get install -y cmake build-essential git libatomic1 zlib1g-dev libhdf5-dev libunwind-dev 

    rm -rf /var/lib/apt/lists/*


    mkdir -p $APPDIR

    tar -C $APPDIR -xzvf /opt/maki-HEAD.tar.gz
    
    rm /opt/maki-HEAD.tar.gz

    
    cd $APPDIR/maki-latest

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

    cmake --build build --config Release --target maki_build_targets -j $(nproc --all) --


    mkdir -p $BINDIR

    cp $APPDIR/maki-latest/build/bin/* $BINDIR


    apt-get purge --auto-remove -y cmake git

%runscript

    exec ${BINDIR}/maki "$@"
