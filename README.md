# maki

## Build using `cmake`

### Prepare mamba environment

These steps ensure an up-to-date version of `cmake`, and `gxx` compiler (at least g++ 10 is required to support C++ Concepts).

```{console}
$ mamba create -n cmake -c anaconda -c conda-forge cmake gcc gxx zlib pybind11
$ mamba activate cmake
```

### Build using cmake

```{console}
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$ cmake --build build --config Release --target maki-build-targets -j $(nproc --all) --
```

#### Note - conda-packaged compilers (VSCode)

To use a compiler installer within a conda environment (through the VSCode interface), add the conda `bin` path (via `conda env list`)  to `Cmake: Additional Compiler Search Dirs` settings.

