# maki (pre-alpha)

Note unit tests for core functionality are still being developed.

## Build from source

These steps ensure an up-to-date version of `cmake`, and `gxx` compiler (at least g++ 10 is required to support C++ Concepts).

```{console}
$ mamba create -n [NAME] -c anaconda -c conda-forge cmake gcc gxx zlib tbb tbb-devel pybind11
$ mamba activate [NAME]
```

Installing with the `setuptools` first builds the core C++ library and bindings with `cmake`, and then installs the Python application (& CLI).

```{console}
$ python -m pip install .
```

## Command Line Interface

For a list of available commands, run `maki --help`.

