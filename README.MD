# BernsteinFEM
An object-oriented approach to the computation of integrals for use in Finite Element Analysis using the Bernstein polynomials as basis function.

Based on the works and code by M. Ainsworth, G. Andriamaro and O. Davydov

There is a overview of the library's functionality [here](docs/Overview.md)

## Requirements

Armadillo library, found at http://arma.sourceforge.net

## Building

This library can be built for Linux, Windows or OSX.
* OSX building not tested yet.

### Linux:
  * You can build this library using the make command, it is recommended to use `make -j` to build multithreaded. This will make a BernsteinFEM.a static library which can be used when compiling your source code.
  * You can use the config file with the command `source config.sh` to set variables to use when compiling your source while in the working directory of this project.
  * Then build executables (e.g. from the source in Tests directory) with `g++ <your_source.cpp> $lib $flags $inc_dir`
  * You can also use CMake for building the library with your preferred IDE.
  
### Windows:
  * Install CMake 3.4.0 or higher;
  * Install Visual Studio;
  * Download and build the Armadillo library above;
  * Make sure you correctly set the Armadillo package variables on CMake project, for it is not done automatically correct on Windows.
  * Also link the BLAS and LAPACK libraries that come within the Armadillo pack, in the CMake project.
  * Ready to build.
