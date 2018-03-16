# Installing `Armadillo` on ubuntu

The precompiled version of `Armadillo` ready for download from Ubuntu repository is quite old, therefore, I downloaded the source code, compiled, and install my own version of the library on our Ubuntu machine.

The process was fairly straightforward, with a few obstacles to mention:

* I had to install LAPACK, ARPACK, SuperLU, OpenBLAS myself.
* CMake configure process, kept finding HDF5 and MKL libraries that came with `Anaconda` that where located in the Anaconda folder. Therefore I did the following things:
  * Used a flag in CMake that was something like -D-No_HDF5 to disable use of HDF5 library.
  * Also, added the path to the MKL library in Anaconda folder using `-Wl,-rpath=/home/amirhossein/anaconda3/lib/` flag in the linking process in the makefile, so that later, it could be found by Armadillo.

Therefore, Armadillo used MKL library instead of LAPACK or OpenBLAS.