/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/



#ifndef CONFIG_H
#define CONFIG_H

#define HAVEfftw3 1
/* #undef HAVEfftw */
#define HAVEeigen3EigenDense 1
/* #undef HAVEEigenDense */
#define FFTW_OMP
#define FFTW_TH
#if defined(FFTW_OMP) || defined(FFTW_TH)
#define FFTW3_THREADS
#endif
#define CMAKE_CXX_COMPILER "/usr/bin/g++-7"
#define CMAKE_CXX_FLAGS " -std=c++11 -fopenmp"
#define CMAKE_SYSTEM "Linux-4.15.0-1021-oem"
#define CMAKE_SYSTEM_PROCESSOR "x86_64"
/* #undef USE_MPI */
/* #undef HAVE_XDRFILE */
#ifdef HAVE_XDRFILE
#define XDR_FILE
#endif








#endif
