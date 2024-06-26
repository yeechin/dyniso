
# Makefile Flags for Linux using ICC compiler 

LOCAL = /opt/local

COMPATH	=
CXX	= g++-12
CC	= gcc-mp-12
FC	= gfortran 

override ADDONS := $(ADDONS)

override SPECIAL := $(SPECIAL) ranf

# dependency check flag

DEP = -M

#==============================================================================
# Third party libraries
#==============================================================================

ifdef USE_OPENMP
OMP = -omp
endif

ifdef USE_FLOAT
FLOAT = -s
endif

# FFTW 

#FFTDIR  = $(LOCAL)/fftw$(OMP)$(FLOAT)
FFTDIR = $(LOCAL)
FFTINC  = -I$(FFTDIR)/include
ifdef USE_OPENMP
FFTLIB  = -L$(FFTDIR)/lib -ldrfftw_threads -ldfftw_threads -ldrfftw -ldfftw
else
FFTLIB  = -L$(FFTDIR)/lib -lrfftw -lfftw
endif

# optimization settings

ifeq ($(USE_OPENMP),1)
  OPENMP = -fopenmp
endif

ifdef USE_FLOAT
  DOUBLE =
else
  DOUBLE = -fdefault-integer-8 -fdefault-real-8 -DR8
endif

OPTXX  = -O2 -m64 $(OPENMP)
MOPTXX = -O3 -m64 $(OPENMP)

OPT  = $(OPTXX)
MOPT = $(MOPTXX)

OPTF  = -m64 -O2 -std=legacy -w $(OPENMP) $(DOUBLE) -DRANF -DFFTW $(MCPU)
MOPTF = -m64 -O3 -w $(OPENMP) $(DOUBLE) -DRANF -DFFTW \
        -ftree-vectorize -mtune=core2 -mfpmath=sse -msse -msse2 -msse3 -msse4 \
        -ftree-vectorizer-verbose=1 -ffast-math $(MCPU)

DEBUG  = -g $(DOUBLE) $(OPENMP)
DEBUGF = -g $(DOUBLE) $(OPENMP)

# loader and loader settings

LD	= $(FC)
DEBUGLD	= $(DEBUGF)
MOPTLD	= $(MOPTF)
OPTLD	= $(OPTF) 

INCS = $(FFTINC)
LIBS = $(FFTLIB)

# X stuff

XLIBDIRS = -L/usr/X11/lib
XINCDIRS = -I/usr/include/X11
