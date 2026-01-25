###############################################################################
# Toolchain via modules (Apolo-3 / OHPC)
# Assumes you did something like:
#   module use /opt/ohpc/pub/moduledeps/gnu14-openmpi5
#   module load gnu14-openmpi5
#   module load cuda/12.5
#   # optional: module load netcdf-fortran/4.6.2 hdf5 ...
###############################################################################

# --- MPI wrappers (these pick correct GCC/OpenMPI from modules) ---
CC      = mpicc
CXX     = mpic++
f77     = mpifort
f90     = mpifort
CLINKER = mpic++

# --- CUDA (from module) ---
CUDA_HOME ?= $(CUDA_HOME)
Cu        = $(CUDA_HOME)/bin/nvcc

# --- Common optimization flags ---
CXXAPPFLAGS  = -O3 -Wno-deprecated -Dfortran3 -Dnewc
f90appflags  = -O3 -x f95-cpp-input
CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc

###############################################################################
# Include paths
# 1) MPI include dirs from wrapper
# 2) CUDA include from CUDA_HOME
# 3) Any other modules that add CPATH (netcdf, hdf5, etc.) are included too
###############################################################################

# Convert "dir1 dir2 ..." -> "-Idir1 -Idir2 ..."
mkI = $(patsubst %,-I%,$(1))

MPI_INCDIRS  := $(shell $(CXX) --showme:incdirs 2>/dev/null)
CPATH_DIRS   := $(subst :, ,$(CPATH))

filein = \
  $(call mkI,$(MPI_INCDIRS)) \
  -I$(CUDA_HOME)/include \
  $(call mkI,$(CPATH_DIRS))

###############################################################################
# Library paths + libs
# 1) MPI lib dirs/libs from wrapper
# 2) CUDA runtime lib
# 3) Any other modules that add LIBRARY_PATH are included too
###############################################################################

# Convert "dir1 dir2 ..." -> "-Ldir1 -Ldir2 ..."
mkL = $(patsubst %,-L%,$(1))

MPI_LIBDIRS := $(shell $(CXX) --showme:libdirs 2>/dev/null)
MPI_LIBS_RAW := $(shell $(CXX) --showme:libs 2>/dev/null)

# If a token is:
#   - already a flag (-lfoo, -Wl,...), keep it
#   - an absolute path (/something), keep it
#   - otherwise, treat it as a library name and prefix -l
MPI_LIBS := $(foreach w,$(MPI_LIBS_RAW),\
  $(if $(filter -%,$(w)),$(w),\
    $(if $(filter /%,$(w)),$(w),-l$(w))\
  )\
)

LIBPATH_DIRS := $(subst :, ,$(LIBRARY_PATH))

LDFLAGS += \
  $(call mkL,$(MPI_LIBDIRS)) \
  -L$(CUDA_HOME)/lib64 -L$(CUDA_HOME)/targets/x86_64-linux/lib \
  $(call mkL,$(LIBPATH_DIRS))

LDLIBS += \
  $(MPI_LIBS) \
  -lcudart \
  -lgfortran -lquadmath

