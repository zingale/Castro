PBOXLIB_HOME= ../../..
FBOXLIB_HOME= ../../../../fParallel

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 1

COMP	   = Intel
FCOMP	   = Intel

USE_MPI    = TRUE

USE_GRAV   = FALSE
USE_HYPRE  = TRUE
USE_RAD    = TRUE
USE_NEUTRINO = FALSE

# has the f90 routine ca_initgroups and group structure in probdata_module
HAS_INITGROUPS = TRUE

EOS_dir     = $(FBOXLIB_HOME)/extern/EOS/gamma_law
Network_dir = $(FBOXLIB_HOME)/extern/networks/gammalaw

# use radiation interpbndry
RAD_INTERP = TRUE

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
include ../Make.hypre

