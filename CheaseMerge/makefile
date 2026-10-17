# Makefile for MARS-F on lohan4 - H. Reimerdes
# 04-21-2006 First try
#

# Macros
FFLAGS = -OPT:Olimit=0 -m64

PROG =	 chease.x
BUILD_DIR ?= ../build

SRCS =   chease.f 

INCS =	       COMBND.inc  COMINT.inc  COMPLO.inc  DECLAR.inc	QUAQDQ.inc \
               COMCON.inc  COMIOD.inc  COMRCH.inc  HERMIT.inc	QUAQQQ.inc \
               BNDIND.inc  COMDAT.inc  COMISO.inc  COMSOL.inc           \
               COMDIM.inc  COMLAB.inc  COMSUR.inc  SOLOV.inc            \
               COMEQD.inc  COMMAP.inc  COMVAC.inc  NEWRUN.inc	        \
               COMERA.inc  COMNUM.inc  COMVEV.inc                       \
               COMBAL.inc  COMESH.inc  COMOPT.inc  CUCCCC.inc           \
               COMBLA.inc     COMETA.inc  COMPHY.inc  CUCDCD.inc  PROCESS.inc

MODS = 

.PRECIOUS:	$(SRCS) $(INCS) $(MODS)

OBJS =	$(SRCS:.f=.o)

#LIBDIR = \
#        -L/usr/pgi/linux86-64/7.0-5/lib

#LIBS = \
#        $(LIBDIR)
LIBS =

#F95 = pgf95
#F95 = ifort 
F95 = mpif90
 
#F95FLAGS = -O2 -r8 -tp amd64e
#F95FLAGS = -O2 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fopenmp -static-libgfortran  
F95FLAGS = -O3 -r8 
LDFLAGS =

all: $(PROG)

.PHONY: gnu ifx

gnu:
	mkdir -p $(BUILD_DIR)
	mpif90 -O2 -fdefault-real-8 -fdefault-double-8 \
		-ffixed-line-length-none -fallow-argument-mismatch \
		-finit-local-zero -mcmodel=medium -o $(BUILD_DIR)/chease-gnu.x $(SRCS)

# CHEASE depends on zero-initialized legacy local storage during startup.
# Without -init=zero the current ifx build segfaults before reading EQDATA.
ifx:
	mkdir -p $(BUILD_DIR)
	mpiifx -O3 -real-size 64 -xHost -init=zero -mcmodel=medium \
		-extend-source 132 -o $(BUILD_DIR)/chease-ifx.x $(SRCS)

$(PROG): $(OBJS)
	$(F95) $(F95FLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.d work.*

test:
	echo $(OBJS) 

tar:
	gtar zcvf ../$(PROG).tar.gz $(SRCS) $(INCS) $(MODS) Makefile_Portland README

.SUFFIXES: $(SUFFIXES) .f

.f.o:
	$(F95) $(F95FLAGS) -c $<
