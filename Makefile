SHELL	      = /bin/tcsh

FC	      = gfortran
CFLAGS	      = -O2 
LD	      = $(FC)

LDFLAGS	      = 

SRCS	      = Main.c \

SRCFOR	      =	corona.f \
		correl.f \
		intcor.f \
		interp.f \
		splinep.f \
		splinet.f \
		splines.f \
		collop.f  \
		collot.f \
		collos.f \
 		collop_sar.f \
 		collop_sar_weight.f \
		differ.f \
		intcman.f \
		intcman2.f \
		corgrid.f \
		preditp.f \
		preditt.f \
		predits.f \
		preditp_sar.f \
		preditp_grid.f

UTIL	       = model.f \
		algemod.f \
		ordimod.f \
		utilmod.f 


EXECS	      =	corona \
		correl \
		intcor \
		interp \
		splinep \
		splinet \
		splines \
		collop  \
		collot \
		collos \
		corgrid \
		intcman \
		intcman2 \
		differ \
		preditp \
		preditt \
		predits 


LIBS	  = model.o algemod.o ordimod.o utilmod.o

default:
		@echo "'make all' will create all executables"
		@echo "'make clean' will delete all execpt executables"
		@echo "'make clear' will delete all"
		@echo "'make <executable>' will crate only that executable"

all:		$(EXECS) 
		@echo "Done, you may run 'make clean'"
clean:
		rm -f *.o *.O core
clear:
		rm -f $(EXECS)
		rm -f *.o *.O core 
util:
		$(FC) $(FFLAGS) -c $(UTIL) 

corona:	 	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

correl:		$(LIBS) 
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

intcor:		$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

interp:		$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

splinep:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

splinet:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

splines:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

collop:		$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

collot:		$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

collos:		$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

preditp:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

preditt:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

predits: 	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

intcman:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f 

intcman2:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f 

differ:		$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)


corgrid:
		$(FC) $(LDFLAGS) -o $@ $@.f

collop_sar:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

preditp_sar:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)

collop_sar_weight:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)
corgr_m:	
		$(FC) $(LDFLAGS) -o $@ $@.f
crosgrid:	
		$(FC) $(LDFLAGS) -o $@ $@.f

preditp_grid:	$(LIBS)
		$(FC) $(LDFLAGS) -o $@ $@.f $(LIBS)
