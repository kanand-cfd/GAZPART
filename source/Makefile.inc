
FC	= mpiifort
LD	= $(FC)
F90	= $(FC)
SWP	= 
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3 -unroll -xAVX -fno-alias -no-prec-div -no-prec-sqrt -ip 
FOPTS	= -r8  # -g -lstdc++
FFLAGS	= $(FOPTS) $(OLEVEL)
FFTWDIR = $(HOME)/p3dfft-2.7.9_lib
LIBS    = $(FFTWDIR)/lib/libp3dfft.a -mkl
INCLUDE =-I$(FFTWDIR)/include
LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=
