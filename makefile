SHELL=/bin/sh

SRCS= \
CausalChargeDiffusion.cpp \
gauss/gauss_quadrature.cpp \
lib.cpp

HDRS= \
gauss/gauss_quadrature.h \
defs.h \
lib.h \
resonance/sfn.h \
asymptotics.h

MAKEFILE=makefile

COMMAND=run_CCD

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS=  -pg -g -O3
#WARNFLAGS= -ansi -pedantic -Wall -W \
#   -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
#   -Wwrite-strings -fshort-enums -fno-common 
WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include
 
$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

CausalChargeDiffusion.o : CausalChargeDiffusion.cpp defs.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c CausalChargeDiffusion.cpp -o CausalChargeDiffusion.o

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
