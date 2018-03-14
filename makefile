SHELL=/bin/sh

SRCS= \
CausalChargeDiffusion.cpp \
gauss/gauss_quadrature.cpp \
lib.cpp \
resonance/lib.cpp \
resonance/decay/readindata.cpp \
resonance/decay/decay.cpp \
resonance/sfn.cpp

HDRS= \
gauss/gauss_quadrature.h \
CCDparams.h \
defs.h \
lib.h \
asymptotics.h \
resonance/lib.h \
resonance/indexers.h \
resonance/decay/decay.h \
resonance/decay/readindata.h \
resonance/decay/parameters.h \
resonance/thermal.h \
resonance/chebyshev.h \
resonance/sfn.h

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

CausalChargeDiffusion.o : CausalChargeDiffusion.cpp defs.h resonance/sfn.h CCDparams.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c CausalChargeDiffusion.cpp -o CausalChargeDiffusion.o

resonance/sfn.o : resonance/sfn.cpp resonance/sfn.h resonance/decay/parameters.h
	$(CC) $(CFLAGS) $(WARNFLAGS) $(LIBS) -c resonance/sfn.cpp -o resonance/sfn.o

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
