SHELL=/bin/sh

SRCS= \
main.cpp \
gauss/gauss_quadrature.cpp \
lib.cpp \
decay/readindata.cpp \
decay/decay.cpp

HDRS= \
gauss/gauss_quadrature.h \
lib.h \
defs.h \
decay/decay.h \
decay/readindata.h \
decay/parameters.h \
asymptotics.h \
thermal.h \
chebyshev.h

MAKEFILE=makefile

COMMAND=run_test

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -pg -g -O3 -fopenmp

WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas
LIBS= -L/sw/lib -I/sw/include

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(COMMAND)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
