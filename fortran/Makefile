#
# Copyright (c) 2022 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

EXTERNAL_OBJ=cbessel.o
OBJS_MID=constants.o invlap.o utility.o 
MAIN=powerlaw-flow.o

TARG=powerlaw

OBJS=$(EXTERNAL_OBJ) $(OBJS_MID) $(MAIN)

SRCS=$(OBJS:.o=.f90)

#### gfortran #################
F90=gfortran
DFLAGS=-Og -Wall -Wextra  -g -fcheck=all 
PFLAGS=-O2 -march=native

# switch between performance or debugging flags
# based on whethere the environment variable PERF is defined
ifdef PERF
  FLAGS=$(PFLAGS)
else
  FLAGS=$(DFLAGS)
endif

# the link step
main: $(OBJS)
	$(F90) $(FLAGS) -o $(TARG) $(OBJS) 

%.o: %.f90
	$(F90) -c $(FLAGS) -o $@ $<

cbessel.o complex_bessel.mod: cbessel.f90
constants.o constants.mod: constants.f90
utility.o utility.mod: utility.f90 constants.mod 
invlap.o invlap.mod: invlap.f90 constants.mod
powerlaw-flow.o: powerlaw-flow.f90 constants.mod \
  utility.mod invlap.mod complex_bessel.mod

clean:
	rm -f *.o *.mod $(TARG)



