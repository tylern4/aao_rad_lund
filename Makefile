# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -Ofast -ffixed-line-length-0 -std=legacy
# link flags
FLFLAGS =

# source files and objects
SRCS = $(patsubst %.F, %.o, $(wildcard *.F))

# program name
PROGRAM = aao_rad

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	gcc -g -Ofast -c unixtime.c
	gcc -g -Ofast -c lenocc.c
	$(FC) $(FLFLAGS) -o $(PROGRAM) $(SRCS) unixtime.o lenocc.o

%.o: %.F
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod
