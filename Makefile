# compiler
FC = gfortran
CC = gcc

# compile flags
FCFLAGS = -g -Ofast -ffixed-line-length-0 -std=legacy
CFLAGS = -g -Ofast -Wno-pointer-to-int-cast

SRC_DIR = src
OBJ_DIR = lib
BIN_DIR = bin

# source files and objects
FSRCS = $(wildcard $(SRC_DIR)/*.f90)
CSRCS = $(wildcard $(SRC_DIR)/*.c)

FOBJ = $(FSRCS:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
COBJ = $(CSRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# program name
PROGRAM = aao_rad

MKDIR_P = mkdir -p
.PHONY: directories


all: directories $(COBJ) $(FOBJ) $(PROGRAM)

directories: ${OBJ_DIR} ${BIN_DIR}

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

$(PROGRAM): $(FOBJ) $(COBJ)
	$(FC) $(FLFLAGS) -o bin/$@_lund $(COBJ) $(FOBJ)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -f lib/*.o bin/aao_rad
