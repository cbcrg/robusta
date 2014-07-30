



CC = gcc

EXECUTABLE = robusta

WFLAGS = -Wall -Wshadow -D_GNU_SOURCE -Wextra

LIB_FLAGS = -lm

OBJ_DIR = obj
SOURCE_DIR = src


_OBJ = main.o parser.o ms_aligners.o pw_aligners.o system_caller.o util.o low_coverage.o classes.o low_coverage_concatenate.o
OBJ = $(patsubst %,$(OBJ_DIR)/%,$(_OBJ))



# Defines additional flags according to user

# uses pedantic mode for compilation
ifeq ($(pedantic),y)
	PEDANTIC_FLAGS := -pedantic -std=c99 -D_XOPEN_SOURCE
else
	PEDANTIC_FLAGS :=
endif


# Turns on release configuration
ifeq ($(mode),release)
	ADD_FLAGS := -O3
else ifeq ($(mode),profile)
	ADD_FLAGS := -O3 -pg
else
	ADD_FLAGS := -g
endif


# Enables parallel compilation (omg.h needed) only applyes for multiple aligners pairwise aligners
ifeq ($(parallel),y)
	USE_PARALLEL := -DPARALLEL -fopenmp
else
	USE_PARALLEL :=
endif


$(EXECUTABLE): $(OBJ)
	$(CC) $(WFLAGS) $(ADD_FLAGS) $(PEDANTIC_FLAGS) $(OBJ)  -o $(EXECUTABLE) $(LIB_FLAGS) $(USE_PARALLEL)
	cp $(EXECUTABLE) $(USER_BIN)/$(EXECUTABLE)


$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.c
	$(CC) $(WFLAGS) $(ADD_FLAGS) $(PEDANTIC_FLAGS) $(USE_PARALLEL) -c -o $@ $<





# Delete object files, documentation files and linux generated backupfiles
.PHONY: clean_all
clean_all: clean clean_doc
	rm -rf *~


# Deletes object files and executable
.PHONY: clean
clean:
	rm -f obj/*.o
	rm -f $(EXECUTABLE)


# Deletes the documentation
.PHONY: clean_doc
clean_doc:
	rm -rf doc/html doc/latex


# Generates a doxygen documentation of the source code
.PHONY: doc
doc:
	doxygen doc/doxygen_config