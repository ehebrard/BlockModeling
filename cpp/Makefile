
MAINDIR = .

COPTIMIZE = -O3 

CCC = g++

BIN=$(MAINDIR)/bin
SRC=$(MAINDIR)/src/lib
MOD=$(MAINDIR)/main
OBJ=$(MAINDIR)/src/obj
INC=$(MAINDIR)/src/include


MODELS = $(wildcard $(MOD)/src/*.cpp)
BINS = $(patsubst $(MOD)/src/%, $(BIN)/%, $(MODELS:.cpp=))


PINCSRC = $(wildcard $(INC)/*.hpp)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))
XLIBSRC = $(wildcard $(XSRC)/*.cc)
XLIBAUX = $(XLIBSRC:.cc=.o)
XLIBOBJ = $(patsubst $(XSRC)/%, $(XOBJ)/%, $(XLIBAUX))


CFLAGS = -std=c++17 -D_LIBCPP_ENABLE_CXX17_REMOVED_AUTO_PTR -I$(INC)
LFLAGS = -L$(OBJ)


## Compile options
%.o:			CFLAGS +=$(COPTIMIZE) 
%.op:			CFLAGS +=$(COPTIMIZE) -pg -ggdb -D NDEBUG
%.od:			CFLAGS +=-O0 -ggdb -D DEBUG -D INVARIANTS #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
%.or:			CFLAGS +=$(COPTIMIZE) -D NDEBUG
%.oc:     CFLAGS +=-O0 -fprofile-arcs -ftest-coverage -ggdb -D DEBUG


#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------

all: lib $(BINS) 

# The library
lib: $(PLIBOBJ) 
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(PLIBOBJ) $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<	
	$(CCC) $(CFLAGS) $(PLIBOBJ) $< -lm -o $(BIN)/$@ 

DATE := $(shell date '+%y-%m-%d')

clean : 
	rm -rf $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~  *~ $(BIN)/* 

