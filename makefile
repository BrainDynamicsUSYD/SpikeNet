EXEC_DIR = ..
EXEC = $(EXEC_DIR)/simulator
SRC_DIR = cpp_sources
SRC_FILES = $(SRC_DIR)/*.cpp
HEADER_FILES = $(SRC_DIR)/*.h
#OBJ = main.o NeuroPop.o NeuroNet.o ElectricalSynapses.o ChemSyn.o SimulatorInterface.o
OBJ = main.o NeuroPop.o NeuroNet.o ChemSyn.o SimuInterface.o

###########################################################################
CXX = g++ # must use version that MATLAB supp

CXXFLAGS = -std=c++11 -Wall -Wextra

CXXDEBUGFLAGS = -Wall -g #-pg #-pg for gprof

CXXOPTIMFLAGS = #-O1 #higher level: -O2 or -O3

CXXINCLUDE = -I/usr/local/include 

CXXLIBS =  -L/usr/local/lib/ -lhdf5_cpp -lhdf5

COMPILE_THIS_ONE = $(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXOPTIMFLAGS) $(CXXINCLUDE) -c $<
###########################################################################
#include includes.mk
all: $(EXEC)
	

# The main idea is that compile each .o separately and then link them
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXOPTIMFLAGS) -O $(OBJ) -o $@  $(CXXLIBS)
	@echo "EXEC compiled"

main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/SimuInterface.h
	$(COMPILE_THIS_ONE)
	@echo "main.o updated"

SimuInterface.o: $(SRC_DIR)/SimuInterface.cpp $(SRC_DIR)/SimuInterface.h $(SRC_DIR)/NeuroNet.h
	$(COMPILE_THIS_ONE)
	@echo "SimuInterface.o updated"

NeuroNet.o: $(SRC_DIR)/NeuroNet.cpp $(SRC_DIR)/NeuroNet.h $(SRC_DIR)/ChemSyn.h
	$(COMPILE_THIS_ONE)
	@echo "NeuroNet.o updated"

NeuroPop.o: $(SRC_DIR)/NeuroPop.cpp $(SRC_DIR)/NeuroPop.h
	$(COMPILE_THIS_ONE)
	@echo "NeuroPop.o updated"

ChemSyn.o: $(SRC_DIR)/ChemSyn.cpp $(SRC_DIR)/ChemSyn.h $(SRC_DIR)/NeuroPop.h
	$(COMPILE_THIS_ONE)
	@echo "ChemSyn.o updated"

#ElectricalSynapses.o: $(SRC_DIR)/ElectricalSynapses.cpp $(SRC_DIR)/ElectricalSynapses.h $(SRC_DIR)/NeuroPop.h
#	$(COMPILE_THIS_ONE)
#	@echo "ElectricalSynapses.o updated"

#documentation
docs: html/index.html
html/index.html: ${SRC_FILES} ${HEADER_FILES} Doxyfile
	doxygen Doxyfile

#$(OBJ): $(SRC_FILES)
#	$(CXX) $(CXXFLAGS) $(arguments) $(CXXOPTIMFLAGS) -c $?
#	@echo "OBJ compiled"

clean: 
	rm *.o
	@echo "clean done"

clean_all:
	rm *.o $(EXEC)
	@echo "clean_all done"

	
