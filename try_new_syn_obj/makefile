EXEC_DIR = ..
EXEC = $(EXEC_DIR)/simulator_new
SRC_DIR = cpp_sources
SRC_FILES = $(SRC_DIR)/*.cpp
#OBJ = main.o Neurons.o NeuronNetwork.o ElectricalSynapses.o ChemicalSynapses.o SimulatorInterface.o
OBJ = main.o Neurons.o NeuronNetwork.o ChemicalSynapses.o SimulatorInterface.o

###########################################################################
CXX = g++ # must use version that MATLAB supp

CXXFLAGS = -std=c++11

CXXDEBUGFLAGS = -Wall -g #-pg #-pg for gprof

CXXOPTIMFLAGS = #-O1 #higher level: -O2 or -O3

CXXLIBS  = 

COMPILE_THIS_ONE = $(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXOPTIMFLAGS) -c $<
###########################################################################
all: $(EXEC)
	

# The main idea is that compile each .o separately and then link them
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXOPTIMFLAGS) -O $(OBJ) -o $@ $(CXXLIBS)
	@echo "EXEC compiled"

main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/SimulatorInterface.h
	$(COMPILE_THIS_ONE)
	@echo "main.o updated"

SimulatorInterface.o: $(SRC_DIR)/SimulatorInterface.cpp $(SRC_DIR)/SimulatorInterface.h $(SRC_DIR)/NeuronNetwork.h
	$(COMPILE_THIS_ONE)
	@echo "SimulatorInterface.o updated"

NeuronNetwork.o: $(SRC_DIR)/NeuronNetwork.cpp $(SRC_DIR)/NeuronNetwork.h $(SRC_DIR)/Neurons.h $(SRC_DIR)/ChemicalSynapses.h $(SRC_DIR)/ElectricalSynapses.h
	$(COMPILE_THIS_ONE)
	@echo "NeuronNetwork.o updated"

Neurons.o: $(SRC_DIR)/Neurons.cpp $(SRC_DIR)/Neurons.h
	$(COMPILE_THIS_ONE)
	@echo "Neurons.o updated"

ChemicalSynapses.o: $(SRC_DIR)/ChemicalSynapses.cpp $(SRC_DIR)/ChemicalSynapses.h $(SRC_DIR)/Neurons.h
	$(COMPILE_THIS_ONE)
	@echo "ChemicalSynapses.o updated"

#ElectricalSynapses.o: $(SRC_DIR)/ElectricalSynapses.cpp $(SRC_DIR)/ElectricalSynapses.h $(SRC_DIR)/Neurons.h
#	$(COMPILE_THIS_ONE)
#	@echo "ElectricalSynapses.o updated"



#$(OBJ): $(SRC_FILES)
#	$(CXX) $(CXXFLAGS) $(arguments) $(CXXOPTIMFLAGS) -c $?
#	@echo "OBJ compiled"

clean: 
	rm *.o
	@echo "clean done"

clean_all:
	rm *.o $(EXEC)
	@echo "clean_all done"

	
