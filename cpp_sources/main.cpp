#include "SimulatorInterface.h"
using namespace std;

int main(int argc, char* argv[]){// arguments should be input file path
	cout << "Number of input files: " << argc-1 << endl;
	for (int i = 1; i < argc; ++i){
		cout << "Processing input file No." << i << " out of " << argc-1 << "..." << endl;
		SimulatorInterface simulator;
		if (simulator.import(argv[i])){ // return true if import is successful
			simulator.simulate();
			simulator.output_results();
			cout << "Input file No." << i << " out of " << argc-1 << " processed." << endl;
		}
	}
	cout << "The planet earth is blue and there's nothing I can do." << endl;
	return 0; 
};
