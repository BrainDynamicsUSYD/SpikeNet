#include "SimuInterface.h"
using namespace std;

int main(int argc, char* argv[]){// arguments should be input file path
	cout << "Number of input files: " << argc-1 << endl;
	for (int i = 1; i < argc; ++i){
		cout << "Processing input file No." << i << " out of " << argc-1 << "..." << endl;
		SimuInterface simulator;
		
#ifdef HDF5
		bool success = simulator.import_HDF5(argv[i]);
#else
		bool success = simulator.import(argv[i]);
#endif
		
		if (success){ // return true if import is successful
			simulator.simulate();
			cout << "Input file No." << i << " out of " << argc-1 << " processed." << endl;
		}
	}
	cout << "The planet earth is blue and there's nothing I can do." << endl;
	return 0; 
};
