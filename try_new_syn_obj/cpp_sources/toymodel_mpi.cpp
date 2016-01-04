#include "mpi.h"
#include <iostream>
#include <unistd.h> 
#include <stdio.h>
#include <vector>
#include <limits.h> // for INT_MAX, INT_MIN

using namespace std;

int main(int argc, char **argv){


	// Initiate MPI environment and do error check
	int ierr;
	ierr = MPI_Init( &argc, &argv );//MPI::Init();
	if (ierr != MPI_SUCCESS) {
		cout << "Error starting MPI program. Terminating." << endl;
		MPI_Abort( MPI_COMM_WORLD, ierr );
	}
	
	// Get communicator size (number of tasks) and current process rank (task id)
	int size, rank;
	MPI_Comm_size( MPI_COMM_WORLD, &size ); //size = MPI::COMM_WORLD.Get_size();
	MPI_Comm_rank( MPI_COMM_WORLD, &rank ); //rank = MPI::COMM_WORLD.Get_rank();
	
	// Extract the original group handle
	MPI_Group orig_group; // group handle
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group); 

	// Get hostname
	char hostname[256], ss[256]; 
	size_t len = 240;
	gethostname(hostname,len);

	// Cout current process rank
	printf("Hello world! I am %d of %d on %s\n", rank, size, hostname);

	/*******************************************************************************/
	/***********************Do some real communications here************************/
	/*******************************************************************************/
	

	// Master process reads input data and create Pop&Syn for each process
	int Num_pop = 2;

	int 
		/*******************************************************************************/
		// pop rank in MPI_COMM_WORLD must be smaller!!!!!!!!!!!
		/*******************************************************************************/
		IsPop[7] = {MPI_UNDEFINED, 1,1, MPI_UNDEFINED,MPI_UNDEFINED,MPI_UNDEFINED,MPI_UNDEFINED},
		IsSyn[7] = {MPI_UNDEFINED, MPI_UNDEFINED,MPI_UNDEFINED, 1, 1, 1, 1},

		PopInd[7]    = {MPI_UNDEFINED, 0,1, MPI_UNDEFINED,MPI_UNDEFINED,MPI_UNDEFINED,MPI_UNDEFINED},
		PrePop[7]    = {MPI_UNDEFINED, MPI_UNDEFINED,MPI_UNDEFINED,  0,1,0,1},  // pre-synaptic population for syn, pop index for pop
		PostPop[7]   = {MPI_UNDEFINED, MPI_UNDEFINED,MPI_UNDEFINED,  0,0,1,1},  // post-synaptic population for syn, pop index for pop

		color_popsyn[7] = {MPI_UNDEFINED,0,0,1,1,1,1},
		color_pre[7]    = {MPI_UNDEFINED,0,1,0,1,0,1},
		color_post[7]   = {MPI_UNDEFINED,0,1,0,0,1,1};

	double 
		self_data[7] = {MPI_UNDEFINED, 0.001, 0.02, 0.3, 4.0, 50.0, 600.0},
		recv_data[7];
		// the above are all local data, here they are global just for convenience in initialisation
		// 7 processors are needed
		


	
	// Group/Communicator management
	/*
	INPUT PARAMTERS
	local_comm
              - Local (intra)communicator
	local_leader
              - Rank in local_comm of leader (often 0)
	peer_comm/bridge_comm
              - a parent communicator directly contains all the local_comm
		Use dedicated bridge_comm to avoid interference
	remote_leader
              - Rank in peer_comm of remote leader, this is a knowledge must be known by all the local_comm
	tag    - Message tag to use in constructing intercommunicator
	*/
	int key = rank;
	int tag_shift = 100; // to distinct from previous tags
	// Pre_comm
	MPI_Comm PopSyn_pre_comm, Pre_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color_pre[rank], key, &PopSyn_pre_comm);
	if (IsPop[rank] == 1){
		MPI_Comm Pop_pre_comm;
		MPI_Comm_split(PopSyn_pre_comm, color_popsyn[rank], key, &Pop_pre_comm);
		int tag_pre = PopInd[rank] + tag_shift;
		MPI_Intercomm_create(Pop_pre_comm, 0, PopSyn_pre_comm, 1, tag_pre, &Pre_comm);
	}
	else if (IsSyn[rank] == 1){
		MPI_Comm Syn_pre_comm;
		MPI_Comm_split(PopSyn_pre_comm, color_popsyn[rank], key, &Syn_pre_comm);
		int tag_pre = PrePop[rank] + tag_shift;
		MPI_Intercomm_create(Syn_pre_comm, 0, PopSyn_pre_comm, 0, tag_pre, &Pre_comm);
	}
	else Pre_comm = MPI_COMM_NULL;	
	// Post_comm
	MPI_Comm PopSyn_post_comm, Post_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color_post[rank], key, &PopSyn_post_comm);
	if (IsPop[rank] == 1){
		MPI_Comm Pop_post_comm;
		MPI_Comm_split(PopSyn_post_comm, color_popsyn[rank], key, &Pop_post_comm);
		int tag_post = PopInd[rank] + tag_shift;
		MPI_Intercomm_create(Pop_post_comm, 0, PopSyn_post_comm, 1, tag_post, &Post_comm);
	}
	else if (IsSyn[rank] == 1){
		MPI_Comm Syn_post_comm;
		MPI_Comm_split(PopSyn_post_comm, color_popsyn[rank], key, &Syn_post_comm);
		int tag_post = PostPop[rank] + tag_shift;
		MPI_Intercomm_create(Syn_post_comm, 0, PopSyn_post_comm, 0, tag_post, &Post_comm);
	}
	else Post_comm = MPI_COMM_NULL;
	

	// Now test the intercommunicator
	if (Pre_comm != MPI_COMM_NULL){
		/*If comm is an intercommunicator, then the call involves all processes in the intercommunicator,
		but with one group (group A) defining the root process.
		All processes in the other group (group B) pass the rank of the root in group A. 
		The root passes the value MPI_ROOT in root.
		All other processes in group A pass the value MPI_PROC_NULL in root.
		Data is broadcast from the root to all processes in group B.
		The receive buffer arguments of the processes in group B must be consistent with the send buffer argument of the root.*/
		if (IsPop[rank] == 1){
			MPI_Bcast(&self_data[rank], 1, MPI_DOUBLE, MPI_ROOT, Pre_comm);
		}
		else if (IsSyn[rank] == 1){
			MPI_Bcast(&recv_data[rank], 1, MPI_DOUBLE, 0, Pre_comm);
			cout << "Syn at " << rank << " received data: " << recv_data[rank] << endl;
		}
	}
	if (Post_comm != MPI_COMM_NULL){
		if (IsPop[rank] == 1){
			MPI_Reduce(&recv_data[rank], &recv_data[rank], 1, MPI_DOUBLE, MPI_SUM, MPI_ROOT, Post_comm);
			cout << "Pop at " << rank << " received data: " << recv_data[rank] << endl;
		}
		else if (IsSyn[rank] == 1){
			MPI_Reduce(&self_data[rank], &self_data[rank], 1, MPI_DOUBLE, MPI_SUM, 0, Post_comm);
		}
	}


	/*******************************************************************************/
	/*******************************************************************************/

	// Close MPI environment
	// No more MPI routines should be called after this
	MPI_Finalize( ); //MPI::Finalize();

	// Master process says goodbye
	if (rank == 0){
		cout << "Master process sends its regards!" << endl;
	}
	return 0;
}


