#ifndef _COMMCPP_
#define _COMMCPP_

#include <mpi.h>
#include <cstring>

#include "lattice.cpp"
#include "util.cpp"


/*=======CONTAINS MESSAGE PASSING INFORMATION========*/
namespace MPI 
{
  int left;
  int right;
	
  void init()
  {
  //Get MPI Information
    MPI_Comm_rank(MPI_COMM_WORLD,&ID);
    MPI_Comm_size(MPI_COMM_WORLD,&nodes);
		
  //Set Neighbours
    left=ID-1;
    right=(ID+1) % nodes;
    if (left==-1)
      left=nodes-1;
  }
}

/*=======CONTAINS SYNCHRONIZATION OF BOUNDARY ROUTINES======*/
namespace Communication 
{
	/*Local MPI status*/
	MPI_Status status;
	
	/*MPI Transfer requests*/
	MPI_Request RecvRL;
	MPI_Request SendRL;
	MPI_Request RecvLR;
	MPI_Request SendLR;
		
	void requestRightLeftReceive() //z=N+1
	{
	  MPI_Irecv(&LocalLattice::Phi[pos(0,0,LocalLattice::N_eta_loc+1)], N_t*N_t*Nc, MPI_PSEUDO_DOUBLE, MPI::right, 21, MPI_COMM_WORLD, &RecvRL);		
	}
	
	void requestRightLeftSend() //z=1
	{
	  MPI_Isend(&LocalLattice::Phi[pos(0,0,1)], N_t*N_t*Nc, MPI_PSEUDO_DOUBLE, MPI::left, 21, MPI_COMM_WORLD, &SendRL);		
	}
	
	void requestLeftRightReceive() // z=0
	{
	  MPI_Irecv(&LocalLattice::Phi[pos(0,0,0)], N_t*N_t*Nc, MPI_PSEUDO_DOUBLE, MPI::left, 12, MPI_COMM_WORLD, &RecvLR);
	}
	
	void requestLeftRightSend() // z=N
	{
	  MPI_Isend(&LocalLattice::Phi[pos(0,0,LocalLattice::N_eta_loc)], N_t*N_t*Nc, MPI_PSEUDO_DOUBLE, MPI::right, 12, MPI_COMM_WORLD, &SendLR);
	}
	
	void exchange()
	{
	  /* Exchange L<-R */
	  requestRightLeftReceive();
	  requestRightLeftSend();
	  MPI_Wait(&RecvRL,&status);
	  MPI_Wait(&SendRL,&status);
	  MPI_Barrier(MPI_COMM_WORLD);

	  /* Exchange L->R */
	  requestLeftRightReceive();
	  requestLeftRightSend();
	  MPI_Wait(&RecvLR,&status);
	  MPI_Wait(&SendLR,&status);	
	  MPI_Barrier(MPI_COMM_WORLD);
	}	
}

/*=======CONTAINS MASTER SLAVE TRANSFER ROUTINES=======*/
namespace MSCommunication 
{
	
	void sendResults(pseudo_double *localData,int localSize){
		
		//PREPARE TRANSFER
		MPI_Status MSstatus;
		MPI_Request MSSend;
		
		//SEND TO MASTER
		MPI_Isend(localData,localSize,MPI_PSEUDO_DOUBLE,0,99,MPI_COMM_WORLD,&MSSend);
		MPI_Wait(&MSSend,&MSstatus);
		
	}
	
	void receiveResults(pseudo_double *localData,int localSize,pseudo_double *globalData){
		
		//PREPARE TRANSFER
		MPI_Status MSstatus;
		MPI_Request *MSRecv=new MPI_Request[nodes];
		
		//DATA FROM LOCAL NODE
		std::memcpy(&globalData[0],localData,localSize*sizeof(pseudo_double));
		
		//DATA FROM OTHER NODES
		//RECEIVE FROM SLAVES
		for(int n=1;n<nodes;n++){
			MPI_Irecv(&globalData[n*localSize],localSize,MPI_PSEUDO_DOUBLE,n,99,MPI_COMM_WORLD,&MSRecv[n]);
		}
		
		for(int n=1;n<nodes;n++){
			MPI_Wait(&MSRecv[n],&MSstatus);
		}
		
		delete[] MSRecv;
		
	}
}

#endif
