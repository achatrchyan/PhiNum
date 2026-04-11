#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include <mpi.h>

#include "comm.cpp"
#include "parameters.cpp"
#include "fieldanalysis.cpp"
#include "clusterfft.cpp"
#include "lattice.cpp"
#include "initialconditions.cpp"
#include "renormalization.cpp"
#include "dynamics.cpp"
#include "momentumspace.cpp"
#include "islands.cpp"
#include "filter.cpp"
#include "visualize.cpp"
#include "statistics.cpp"
int main(int argc, char **argv)
{
// start MPI and put out the number of processes involved
  MPI_Init(&argc,&argv);
  
  MPI::init();

  if (ID==0) 
    std::cout << "\nHi, Alex. Here are " << nodes << " process(es).\n";

//Input & Preparation
  GetParameters(argc, argv);
  CompleteParameters();

  Prepare_Output();
  Print_Info_Console();
  Print_Info_File();
  MPI_Barrier(MPI_COMM_WORLD);
  
  LocalLattice::init();
  MPI_Barrier(MPI_COMM_WORLD);
  
  FFT::init();
  MPI_Barrier(MPI_COMM_WORLD);

  MomentumSpace::init();
  MPI_Barrier(MPI_COMM_WORLD);  

  Statistic::init();
  MPI_Barrier(MPI_COMM_WORLD);  
  
  for (runID=1; runID<=runs; runID++)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (ID==0 && runs!=1)
      std::cout << "Run: " << runID << std::endl;
    
    fix_tadpole();
#if (islanddetect == 1)  
  island_analysis::init();
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    SetInitialConditions();
    MPI_Barrier(MPI_COMM_WORLD);
	
#if (islanddetect == 1)      
    CalculateVolumeAverages();
    island_analysis::thres = 2.57583 * sqrt(LocalLattice::avgphi2);
    Full_Analysis();
    add_stat_new();
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  //Dynamics  
    for (reset(); LocalLattice::tau<TMax; evolve() )
    {
      Print_AveragesExtrema();
       
#if (islanddetect == 1)
      if (numeric::total_stepcount%5==0)
      	Full_Analysis();
#endif
      if (false_vacuum_decayed())
      {
          add_stat_time();
          break;
      }
      Print_fp();
    }   
#if (islanddetect == 1)  
  island_analysis::fin();
  MPI_Barrier(MPI_COMM_WORLD);
#endif 
  }
  
#if (islanddetect == 1)  
  output_stat();
#endif
  output_stat_time();


  MomentumSpace::fin();
  MPI_Barrier(MPI_COMM_WORLD);
  
  FFT::fin();
  MPI_Barrier(MPI_COMM_WORLD); 

  LocalLattice::fin();
  MPI_Barrier(MPI_COMM_WORLD);

  if (ID==0) 
    std::cout << "Done!\n";
  
  MPI_Finalize();
}
