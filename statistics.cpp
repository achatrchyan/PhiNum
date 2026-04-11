#ifndef _STATCPP_
#define _STATCPP_

#include <iostream>
#include <fstream> 
#include <iomanip>
#include <mpi.h>

#include "lattice.cpp"
#include "islands.cpp"

namespace Statistic
{
    pseudo_double T_decay;
    pseudo_double *T_decays;
    
    const int nnnnn = 200;
    
    long long islands_phi0_r[nnnnn][nnnnn];
    long long islands_phi0_r0[nnnnn][nnnnn];
 
    long long totislands_phi0_r[nnnnn][nnnnn];
    long long totislands_phi0_r0[nnnnn][nnnnn];

    void init()
    {
        T_decay=0;
        T_decays=new pseudo_double[runs];
		
        for (int j=0;j<nnnnn;j++)
            for (int k=0;k<nnnnn;k++)
            {
                islands_phi0_r[j][k]=0;
                islands_phi0_r0[j][k]=0;
                totislands_phi0_r[j][k]=0;
                totislands_phi0_r0[j][k]=0;                
            }
    }
}

using namespace Statistic;

void add_stat_new()
{
    using namespace LocalLattice;
    using namespace island_analysis;

    for (int i=0;i<found_islands.size();i++)
    {
        if (found_islands[i].disp)
        {
            pseudo_double fffff = fabs(found_islands[i].Peak_value.back());
            pseudo_double rrrrr = pow(found_islands[i].V.back()*3./4./acos(-1.), 1./3.);
            pseudo_double rrrrr0 = pow(found_islands[i].V.back()*3./4./acos(-1.), 1./3.)/sqrt( log( fffff* fffff / thres / thres) );

            if ((fffff<5*thres) && (rrrrr<5./QS))
                islands_phi0_r[int(fffff/5/thres*nnnnn)][int(rrrrr/5.*QS*nnnnn)]++;
            if ((fffff<5*thres) && (rrrrr0<5./QS))
                islands_phi0_r0[int(fffff/5/thres*nnnnn)][int(rrrrr0/5.*QS*nnnnn)]++;
        }
	} 

    MPI_Barrier(MPI_COMM_WORLD);
}

void output_stat()
{
    using namespace LocalLattice;
    using namespace island_analysis;

	MPI_Reduce(islands_phi0_r, totislands_phi0_r, nnnnn*nnnnn, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(islands_phi0_r0, totislands_phi0_r0, nnnnn*nnnnn, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  	MPI_Barrier(MPI_COMM_WORLD);
  
  	if (ID==0)
	{
    	std::ofstream fout;
    	char *fname=new char[256];
    	snprintf(fname,256,"./output/%i/islands_stat_2d.txt",LASTRUN);
    	fout.open(fname);
    	for (int i=0;i<nnnnn;i++)
        	for (int j=0;j<nnnnn;j++)
            	fout<<i*(5.*thres/nnnnn)<<' '<<j*(5./QS/nnnnn)<<' '<<totislands_phi0_r[i][j]/( 5.*thres/(pseudo_double)(nnnnn) )/ (5./QS/nnnnn) /(long double)(runs)/(long double)(N_t)/(long double)(N_t)/(long double)(N_t)/a_t/a_t/a_t<<' '<<totislands_phi0_r0[i][j]/(5.*thres/(pseudo_double)(nnnnn) )/(5./QS/nnnnn)/(long double)(runs)/(long double)(N_t)/(long double)(N_t)/(long double)(N_t)/a_t/a_t/a_t<<endl;
    	fout.close();

    	snprintf(fname,256,"./output/%i/islands_stat_1d.txt",LASTRUN);
    	fout.open(fname);
    	for (int i=0;i<nnnnn;i++)
    	{
       		long long s=0;
        	for (int j=0;j<nnnnn;j++)
            	s+=totislands_phi0_r[i][j];
        	fout<<i*(5*thres/nnnnn)<<' '<<(long double)(s)/(5*thres/nnnnn)/(long double)(runs)/(long double)(N_t)/(long double)(N_t)/(long double)(N_t)/a_t/a_t/a_t<<endl;
    	}
    	fout.close();

	    long long s =0;
   	 	for (int i=0;i<nnnnn;i++)
        	for (int j=0;j<nnnnn;j++)
           	 	s+=totislands_phi0_r[i][j];
    	cout<<s<<' '<<(long double)(s)/(long double)(runs)/(long double)(N_t)/(long double)(N_t)/(long double)(N_t)/a_t/a_t/a_t<<endl;
    
    	delete[] fname;
	}
}
 
void add_stat_time()
{
    T_decay+=LocalLattice::tau;
    T_decays[runID-1] = LocalLattice::tau;

    if (ID==0)
	cout<<LocalLattice::tau<<' '<<T_decay/runID<<endl;
}

void output_stat_time()
{
    sort(T_decays, T_decays+runs);
    if (ID==0)
    {
        std::ofstream fout;
    	char *fname=new char[256];
    	snprintf(fname,256,"./output/%i/decay_times.txt",LASTRUN);
    	fout.open(fname);
    	for (int i=0;i<runs;i++)
        	fout<<T_decays[i]<<' '<<runs - i<<endl;
    	fout.close();
        delete [] fname;

    	fname=new char[256];
    	snprintf(fname,256,"./output/%i/mean_decay_time.txt",LASTRUN);
    	fout.open(fname);
    	fout<<T_decay/runs<<endl;
    	fout.close();
        delete [] fname;
    }
}
 
#endif
