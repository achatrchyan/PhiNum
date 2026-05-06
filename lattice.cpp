#ifndef _LATTICECPP_
#define _LATTICECPP_

#include <stdlib.h>

#include "parameters.cpp"

namespace LocalLattice 
{
//Time
  pseudo_double tau;
  pseudo_double dtau;

#if (expansion==1)
  pseudo_double tau_conformal;
  pseudo_double scale_factor;
#endif
//Space  
  int N_t_loc;
  int N_eta_loc;
  int vol;
	
//Classical Fields
  pseudo_double *Phi;
  pseudo_double *Pi;
  
#if (withEnergyPS==1)
  pseudo_double *Delta;		//Density constant
#endif

  pseudo_double *Phi_backup;
  pseudo_double *Pi_backup;
  
  
//Averaged over local lattice quantities
  pseudo_double *phi, *pi; 
  pseudo_double pi2, phi2;
  pseudo_double phi4, phi6;
  pseudo_double phipi;
  
  pseudo_double Energy, Interaction, Charge, Charge2;
  
//Averaged over nodes
  pseudo_double *avgphi, *avgpi; 
  pseudo_double avgpi2, avgphi2;
  pseudo_double avgphi4, avgphi6;
  pseudo_double avgphipi;
  
  pseudo_double avgEnergy, avgInteraction, avgCharge, avgCharge2;
  
//Extremal values of the classical Field - over local lattice
  pseudo_double maxPhi,minPhi,maxAbsPhi;
  int maxpos,minpos,maxAbspos;

//Extremal values of the classical Field - over all nodes
  pseudo_double GlmaxPhi,GlminPhi,GlmaxAbsPhi;
  int GlmaxID,GlminID,GlmaxAbsID;
  int Glmaxpos, Glminpos, GlmaxAbspos;
 
//Counter terms
  pseudo_double counterMsqr_long;
  pseudo_double counterMsqr_trans;
  pseudo_double BareM_sqr_long;
  pseudo_double BareM_sqr_trans;
#if (theory==3)
  pseudo_double BareM_mono_sqr_long;
  pseudo_double BareM_mono_sqr_trans;
#endif
	
  void init()
  {
    N_t_loc=N_t;
    N_eta_loc=N_eta/nodes;
    vol=N_t_loc*N_t_loc*N_eta_loc;
		
//Phi and Pi will contain 2 additional border slices in eta direction
    Phi=new pseudo_double[Nc*N_t*N_t*(N_eta_loc+2)];
    Pi=new pseudo_double[Nc*N_t*N_t*(N_eta_loc+2)];
#if (withEnergyPS==1)
    Delta=new pseudo_double[Nc*N_t*N_t*(N_eta_loc+2)];
#endif
    for (int k = 0; k < Nc*N_t*N_t*(N_eta_loc+2); k++)
    {
      Phi[k]=0;
      Pi[k]=0;
#if (withEnergyPS==1)
      Delta[k]=0.;
#endif
    }
  
    phi=new pseudo_double[Nc];
    pi=new pseudo_double[Nc];
    
    avgphi=new pseudo_double[Nc];
    avgpi=new pseudo_double[Nc];
    
    
    counterMsqr_long = 0., BareM_sqr_long = M_sqr;
    counterMsqr_trans = 0., BareM_sqr_trans = M_sqr;
#if (relicpockets==1)
    BareM_sqr_trans = 0;
#endif
#if (theory==3)
    BareM_mono_sqr_long=M_mono_sqr;
    BareM_mono_sqr_trans=M_mono_sqr;
#endif
  }
  
	
  void fin()
  {
    delete[] Phi;
    delete[] Pi;
#if (withEnergyPS==1)
    delete[] Delta;
#endif
    
    delete[] phi;
    delete[] pi;
  }
}

using namespace LocalLattice;

void BackupConfig()
{
    Pi_backup = new pseudo_double[Nc*N_t*N_t*(N_eta_loc+2)];
    Phi_backup = new pseudo_double[Nc*N_t*N_t*(N_eta_loc+2)];
    
    for (int k = 0; k < Nc*N_t*N_t*(N_eta_loc+2); k++)
    {
        Phi_backup[k] = Phi[k];
        Pi_backup[k] = Pi[k];
    }
}


void RestoreConfig()
{
    for (int k = 0; k < Nc*N_t*N_t*(N_eta_loc+2); k++)
    {
        Phi[k] = Phi_backup[k];
        Pi[k] = Pi_backup[k];
    }
    delete[] Phi_backup;
    delete[] Pi_backup;
}
#endif
