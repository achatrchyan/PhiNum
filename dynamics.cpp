#ifndef _DYNAMICSCPP_
#define _DYNAMICSCPP_

#include "lattice.cpp"
#include "comm.cpp"
#include "parameters.cpp"
#include "fieldanalysis.cpp"

namespace numeric
{
  int regime;		//0 - constant time step, 1-adaptive time step
  
  int total_stepcount;
  int normal_stepcount, adaptive_stepcount;
    
  pseudo_double Last_Field_update,Last_Momenta_update;
#if (expansion==1)
  pseudo_double Last_Field_update_conformal,Last_Momenta_update_conformal;
#endif
  
  pseudo_double LastMOoutput, LastSOoutput;
#if (expansion==1)
  pseudo_double LastMOoutput_conformal,LastSOoutput_conformal;
#endif
}

using namespace numeric;

void reset()
{
  regime=0;
  
  total_stepcount=0;
  normal_stepcount=0;
  adaptive_stepcount=0;
  
  LocalLattice::tau=T0;
#if (expansion==1)
  LocalLattice::tau=T0;
  LocalLattice::tau_conformal=( sqrt(1.+2.*H_0*T0) - 1. )/H_0;
  LocalLattice::scale_factor = sqrt(1.+2.*H_0*T0);
#endif
  LocalLattice::dtau= dtBaseTrans * a_t;
  Last_Field_update=LocalLattice::tau;
  Last_Momenta_update=LocalLattice::tau;  
  
#if (expansion==1)
  Last_Field_update_conformal=LocalLattice::tau_conformal;
  Last_Momenta_update_conformal=LocalLattice::tau_conformal;  
#endif
}

void updateMomenta()
{
  int Kooo=pos(0,0,1);
  int Kuoo,Kdoo,Kouo,Kodo,Koou,Kood;

  using namespace LocalLattice;
	
  for(int z=1;z<N_eta_loc+1;z++)
    for(int y=0;y<N_t_loc;y++)
      for(int x=0;x<N_t_loc;x++)
      {
//Calculate phi^2
      	pseudo_double locPhiSqr=0.0;
	      for (int a=0;a<Nc;a++)
	        locPhiSqr+=sqr( Phi[ Kooo+a ] );
	
      //Update neighbouring points for laplacian
	      Kuoo=pos(mod(x+1,N_t_loc), y, z);
	      Kouo=pos(x, mod(y+1,N_t_loc), z);
	      Koou=pos(x, y, z+1);

	      Kdoo=pos(mod(x-1,N_t_loc), y, z);
	      Kodo=pos(x, mod(y-1,N_t_loc), z);
	      Kood=pos(x, y, z-1);

      //Update Pi
	      for (int a=0;a<Nc;a++)
	      {
	        pseudo_double LapPhi = (Phi[Kuoo] + Phi[Kdoo] + Phi[Kouo] + Phi[Kodo] + Phi[Koou] + Phi[Kood] - 6.*Phi[Kooo])/sqr(a_t);
#if (theory==0)		  
#if (expansion==0)
	        if (a==0)
	          Pi[Kooo]+= ( LapPhi - ( Mass_Mode * BareM_sqr_long + Lambda_Mode * (locPhiSqr/(6.*Nc)) + (sqr(locPhiSqr/Kappa)/(120.*Nc*Nc)) ) * Phi[Kooo] )    *   (tau-Last_Momenta_update);
	        else
	          Pi[Kooo]+= ( LapPhi - ( Mass_Mode * BareM_sqr_trans + Lambda_Mode * (locPhiSqr/(6.*Nc)) + ( sqr( locPhiSqr/Kappa ) / (120.*Nc*Nc) ) ) * Phi[Kooo] )   *   (tau-Last_Momenta_update);
#endif
#if (expansion==1)
	        if (a == 0)
	          Pi[Kooo]+= ( LapPhi - ( Mass_Mode * BareM_sqr_long*scale_factor*scale_factor + Lambda_Mode * (locPhiSqr/(6.*Nc)) + (sqr(locPhiSqr/Kappa/scale_factor)/(120.*Nc*Nc)) ) * Phi[Kooo] )    *   (tau_conformal-Last_Momenta_update_conformal);
	        else
	          Pi[Kooo]+= ( LapPhi - ( Mass_Mode * BareM_sqr_trans*scale_factor*scale_factor + Lambda_Mode * (locPhiSqr/(6.*Nc)) + ( sqr( locPhiSqr/Kappa/scale_factor ) / (120.*Nc*Nc) ) ) * Phi[Kooo] )     *   (tau_conformal-Last_Momenta_update_conformal);
#endif
#endif
#if (theory==1)
#if (expansion==0)
	        if (a==0)
	          Pi[Kooo]+= ( LapPhi - h - ( Mass_Mode * BareM_sqr_long + Lambda_Mode * (locPhiSqr/(6.*Nc)) ) * Phi[Kooo] )    *   (tau-Last_Momenta_update);
	        else
	          Pi[Kooo]+= ( LapPhi - ( Mass_Mode * BareM_sqr_trans + Lambda_Mode * (locPhiSqr/(6.*Nc)) ) * Phi[Kooo] )   *   (tau-Last_Momenta_update);
#endif
#if (expansion==1)
	        if (a == 0)
	          Pi[Kooo] += (LapPhi - h*scale_factor*scale_factor*scale_factor - ( Mass_Mode * BareM_sqr_long * scale_factor*scale_factor + Lambda_Mode * (locPhiSqr/(6.*Nc)) ) * Phi[Kooo] )     *   (tau_conformal-Last_Momenta_update_conformal);
	        else
	          Pi[Kooo] += (LapPhi - ( Mass_Mode * BareM_sqr_trans* scale_factor*scale_factor  + Lambda_Mode * (locPhiSqr/(6.*Nc)) ) * Phi[Kooo])    *   (tau_conformal-Last_Momenta_update_conformal);
#endif
#endif
#if (theory==2)
#if (expansion==0)
	        if (a == 0)
	          Pi[Kooo] += (LapPhi - BareM_sqr_long*sqrt(BareM_sqr_long) * (Phi[Kooo]/sqrt(locPhiSqr))* sin(sqrt(locPhiSqr) / sqrt(BareM_sqr_long)))    *   (tau - Last_Momenta_update);
	        else
	          Pi[Kooo] += (LapPhi - BareM_sqr_trans*sqrt(BareM_sqr_trans)* (Phi[Kooo]/sqrt(locPhiSqr)) * sin(sqrt(locPhiSqr) / sqrt(BareM_sqr_trans)))   *   (tau - Last_Momenta_update);  
#endif
#if (expansion==1)
	        if (a == 0)
	          Pi[Kooo] += (LapPhi - BareM_sqr_long*sqrt(BareM_sqr_long) *scale_factor*scale_factor*scale_factor * (Phi[Kooo]/sqrt(locPhiSqr))* sin(sqrt(locPhiSqr)/scale_factor / sqrt(BareM_sqr_long)))     *   (tau_conformal-Last_Momenta_update_conformal);
	        else
	          Pi[Kooo] += (LapPhi - BareM_sqr_trans*sqrt(BareM_sqr_trans)*scale_factor*scale_factor*scale_factor * (Phi[Kooo]/sqrt(locPhiSqr)) * sin(sqrt(locPhiSqr)/scale_factor / sqrt(BareM_sqr_trans)))    *   (tau_conformal-Last_Momenta_update_conformal);
#endif
#endif
#if (theory==3)
#if (expansion==0)
	        if (a == 0)
	          Pi[Kooo] += (LapPhi - BareM_mono_sqr_long * Phi[Kooo] - BareM_sqr_long*sqrt(BareM_sqr_long) * (Phi[Kooo]/sqrt(locPhiSqr)) * sin(sqrt(locPhiSqr) / sqrt(BareM_sqr_long)))    *   (tau - Last_Momenta_update);
	        else
	          Pi[Kooo] += (LapPhi - BareM_mono_sqr_trans * Phi[Kooo] - BareM_sqr_trans*sqrt(BareM_sqr_trans) * (Phi[Kooo]/sqrt(locPhiSqr)) *sin(sqrt(locPhiSqr) / sqrt(BareM_sqr_trans)))   *   (tau - Last_Momenta_update);  
#endif
#if (expansion==1)
	        if (a == 0)
	          Pi[Kooo] += (LapPhi - BareM_mono_sqr_long * scale_factor*scale_factor * Phi[Kooo] - BareM_sqr_long*sqrt(BareM_sqr_long) * scale_factor*scale_factor*scale_factor * sin(Phi[Kooo] / (scale_factor*sqrt(BareM_sqr_long))))     *   (tau_conformal-Last_Momenta_update_conformal);
	        else
	          Pi[Kooo] += (LapPhi - BareM_mono_sqr_trans * scale_factor*scale_factor * Phi[Kooo] - BareM_sqr_trans*sqrt(BareM_sqr_trans) *scale_factor*scale_factor*scale_factor * sin(Phi[Kooo] /  (scale_factor*sqrt(BareM_sqr_long))))    *   (tau_conformal-Last_Momenta_update_conformal);
#endif
#endif
#if (theory==4)
#if (expansion==0)
	        if (a == 0)
		        Pi[Kooo] += (LapPhi - ( Mass_Mode * BareM_sqr_long - g* Phi[Kooo] / 2. + Lambda_Mode * (locPhiSqr / (6.*Nc))) * Phi[Kooo] )    *   (tau - Last_Momenta_update);
	        else
		        Pi[Kooo] += (LapPhi - ( Mass_Mode * BareM_sqr_trans - g* Phi[Kooo] / 2. + Lambda_Mode * (locPhiSqr / (6.*Nc))) * Phi[Kooo] )   *   (tau - Last_Momenta_update);
#if (relicpockets == 1)
          if (a == 0)
		        Pi[Kooo] += (LapPhi - ( BareM_sqr_long - g* Phi[Kooo] / 2. + (Phi[Kooo]*Phi[Kooo] / 6.) ) * Phi[Kooo] - 0.25*M_axion_sqr/Deltaphi/cosh((Phi[Kooo]-phimax)/Deltaphi)/cosh((Phi[Kooo]-phimax)/Deltaphi)*Phi[Kooo + 1]*Phi[Kooo+1] )    *   (tau - Last_Momenta_update);
//		        Pi[Kooo] += (LapPhi - ( BareM_sqr_long - g* Phi[Kooo] / 2. + (Phi[Kooo]*Phi[Kooo] / 6.) ) * Phi[Kooo] - exp(a_xi*(Phi[Kooo]-2.))*a_xi*Phi[Kooo + 1]*Phi[Kooo+1] )    *   (tau - Last_Momenta_update);
	        else
		        Pi[Kooo] += (LapPhi - M_axion_sqr * 0.5*(1+tanh((Phi[Kooo-1]-phimax)/Deltaphi)) * Phi[Kooo] )   *   (tau - Last_Momenta_update);
//		        Pi[Kooo] += (LapPhi - (2. * exp(a_xi*(Phi[Kooo-1]-2.))) * Phi[Kooo] )   *   (tau - Last_Momenta_update);
#endif 
#endif
#if (expansion==1)
	        if (a == 0)
		        Pi[Kooo] += (LapPhi - ( Mass_Mode * BareM_sqr_long *scale_factor*scale_factor - g*scale_factor* Phi[Kooo] / 2. + Lambda_Mode * (locPhiSqr / (6.*Nc))) * Phi[Kooo] )     *   (tau_conformal - Last_Momenta_update_conformal);
	        else
		        Pi[Kooo] += (LapPhi - ( Mass_Mode * BareM_sqr_long *scale_factor*scale_factor - g*scale_factor* Phi[Kooo] / 2. + Lambda_Mode * (locPhiSqr / (6.*Nc))) * Phi[Kooo] )     *   (tau_conformal - Last_Momenta_update_conformal);
#endif
#endif
#if (theory==5)
#if (expansion==0)
	        if (a == 0)
		        Pi[Kooo] += (LapPhi -  BareM_sqr_long * Phi[Kooo] * pow((1.+ locPhiSqr/ BareM_sqr_long), (p_mono - 1)) )    *   (tau - Last_Momenta_update);
	        else
		        Pi[Kooo] += (LapPhi -  BareM_sqr_trans * Phi[Kooo] * pow((1. + locPhiSqr / BareM_sqr_long), (p_mono - 1)) )   *   (tau - Last_Momenta_update);
#endif
#if (expansion==1)
	        if (a == 0)
		        Pi[Kooo] += (LapPhi - BareM_sqr_long * scale_factor * scale_factor * Phi[Kooo] * pow((1. + locPhiSqr / BareM_sqr_long / scale_factor / scale_factor), (p_mono - 1)))     *   (tau_conformal - Last_Momenta_update_conformal);
	        else
		        Pi[Kooo] += (LapPhi - BareM_sqr_trans * scale_factor * scale_factor * Phi[Kooo] * pow((1. + locPhiSqr / BareM_sqr_long / scale_factor / scale_factor), (p_mono - 1)))    *   (tau_conformal - Last_Momenta_update_conformal);
#endif
#endif
	  Kooo++; Kuoo++; Kouo++; Koou++; Kdoo++; Kodo++; Kood++;
	} 


      }
      
  Last_Momenta_update=tau;
#if (expansion==1)
  Last_Momenta_update_conformal=tau_conformal;
#endif
}

void updateFields()
{
  using namespace LocalLattice;
  
#if (expansion==0)
  for (int Kooo=pos(0,0,1); Kooo < Nc*sqr(N_t)*(N_eta_loc+1); Kooo++ )
    Phi[Kooo]+=Pi[Kooo]*(tau-Last_Field_update);
#endif
#if (expansion==1)
  for (int Kooo=pos(0,0,1); Kooo < Nc*sqr(N_t)*(N_eta_loc+1); Kooo++ )
    Phi[Kooo]+=Pi[Kooo]*(tau_conformal-Last_Field_update_conformal);
#endif  
   
  Last_Field_update=tau;
#if (expansion==1)
  Last_Field_update_conformal=tau_conformal;
#endif
}

void evolve()
{
  using namespace LocalLattice;
   
  if (regime==0)
  {
    tau+=0.5*dtau;
#if (expansion==1)
    tau_conformal = ( sqrt(1.+2.*H_0*tau) - 1. )/H_0;
    scale_factor = sqrt(1.+2.*H_0*tau);
#endif
    updateMomenta();
    tau+=0.5*dtau;
#if (expansion==1)
    tau_conformal = ( sqrt(1.+2.*H_0*tau) - 1. )/H_0;
    scale_factor = sqrt(1.+2.*H_0*tau);
#endif
    updateFields(); 
    
    normal_stepcount++;    
  }
 /* if (regime==1)
  {
    pseudo_double E0=analysis::E;
    
    tau+=0.5*dtau;
    updateMomenta();
    tau+=0.5*dtau;
    updateFields(); 
        
    adaptive_stepcount++;
    
    Calculate_Obs();
    pseudo_double dE=fabs(E0-analysis::E);
    if ( (dE/E0)>=1e-5 )
      dtau/=2.;   
  }*/
  
  Communication::exchange();
  MPI_Barrier(MPI_COMM_WORLD);
  total_stepcount++; 
}
  
#endif  
