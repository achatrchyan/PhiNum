#ifndef _FIELDANALYSISCPP_
#define _FIELDANALYSISCPP_

#include <iostream>
#include <fstream> 
#include <iomanip>
#include <mpi.h>

#include "parameters.cpp"
#include "lattice.cpp"
#include "dynamics.cpp"

std::ofstream SpatialObs;   //this is to output the spatial observables



//A useful structure, spatial point on the grid
struct point
{
  int pos;
  int id;
  int z()
  {
    return id*LocalLattice::N_eta_loc-1+pos/(Nc*N_t*N_t);
  }
  int y()
  {
    return ( pos%(Nc*N_t*N_t) )/(Nc*N_t);
  }
  int x()
  {
    return ( ( pos%(Nc*N_t*N_t) ) % (Nc*N_t) )/Nc;
  }
};

point quick_point()
{
   point P;
   P.id =0 ;
   P.pos = (N_t/2)*(N_t/2)*(N_t/2);
   return P;
}



//Determine which folder to create, to avoid overwriting
void Prepare_Output()
{
  if (ID==0)
  {
    std::ifstream fin ("lastrun.txt");
    fin>>LASTRUN;
    fin.close();
    std::ofstream fout ("lastrun.txt");
    fout<<++LASTRUN<<endl;
    fout.close();
    
    char* cmd=new char[256];
    snprintf(cmd,256,"mkdir ./%i",LASTRUN);
    system(cmd);
    delete[] cmd;
    
 //   std::cerr << "Created a folder "<< LASTRUN<<std::endl; 
  }  
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&LASTRUN, 1, MPI_INT, 0, MPI_COMM_WORLD); 
}

void Print_Info_Console()
{
  if (ID==0)
  {
    std::cout << "Parameters" <<std::endl;
   
#if (theory==0)  
    std::cout << "Nc=" <<  Nc <<std::endl; 
    std::cout << "m^2=" <<  Mass_Mode*M_sqr <<std::endl; 
    std::cout << "Lambda=" <<  Lambda_Mode*Lambda <<std::endl; 
    std::cout << "Nyu=" <<  (Lambda/Kappa)*(Lambda/Kappa) <<std::endl;
#endif

#if (theory==1)  
    std::cout << "Nc=" <<  Nc <<std::endl; 
    std::cout << "m^2=" <<  Mass_Mode*M_sqr <<std::endl; 
    std::cout << "Lambda=" <<  Lambda_Mode*Lambda <<std::endl; 
    std::cout << "h=" <<  h/sqrt(Lambda) <<std::endl; 
#endif
    
#if (theory==2)  
    std::cout << "Nc=" <<  Nc <<std::endl; 
    std::cout << "m^2=" <<  M_sqr <<std::endl; 
    std::cout << "f=" <<  F_SSB <<std::endl; 
    std::cout << "Lambda=" <<  Lambda <<std::endl; 
#endif
    
#if (theory==3)  
    std::cout << "Nc=" <<  Nc <<std::endl; 
    std::cout << "m_mono^2=" <<  M_mono_sqr <<std::endl; 
    std::cout << "m^2=" <<  M_sqr <<std::endl; 
    std::cout << "f=" <<  F_SSB <<std::endl; 
    std::cout << "alpha=" <<  alpha <<std::endl; 
    std::cout << "Lambda=" <<  Lambda <<std::endl; 
#endif

#if (theory==4)  
	std::cout << "Nc=" << Nc << std::endl;
	std::cout << "m^2=" << Mass_Mode*M_sqr << std::endl;
	std::cout << "Lambda=" << Lambda_Mode*Lambda << std::endl;
	std::cout << "g=" << g << std::endl;
#endif
    
#if (theory==5)  
    std::cout << "Nc=" <<  Nc <<std::endl; 
    std::cout << "m^2=" <<  M_sqr <<std::endl; 
    std::cout << "f=" <<  F_SSB <<std::endl; 
    std::cout << "p=" <<  p_mono <<std::endl; 
    std::cout << "Lambda=" <<  Lambda <<std::endl; 
#endif
    std::cout << "N_s=" << N_t << std::endl;
    std::cout << "a_s=" << a_t << std::endl;
    std::cout << "T0=" << T0 << std::endl;
    std::cout << "TMax=" << TMax << std::endl;
    std::cout << "dt=" << a_t*dtBaseTrans << std::endl;
    std::cout << "TSO=" << TSO << std::endl;
    std::cout << "TMO=" << TMO << std::endl;

#if (initialize==2)
    std::cout << "f0(p)="<<Amplitude<<"/Lambda, for p<"<<QS<<std::endl;
    std::cout << "phi0=" << phi0<<std::endl;
    std::cout << "pi0=" << pi0<<std::endl;    
#endif    
    
#if (initialize==3)
    std::cout << "f0(p)="<<Amplitude<<"/Lambda, for p<"<<QS<<std::endl;
    std::cout << "phi0_0=" << phi0_0<<std::endl;
    std::cout << "phi0_1=" << phi0_1<<std::endl;
    std::cout << "pi0_0=" << pi0_0<<std::endl;
    std::cout << "pi0_1=" << pi0_1<<std::endl;        
#endif    
    
    std::cout<< std::endl<< runs;
    (runs==1)?(std::cout<<" run with"):(std::cout<<" runs with");
    (renorm_tadpole==1)?(std::cout<<" "):(std::cout<<"out ");
    std::cout<<"renormalization."<< std::endl;
  }
}

void Print_Info_File()
{
  if (ID==0)
  {
    std::ofstream out;
	
    char *fname=new char[256];
    snprintf(fname,256,"./%i/info.txt",LASTRUN);
    out.open(fname);
    delete[] fname;
  
    out << "Parameters" <<std::endl;
#if (theory==0)  
    out << "Nc=" <<  Nc <<std::endl; 
    out << "m^2=" <<  Mass_Mode*M_sqr <<std::endl; 
    out << "Lambda=" <<  Lambda_Mode*Lambda <<std::endl; 
    out << "Nyu=" <<  (Lambda/Kappa)*(Lambda/Kappa) <<std::endl; 
#endif

#if (theory==1)  
    out << "Nc=" <<  Nc <<std::endl; 
    out << "m^2=" <<  Mass_Mode*M_sqr <<std::endl; 
    out << "Lambda=" <<  Lambda_Mode*Lambda <<std::endl; 
    out << "h=" << h/sqrt(Lambda) <<std::endl; 
#endif
    
#if (theory==2)  
    out << "Nc=" <<  Nc <<std::endl; 
    out << "m^2=" <<  M_sqr <<std::endl; 
    out << "f=" <<  F_SSB <<std::endl; 
    out << "Lambda=" <<  Lambda <<std::endl; 
#endif
    
#if (theory==3)  
    out << "Nc=" <<  Nc <<std::endl; 
    out << "m_mono^2=" <<  M_mono_sqr <<std::endl; 
    out << "m^2=" <<  M_sqr <<std::endl; 
    out << "f=" <<  F_SSB <<std::endl; 
    out << "alpha=" <<  alpha <<std::endl; 
    out << "Lambda=" <<  Lambda <<std::endl; 
#endif

#if (theory==4)  
	out << "Nc=" << Nc << std::endl;
	out << "m^2=" << Mass_Mode*M_sqr << std::endl;
	out << "Lambda=" << Lambda_Mode*Lambda << std::endl;
	out << "g=" << g << std::endl;
#endif
    
#if (theory==5)  
    out << "Nc=" <<  Nc <<std::endl; 
    out << "m^2=" <<  M_sqr <<std::endl; 
    out << "f=" <<  F_SSB <<std::endl; 
    out << "p=" << p_mono <<std::endl;
    out << "Lambda=" <<  Lambda <<std::endl; 
#endif
    out << "N_s=" << N_t << std::endl;
    out << "a_s=" << a_t << std::endl;
    out << "T0=" << T0 << std::endl;
    out << "TMax=" << TMax << std::endl;
    out << "dtBT=" << dtBaseTrans << std::endl;
    out << "dtBL=" << dtBaseLong<< std::endl;
    out << "TSO=" << TSO << std::endl;
    out << "TMO=" << TMO << std::endl;

#if (initialize==2)
    out << "f0(p)="<<Amplitude<<"/Lambda, for p<"<<QS<<std::endl;
    out << "phi0=" << phi0<<std::endl;
    out << "pi0=" << pi0<<std::endl;    
#endif    
    
#if (initialize==3)
    out << "f0(p)="<<Amplitude<<"/Lambda, for p<"<<QS<<std::endl;
    out << "phi0_0=" << phi0_0<<std::endl;
    out << "phi0_1=" << phi0_1<<std::endl;
    out << "pi0_0=" << pi0_0<<std::endl;    
    out << "pi0_1=" << pi0_1<<std::endl;    
#endif    
   
  
    out<< std::endl<< runs;
    (runs==1)?(out<<" run with"):(out<<" runs with");
    (renorm_tadpole==1)?(out<<" "):(out<<"out ");
    out<<"renormalization."<< std::endl;
  
    out.close();
  }
}




//Make a copy of the field configuration
void SaveConfig(int Name)
{
  FILE * of;
  
  char *fname=new char[256];
  snprintf(fname,256,"./field-%i-%i.out",Name,ID);
  of = fopen (fname, "w");
  delete[] fname;
  
  for (int Kooo=pos(0,0,1); Kooo < Nc*sqr(N_t)*(LocalLattice::N_eta_loc+1); Kooo++ )
  {
#if (longdouble == 0)      
    fprintf (of,"%e %e\n", LocalLattice::Phi[Kooo], LocalLattice::Pi[Kooo]);
#elif (longdouble == 1)
    fprintf (of,"%Le %Le\n", LocalLattice::Phi[Kooo], LocalLattice::Pi[Kooo]);
#endif
  }
  fclose(of);

  if (ID==0)
  {
    fname=new char[256];
    snprintf(fname,256,"./moment-%i.out",Name);
    ofstream fout (fname);
    delete[] fname;
    fout<<LocalLattice::tau<<endl;
  }
}

void ExtractConfig(int Name)
{
  using namespace LocalLattice;
  FILE * of;
  
  char *fname=new char[256];
  snprintf(fname,256,"./field-%i-%i.out",Name,ID);
  of = fopen (fname, "r");
  delete[] fname;
  
  for (int Kooo=pos(0,0,1); Kooo < Nc*sqr(N_t)*(N_eta_loc+1); Kooo++ )
  {
#if (longdouble == 0)      
	  fscanf(of, "%lf %lf\n", Phi + Kooo, Pi + Kooo);
#elif (longdouble == 1)
    fscanf (of,"%Le %Le\n", Phi+Kooo, Pi+Kooo);
#endif
  }
  fclose(of);  
  Communication::exchange();
  MPI_Barrier(MPI_COMM_WORLD);

  fname=new char[256];
  snprintf(fname,256,"./moment-%i.out",Name);
  ifstream fin (fname);
  delete[] fname;
  fin>>T0;
}



//------------------------------------------------------------------------------------------------------------------------

void CalculateExtrema()
{
  using namespace LocalLattice;
  
//reset
  maxPhi=-10000000.0, minPhi=10000000.0, maxAbsPhi=0.0;
  maxpos=-1, minpos=-1, maxAbspos=-1;
  GlmaxPhi=0.0, GlminPhi=0.0, GlmaxAbsPhi=0.0;
  GlmaxID=-1, GlminID=-1, GlmaxAbsID=-1;  
  
//locally
  pseudo_double AbsPhi=0.;
  for (int k=pos(0,0,1); k < Nc*sqr(N_t)*(N_eta_loc+1); k++ )
  {
    AbsPhi+=((Phi[k])*(Phi[k]));
    if ((k+1)%Nc==0)
    {
      AbsPhi=sqrt(AbsPhi);
      if (AbsPhi>maxAbsPhi)
      {
      	maxAbsPhi=AbsPhi;
	      maxAbspos=k+1-Nc;	    
      }
      AbsPhi=0;
    }    
    
    if (Phi[k]>maxPhi)
    {
      maxPhi=Phi[k];
      maxpos=k;	    
    }
    if (Phi[k]<minPhi)
    {
      minPhi=Phi[k];
      minpos=k;	    
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Reduce(&maxPhi, &GlmaxPhi, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&minPhi, &GlminPhi, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxAbsPhi, &GlmaxAbsPhi, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Bcast(&GlmaxPhi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GlminPhi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GlmaxAbsPhi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);		
  
  int lcminID=-1, lcmaxID=-1, lcmaxAbsID=-1;
  if (minPhi==GlminPhi)
    lcminID=ID;
  if (maxPhi==GlmaxPhi)
    lcmaxID=ID;
  if (maxAbsPhi==GlmaxAbsPhi)
    lcmaxAbsID=ID;
   
  MPI_Barrier(MPI_COMM_WORLD);		
  
  MPI_Reduce(&lcmaxID, &GlmaxID, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&lcminID, &GlminID, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&lcmaxAbsID, &GlmaxAbsID, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);		
    
  MPI_Bcast(&GlmaxID, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GlminID, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&GlmaxAbsID, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
  MPI_Bcast(&minpos, 1, MPI_INT, GlminID, MPI_COMM_WORLD);
  MPI_Bcast(&maxpos, 1, MPI_INT, GlmaxID, MPI_COMM_WORLD);
  MPI_Bcast(&maxAbspos, 1, MPI_INT, GlmaxAbsID, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  return ;
}

void CalculateVolumeAverages()
{
  using namespace LocalLattice;


//calculate on each node separately  
  for (int a=0;a<Nc;a++)
  {
    phi[a]=0.0;
    pi[a]=0.0;
  }
  phi2=0.; pi2=0.; 
  phi4=0.; phi6=0.;
  phipi=0.;
  
  Interaction = 0.;
  Charge=0.;

  pseudo_double GradPhiSqr=0.0;

  int Kooo=pos(0,0,1);
  
  for(int z=1;z<N_eta_loc+1;z++)
    for(int y=0;y<N_t;y++)
      for(int x=0;x<N_t;x++)
      {
      	int Kuoo=pos(mod(x+1,N_t), y, z);
	      int Kouo=pos(x, mod(y+1,N_t), z);
	      int Koou=pos(x, y, z+1);

#if (withcharge==1)
        if (Nc==2)
	        Charge+=Phi[Kooo]*Pi[Kooo+1]-Phi[Kooo+1]*Pi[Kooo];
#endif	    
	      pseudo_double locPhiSqr=0.0;
        pseudo_double locPhi = 0.;  
	      for(int a=0;a<Nc;a++)
	      {
	        phi[a]+=Phi[Kooo];	pi[a]+=Pi[Kooo];
	  
	        phi2+=sqr(Phi[Kooo]);	pi2+=sqr(Pi[Kooo]);
	        phipi+=Phi[Kooo]*Pi[Kooo];
	  
	        locPhiSqr+=sqr(Phi[Kooo]);	//it is calculated separately for phi^4 and phi^6
          locPhi+=Phi[Kooo];
	  
	        GradPhiSqr+= ( sqr( Phi[Kuoo]-Phi[Kooo] )/sqr(a_t) + sqr( Phi[Kouo]-Phi[Kooo] )/sqr(a_t) + sqr( Phi[Koou]-Phi[Kooo] )/sqr(a_t) );
#if (withEnergyPS==1)
          Delta[Kooo] = (0.5*sqr(Pi[Kooo]) +  sqr( Phi[Kuoo]-Phi[Kooo] )/sqr(a_t) + sqr( Phi[Kouo]-Phi[Kooo] )/sqr(a_t) + sqr( Phi[Koou]-Phi[Kooo] )/sqr(a_t) );
#if (expansion==1)
          Delta[Kooo] += (0.5*sqr(H_0 / scale_factor)*sqr(Phi[Kooo]) - Phi[Kooo] * Pi[Kooo] * H_0 / scale_factor);
#endif
#endif

//mass term
#if (theory==0)
#if (expansion==0)
          if (a==0)
	          Interaction+= (Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo]);
	        else
	          Interaction+=(Mass_Mode)*0.5*BareM_sqr_trans*sqr(Phi[Kooo]);	 
#endif
#if (expansion==1)
          if (a==0)
	          Interaction+= (Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
          else
            Interaction+=(Mass_Mode)*0.5*BareM_sqr_trans*scale_factor*scale_factor*sqr(Phi[Kooo]);
#endif
#endif	  
#if (theory==1)
#if (expansion==0)
	        if (a==0)
	          Interaction+=( h*Phi[Kooo] + (Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo]) );
	        else
	          Interaction+=(Mass_Mode)*0.5*BareM_sqr_trans*sqr(Phi[Kooo]);	  
#endif
#if (expansion==1)
          if (a==0)
            Interaction+= ( h*Phi[Kooo]*scale_factor*scale_factor*scale_factor + (Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]) );
	        else
            Interaction+=(Mass_Mode)*0.5*BareM_sqr_trans*scale_factor*scale_factor*sqr(Phi[Kooo]);
#endif
#endif
#if (theory==3)
#if (expansion==0)
	        if (a==0)
	          Interaction+=0.5*BareM_mono_sqr_long*sqr(Phi[Kooo]);
	        else
	          Interaction+=0.5*BareM_mono_sqr_trans*sqr(Phi[Kooo]);
#endif
#if (expansion==1)
	        if (a==0)
	          Interaction+=(0.5*BareM_mono_sqr_long*( sqr(scale_factor)  )*sqr(Phi[Kooo])   );
	        else
	          Interaction+=(0.5*BareM_mono_sqr_trans*( sqr(scale_factor)  )*sqr(Phi[Kooo])      );
#endif
#endif
#if (theory==4)
#if (expansion==0)
	        if (a==0)
	          Interaction+=(Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo]);
	        else
	          Interaction+=(Mass_Mode)*0.5*BareM_sqr_trans*sqr(Phi[Kooo]);	  
#endif
#if (expansion==1)
          if (a==0)
            Interaction+=(Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
	        else
            Interaction+=(Mass_Mode)*0.5*BareM_sqr_trans*scale_factor*scale_factor*sqr(Phi[Kooo]);
#endif
#endif

#if (withEnergyPS==1)
#if (expansion == 0)
#if (theory==0)
          Delta[Kooo] += 0.5*(Mass_Mode)*BareM_mono_sqr_long*sqr(Phi[Kooo]);
          Delta[Kooo] += (Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc)) +(sqr(sqr(Phi[Kooo]))*sqr(Phi[Kooo])/(720.0*Nc*Nc*Kappa*Kappa)) ;
#endif
#if (theory==1)
          Delta[Kooo] += 0.5*(Mass_Mode)*BareM_mono_sqr_long*sqr(Phi[Kooo]);
          Delta[Kooo] += (  h*Phi[Kooo] + (Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc)));
#endif
#if (theory==2)
          Delta[Kooo]+=sqr(BareM_sqr_long)*(1.-cos(  Phi[Kooo] /(sqrt(BareM_sqr_long)) ));
#endif
#if (theory==3)
          Delta[Kooo] +=0.5*BareM_mono_sqr_long*sqr(Phi[Kooo]);
          Delta[Kooo] += sqr(BareM_sqr_long)*(1.-cos(  Phi[Kooo] /(sqrt(BareM_sqr_long)) ));
#endif
#if (theory==4)
          Delta[Kooo] += (Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo]);
          Delta[Kooo] += (((-g*sqr(Phi[Kooo])*Phi[Kooo]) / (6.0)) + (Lambda_Mode)*(sqr(sqr(Phi[Kooo])) / (24.0*Nc)));
#endif
#if (theory==5)
          Delta[Kooo] += sqr(BareM_sqr_long)*0.5 / p_mono * (pow((1. + Phi[Kooo] * Phi[Kooo] / BareM_sqr_long), p_mono) - 1.);
#endif
#endif
#if (expansion == 1)
#if (theory==0)
          Delta[Kooo] += (Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
          Delta[Kooo] += (Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc))+(sqr(sqr(Phi[Kooo]))*sqr(Phi[Kooo])/(720.0*Nc*Nc*Kappa*Kappa*scale_factor*scale_factor));
#endif
#if (theory==1)
          Delta[Kooo] += (Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
          Delta[Kooo] += (  h*Phi[Kooo]*scale_factor*scale_factor*scale_factor + (Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc)) );
#endif
#if (theory==2)
          Delta[Kooo]+=sqr(scale_factor*scale_factor*BareM_sqr_long)*(1.-cos(  Phi[Kooo] /(scale_factor*sqrt(BareM_sqr_long)) ));
#endif
#if (theory==3)
          Delta[Kooo] +=(0.5*BareM_mono_sqr_long*( sqr(scale_factor)  )*sqr(Phi[Kooo])   );
          Delta[Kooo]+=sqr(scale_factor*scale_factor*BareM_sqr_long)*(1.-cos(  Phi[Kooo] /(scale_factor*sqrt(BareM_sqr_long)) ));
#endif
#if (theory==4)
          Delta[Kooo] += (Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
          Delta[Kooo] += (((-g*scale_factor*sqr(Phi[Kooo])*Phi[Kooo]) / (6.0)) + (Lambda_Mode)*(sqr(sqr(Phi[Kooo])) / (24.0*Nc)));
#endif
#if (theory==5)
          Delta[Kooo] += sqr(scale_factor*scale_factor*BareM_sqr_long)*0.5/p_mono * (pow((1. + Phi[Kooo] * Phi[Kooo]/BareM_sqr_long/scale_factor/scale_factor), p_mono) - 1.); 
#endif
#endif
#endif
	        Kooo++;
	        Kuoo++; Kouo++; Koou++;
	      }

	      phi4+=sqr(locPhiSqr);
	      phi6+=sqr(locPhiSqr)*locPhiSqr;

//interaction terms
#if (theory==0)
#if (expansion==0)
	      Interaction+=(Lambda_Mode)*(sqr(locPhiSqr)/(24.0*Nc))+(sqr(locPhiSqr)*locPhiSqr/(720.0*Nc*Nc*Kappa*Kappa));
#endif
#if (expansion==1)
	      Interaction+=( (Lambda_Mode)*(sqr(locPhiSqr)/(24.0*Nc)) +(sqr(locPhiSqr)*locPhiSqr/(720.0*Nc*Nc*Kappa*Kappa*scale_factor*scale_factor)) );
#endif
#endif
#if (theory==1)
#if (expansion==0)
	      Interaction+=(Lambda_Mode)*(sqr(locPhiSqr)/(24.0*Nc));
#endif
#if (expansion==1)
	      Interaction+= (Lambda_Mode)*(sqr(locPhiSqr)/(24.0*Nc));
#endif
#endif
#if (theory==2)
#if (expansion==0)
	      Interaction+= sqr(BareM_sqr_long)*(1.-cos( sqrt(locPhiSqr/BareM_sqr_long) ));
#endif
#if (expansion==1)
	      Interaction+= sqr(scale_factor*scale_factor*BareM_sqr_long)*(1.-cos( sqrt(locPhiSqr)/(scale_factor*sqrt(BareM_sqr_long)) ));
#endif
#endif
#if (theory==3)
#if (expansion==0)
	      Interaction+= sqr(BareM_sqr_long)*(1.-cos( sqrt(locPhiSqr/BareM_sqr_long) ));
#endif
#if (expansion==1)
	      Interaction+= sqr(scale_factor*scale_factor*BareM_sqr_long)*(1.-cos( sqrt(locPhiSqr)/(scale_factor*sqrt(BareM_sqr_long)) ));
#endif
#endif
#if (theory==4)
#if (expansion==0)
	      Interaction += (((-g*locPhi*locPhiSqr) / (6.0)) + (Lambda_Mode)*(sqr(locPhiSqr) / (24.0*Nc)));
#endif
#if (expansion==1)
	      Interaction+=(((-g*scale_factor*locPhi*locPhiSqr) / (6.0)) + (Lambda_Mode)*(sqr(locPhiSqr) / (24.0*Nc)));
#endif
#endif
#if (theory==5)
#if (expansion==0)
	      Interaction+= sqr(BareM_sqr_long)*0.5/p_mono * (pow((1. + locPhiSqr/BareM_sqr_long), p_mono) - 1.);
#endif
#if (expansion==1)
	      Interaction+= sqr(scale_factor*scale_factor*BareM_sqr_long)*0.5/p_mono * (pow((1. + locPhiSqr/BareM_sqr_long/scale_factor/scale_factor), p_mono) - 1.); 
#endif
#endif
      }

  for (int a=0;a<Nc;a++)
  {
    phi[a]/=(vol);	pi[a]/=(vol);
  }

  Charge = (0.5*pi2 - GradPhiSqr/6. - Interaction);
  Energy = (0.5*pi2 + 0.5*GradPhiSqr + Interaction);
#if (expansion==1)
  Energy +=(0.5*phi2*sqr(H_0/scale_factor) - phipi*H_0/scale_factor );
  Charge +=(0.5*phi2*sqr(H_0/scale_factor) - phipi*H_0/scale_factor );
#endif

  Energy /= (vol);
  Interaction /= (vol);
  Charge/=(vol);

  phi2/=(vol);		pi2/=(vol);
  phi4/=(vol);
  phi6/=(vol);
  phipi/=(vol);

//average over all nodes via communication  
  avgEnergy=0.;
  avgInteraction=0.;
  avgpi2=0.;
  avgphi2=0.;
  avgphi4=0.;
  avgphi6=0.;
  avgphipi=0;
  avgCharge=0.;
  
  for (int a=0; a<Nc; a++)
  {
    avgphi[a]=0.;
    avgpi[a]=0.;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(&Energy, &avgEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Interaction, &avgInteraction, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Charge, &avgCharge, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(phi, avgphi, Nc, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(pi, avgpi, Nc, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&phi2, &avgphi2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&phi4, &avgphi4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&phi6, &avgphi6, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&pi2, &avgpi2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&phipi, &avgphipi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  	
  if (ID==0)
  {
    avgEnergy/=(nodes);
    avgInteraction/=(nodes);
    avgCharge/=(nodes);
    avgpi2/=(nodes);
    avgphi2/=(nodes);
    avgphi4/=(nodes);
    avgphi6/=(nodes);
    avgphipi/=(nodes);
    
    for (int a=0; a<Nc; a++)
    {
      avgphi[a]/=(nodes);
      avgpi[a]/=(nodes);
    }
  } 

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&avgEnergy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgInteraction, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgCharge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(avgphi, Nc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(avgpi, Nc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgphi2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgphi4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgphi6, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgpi2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&avgphipi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Print_Extrema()
{
  using namespace std;
  using namespace LocalLattice;
  
  if (numeric::total_stepcount>0 && tau<numeric::LastSOoutput+TSO)
    return;
    
  if (numeric::total_stepcount==0)
  {
    if (ID==0)
    {
      char *fname=new char[256];
      snprintf(fname,256,"./%i/SO-%i.txt",LASTRUN, runID);
      SpatialObs.open(fname, fstream::app);
      delete[] fname;    
    }
    numeric::LastSOoutput=tau;
#if (expansion==1)
    numeric::LastSOoutput_conformal=( sqrt(1.+2.*H_0*tau) - 1. )/H_0;
#endif
  }
  
  else 
    if (tau>numeric::LastSOoutput+TSO)
    {
      numeric::LastSOoutput+=TSO;
#if (expansion==1)
      numeric::LastSOoutput_conformal=( sqrt(1.+2.*H_0*numeric::LastSOoutput) - 1. )/H_0;
#endif
    }

  CalculateExtrema();
  
  if (ID==0)
  {
#if (expansion ==0)
  	SpatialObs<<numeric::LastSOoutput<<' ';				//time
#endif
#if (expansion==1)
  	SpatialObs<<scale_factor<<' '; //scale factor
#endif
    
    SpatialObs<<GlminPhi<<' '<<GlmaxPhi<<' '<<GlmaxAbsPhi<<endl;
  }
}

void Print_AveragesExtrema()
{
  using namespace std;
  using namespace LocalLattice;
  
  if (numeric::total_stepcount>0 && tau<numeric::LastSOoutput+TSO)
    return;
    
  if (numeric::total_stepcount==0)
  {
    if (ID==0)
    {
      if (runID>1)
        SpatialObs.close();
      char *fname=new char[256];
      snprintf(fname,256,"./%i/SO-%i.txt",LASTRUN, runID);
      SpatialObs.open(fname, fstream::app);
      delete[] fname;    
    }
    numeric::LastSOoutput=tau;
#if (expansion==1)
    numeric::LastSOoutput_conformal=( sqrt(1.+2.*H_0*tau) - 1. )/H_0;
#endif
  }
  
  else 
    if (tau>=numeric::LastSOoutput+TSO)
    {
      numeric::LastSOoutput+=TSO;
#if (expansion==1)
      numeric::LastSOoutput_conformal=( sqrt(1.+2.*H_0*numeric::LastSOoutput) - 1. )/H_0;
#endif
    }  
  CalculateExtrema();
  CalculateVolumeAverages();
  
  if (ID==0)
  {
#if (expansion ==0)
	  SpatialObs<<numeric::LastSOoutput<<' ';				//time
#endif
#if (expansion==1)
	  SpatialObs<<scale_factor<<' '; //scale factor
#endif
    for (int i=0;i<Nc;i++)
      SpatialObs<<avgphi[i]<<' '<<avgpi[i]<<' ';
    SpatialObs<<avgphi2<<' '<<avgpi2<<' '<<avgphipi<<' '<<avgphi4<<' '<<avgphi6<<' ';		//average field squared
    SpatialObs<<avgEnergy<<' '<<avgInteraction<<' '<<avgCharge<<' ';
    
    SpatialObs<<GlminPhi<<' '<<GlmaxPhi<<' '<<GlmaxAbsPhi<<endl;
  }
}

void Print_Averages()
{
  using namespace std;
  using namespace LocalLattice;
  
  if (numeric::total_stepcount>0 && tau<numeric::LastSOoutput+TSO)
    return;
    
  if (numeric::total_stepcount==0)
  {
    if (ID==0)
    {
      char *fname=new char[256];
      snprintf(fname,256,"./%i/SO-%i.txt",LASTRUN, runID);
      SpatialObs.open(fname, fstream::app);
      delete[] fname;    
    }
    numeric::LastSOoutput=tau;
#if (expansion==1)
    numeric::LastSOoutput_conformal=( sqrt(1.+2.*H_0*tau) - 1. )/H_0;
#endif
  }
  
  else 
    if (tau>=numeric::LastSOoutput+TSO)
    {
      numeric::LastSOoutput+=TSO;
#if (expansion==1)
      numeric::LastSOoutput_conformal=( sqrt(1.+2.*H_0*numeric::LastSOoutput) - 1. )/H_0;
#endif
    }  
  CalculateVolumeAverages();
  
  if (ID==0)
  {
#if (expansion ==0)
	  SpatialObs<<numeric::LastSOoutput<<' ';				//time
#endif
#if (expansion==1)
	  SpatialObs<<scale_factor<<' '; //scale factor
#endif
    for (int i=0;i<Nc;i++)
      SpatialObs<<avgphi[i]<<' '<<avgpi[i]<<' ';	//average field components
    SpatialObs<<avgphi2<<' '<<avgpi2<<' '<<avgphipi<<' '<<avgphi4<<' '<<avgphi6<<' ';		//average field squared
    SpatialObs<<avgEnergy<<' '<<avgInteraction<<' '<<avgCharge<<endl;
  }
}



//------------------------------------------------------------------------------------------------------------------------

void field_statistics()
{
	if (ID!=0)
		return;
	using namespace LocalLattice;

	int counter[2][500];
	for (int i = 0;i<500;i++)
		for (int j = 0;j<2;j++)
			counter[j][i] = 0;

	CalculateExtrema();
	pseudo_double dphi = 2.*GlmaxAbsPhi / 500;

	for (int K = pos(0, 0, 1); K < Nc*sqr(N_t)*(N_eta_loc + 1); K++)
		counter[0][int(fabs(Phi[K]) / dphi)]++;

	for (int i = 498;i >= 0;i--)
		counter[1][i] = counter[1][i + 1] + counter[0][i];

	ofstream fout;
	char *fname=new char[256];
	snprintf(fname,256,"./%i/field_statistics.txt",LASTRUN);
	fout.open(fname, std::fstream::app);
	delete[] fname;
	
	for (int i = 0;i<500;i++)
		fout << double(i)*dphi << ' ' << double(counter[0][i]) / dphi / double(N_t * N_t * N_eta_loc * Nc) << ' ' << double(counter[1][i]) / double(N_t * N_t * N_eta_loc * Nc) << endl;

	fout.close();

	return ;
}

bool false_vacuum_decayed()
{
#if (theory == 0 || theory==4 || theory == 1)
  return (LocalLattice::avgInteraction<0 || LocalLattice::avgInteraction!=LocalLattice::avgInteraction || LocalLattice::avgInteraction*0.!=0.);
#endif
  return 0;
}


#endif
