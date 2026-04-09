#include <cmath>

#include "parameters.cpp"
#include "lattice.cpp"
#include "dynamics.cpp"

std::ofstream Condensation;

void InitBCounts();

namespace MomentumSpace
{
  std::complex<pseudo_double> *PhiP, *PiP;
#if (withEnergyPS==1)
  std::complex<pseudo_double> *DeltaP;
#endif
  
  pseudo_double *CorrF,*TotalCorrF;
  pseudo_double *CorrK,*TotalCorrK;
  pseudo_double *CorrQ,*TotalCorrQ;
#if (withEnergyPS==1)
  pseudo_double *CorrD,*TotalCorrD;
#endif
  
  int *BCounts,*TotalBCounts;
  
  pseudo_double pTmin,pTmax,binsizeT;
  int nTBins;
  
  pseudo_double *FAbsP;
  
  void init()
  {
    PhiP = new std::complex<pseudo_double>[(N_t/2+1)*N_t*LocalLattice::N_eta_loc];
    PiP = new std::complex<pseudo_double>[(N_t/2+1)*N_t*LocalLattice::N_eta_loc];
#if (withEnergyPS==1)
    DeltaP = new std::complex<pseudo_double>[(N_t/2+1)*N_t*LocalLattice::N_eta_loc];
#endif

    pTmin=(2.0/(a_t)) * sin(acos(-1.)/N_t);
    pTmax=sqrt( 3.0 * ( 4.0/sqr(a_t) ) );
    binsizeT=pTmin;
#if (fine_bin==1)
	binsizeT = pTmin/8.;
#endif

    nTBins = int(pTmax/binsizeT+0.5)+2;	// +1 would be sufficient but +2 should avoid any allocation errors

    BCounts=new int[nTBins];
    TotalBCounts = new int [nTBins];

    InitBCounts();
    
    CorrF = new  pseudo_double [Nc*nTBins];
    CorrK = new  pseudo_double [Nc*nTBins];
    CorrQ = new  pseudo_double [Nc*nTBins];
    
    TotalCorrF = new  pseudo_double [Nc*nTBins];
    TotalCorrK = new  pseudo_double [Nc*nTBins];
    TotalCorrQ = new  pseudo_double [Nc*nTBins];
#if (withEnergyPS==1)
    CorrD = new pseudo_double [Nc*nTBins];
    TotalCorrD = new pseudo_double [Nc*nTBins];
#endif
  }
  
  void fin()
  {
    delete[] PhiP;
    delete[] PiP;
    delete[] BCounts;
    delete[] CorrF;
    delete[] CorrK;
    delete[] CorrQ;
    delete[] TotalBCounts;
    delete[] TotalCorrF;
    delete[] TotalCorrK;
    delete[] TotalCorrQ;
#if (withEnergyPS==1)
    delete[] CorrD;
    delete[] TotalCorrD;
    delete[] DeltaP;
#endif
  }
}

using namespace MomentumSpace;

void InitBCounts()
{
  pseudo_double ppx=0., ppy=0., nu=0., pt=0., psqr=0.;
  int pzGlobal=0,posArrayT=0;

  for (int k=0; k<nTBins; k++)
    BCounts[k] = 0;

  for(int pz=0;pz<LocalLattice::N_eta_loc;pz++)
  {
    pzGlobal=pz+ID*LocalLattice::N_eta_loc;

//setting nu
    if (pzGlobal==0) 
      nu=0;
    else
    {	//nu!=0
      if (pzGlobal<N_eta/2+1) 
	nu=2.0/(a_eta)*sin((pzGlobal)*acos(-1.)/N_eta);
      else 
	nu=-2.0/(a_eta)*sin((pzGlobal)*acos(-1.)/N_eta);
    }
      
    for(int px=0;px<N_t;px++)
    {
//setting ppx
      if (px<N_t/2+1) 
	ppx=2.0/(a_t)*sin(px*acos(-1.)/N_t);
      else 
	ppx=-2.0/(a_t)*sin(px*acos(-1.)/N_t);
	
      for(int py=0;py<N_t/2+1; py++)
      {
//setting ppy
	ppy=2.0/(a_t)*sin(py*acos(-1.)/N_t);
	
	psqr = sqr(ppx)+sqr(ppy)+sqr(nu);
	posArrayT = int( sqrt(psqr)/binsizeT + 0.5 );

	if (posArrayT >= nTBins) 
	  std::cerr << "FATAL ERROR: bin array too small" << std::endl;
		
	BCounts[posArrayT]++;
	if (py!=0 && py!=N_t/2)
	  BCounts[posArrayT]++;
	 
      }
    }
  }
}

void Calculate_CorrelationFunctions()
{
  int posLOC=0, posARRAY=0, numCounts=0, pzGlobal=0;
  pseudo_double ppx,ppy,nu,psqr;

  for (int k = 0; k<nTBins*Nc; k++)
  {
    CorrF[k]=0.;
    CorrK[k]=0.;
    CorrQ[k]=0.;
#if (withEnergyPS==1)
    CorrD[k]=0.;
#endif
  }

// go through all colors
  for (int a = 0; a < Nc ; a++)
  {
    for (int k=0;k<(N_t/2+1)*N_t*LocalLattice::N_eta_loc;k++)
    {
      PhiP[k]=0.; PiP[k]=0.;
#if (withEnergyPS==1)
      DeltaP[k]=0.;
#endif 
    }
#if (withEnergyPS==1)
  pseudo_double avden=0., avden2=0., Avden=0., Avden2=0., aavden=0., Aavden=0.;
  //calculate the density contrast
  for (int Kooo=pos(0,0,1); Kooo < Nc*sqr(N_t)*(LocalLattice::N_eta_loc+1); Kooo++ )
    aavden+=LocalLattice::Delta[Kooo];
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&aavden, &Aavden, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&Aavden, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  Aavden/= ((pseudo_double) (N_t*N_t*N_t));
  for (int Kooo=pos(0,0,1); Kooo < Nc*sqr(N_t)*(LocalLattice::N_eta_loc+1); Kooo++ )
  {
    LocalLattice::Delta[Kooo]-= Aavden;
    LocalLattice::Delta[Kooo]/= Aavden;
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

//FT
    MPI_Barrier(MPI_COMM_WORLD);
    FFT::fftXtoP_direct(a,LocalLattice::Phi,PhiP);
    MPI_Barrier(MPI_COMM_WORLD);
    FFT::fftXtoP_direct(a,LocalLattice::Pi,PiP);
    MPI_Barrier(MPI_COMM_WORLD);
#if (withEnergyPS==1)
    FFT::fftXtoP_direct(a,LocalLattice::Delta,DeltaP);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
//sum to correlation functions
    for(int pz=0;pz<LocalLattice::N_eta_loc;pz++)
    {
	pzGlobal=pz+ID*LocalLattice::N_eta_loc;
	if (pzGlobal==0) 
	  nu=0;
	else
	{
	  if (pzGlobal<N_eta/2+1) 
	    nu=2.0/(a_eta)*sin((pzGlobal)*acos(-1.)/N_eta);
	  else 
	    nu=-2.0/(a_eta)*sin((pzGlobal)*acos(-1.)/N_eta);
	}

	for(int px=0;px<N_t;px++)
	{
	  if (px<N_t/2+1) 
	    ppx=2.0/(a_t)*sin(px*acos(-1.)/N_t);
	  else 
	    ppx=-2.0/(a_t)*sin(px*acos(-1.)/N_t);
		
	  for(int py=0;py<N_t/2+1; py++)
	  {
	    ppy=2.0/(a_t)*sin(py*acos(-1.)/N_t);
	
	    psqr = sqr(ppx)+sqr(ppy)+sqr(nu);
//position of fields in momentum space
	    posLOC=py+px*(N_t/2+1)+pz*N_t*(N_t/2+1);
	
//	    if (posLOC==0 && ID==0)
//	       	      std::cout<<(PhiP[posLOC])<<' '<<(PiP[posLOC])<<std::endl;

//position in the binned correlation functions
	    posARRAY = int(sqrt(psqr)/binsizeT+0.5)*Nc + a;
	    
	    CorrF[posARRAY]+=sqr(abs(PhiP[posLOC]));
	    CorrK[posARRAY]+=sqr(abs(PiP[posLOC]));
	    CorrQ[posARRAY]+=0.5*real(PhiP[posLOC]*conj(PiP[posLOC]) + conj(PhiP[posLOC])*PiP[posLOC]);
#if (withEnergyPS==1)
	    CorrD[posARRAY]+=sqr(abs(DeltaP[posLOC]));
#endif	   
	    if (py!=0 && py!=N_t/2)
	    {
	      CorrF[posARRAY]+=sqr(abs(PhiP[posLOC]));
	      CorrK[posARRAY]+=sqr(abs(PiP[posLOC]));
	      CorrQ[posARRAY]+=0.5*real(PhiP[posLOC]*conj(PiP[posLOC]) + conj(PhiP[posLOC])*PiP[posLOC]);
#if (withEnergyPS==1)
	      CorrD[posARRAY]+=sqr(abs(DeltaP[posLOC]));
#endif
	    }
	  }
	}
      }
    }
    
    //reset total arrays
    for (int k = 0; k<nTBins; k++)
    {
      TotalBCounts[k] = 0;

      for (int a = 0; a<Nc; a++)
      {
	TotalCorrF[k*Nc + a]=0.;
	TotalCorrK[k*Nc + a]=0.;
	TotalCorrQ[k*Nc + a]=0.;
#if (withEnergyPS==1)
	TotalCorrD[k*Nc + a]=0.;
#endif
      }
    }
        
    //sum up all the spectra
    MPI_Reduce(CorrF, TotalCorrF, Nc*nTBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(CorrK, TotalCorrK, Nc*nTBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(CorrQ, TotalCorrQ, Nc*nTBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#if (withEnergyPS==1)
    MPI_Reduce(CorrD, TotalCorrD, Nc*nTBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif  
        
    MPI_Reduce(BCounts, TotalBCounts, nTBins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
    MPI_Barrier(MPI_COMM_WORLD);
  
}

void Print_fp()
{
  using namespace LocalLattice;
  
  if (numeric::total_stepcount>0 && tau<numeric::LastMOoutput+TMO)
    return;
    
  if (numeric::total_stepcount==0)
  {
    numeric::LastMOoutput=tau;
#if (expansion==1)
    numeric::LastMOoutput_conformal=( sqrt(1.+2.*H_0*tau) - 1. )/H_0;
#endif
  }
  
  else 
    if (tau>=numeric::LastMOoutput+TMO)
    {
      numeric::LastMOoutput+=TMO;
#if (expansion==1)
      numeric::LastMOoutput_conformal=( sqrt(1.+2.*H_0*numeric::LastMOoutput) - 1. )/H_0;   
#endif
    }
  
  
  Calculate_CorrelationFunctions();
  
  if (ID==0)
  {
    std::ofstream specOut;
    char *fname=new char[256];
    snprintf(fname,256,"./%i/Distrib-%i-%g.txt",LASTRUN,  runID, really_double(numeric::LastMOoutput));
    specOut.open(fname);
    delete[] fname;
      
//      modesout<<tauLOC<<' ';    
    
    pseudo_double pabs=0.;
    pseudo_double NormCorr=1./(sqr(N_t*a_t)*N_eta*a_eta);	//FT factor
    pseudo_double norm2=0.;		//number of items in the bin

    for (int k = 0; k < nTBins; k++)
    {
//divide by number of total counts
      norm2 = (TotalBCounts[k]==0 ? 0 : 1/(pseudo_double)(TotalBCounts[k]));
      pabs = binsizeT * k;
	
      if (norm2!=0 && TotalCorrF[k*Nc]!=0)
      {
	if (k!=0)
	  specOut << std::endl;;
	specOut << pabs << " ";

	specOut << TotalBCounts[k] << " ";
	
	pseudo_double FF=0,KK=0,QQ=0;
#if (withEnergyPS==1)
	pseudo_double DD=0;
#endif
	for (int a=0;a<Nc;a++)
	{
//	  specOut << TotalCorrF[k*Nc+a]*NormCorr*norm2 << " ";
	  FF+=TotalCorrF[k*Nc+a];
//	  specOut << TotalCorrK[k*Nc+a]*NormCorr*norm2 << " ";
	  KK+=TotalCorrK[k*Nc+a];
//	  specOut << TotalCorrQ[k*Nc+a]*NormCorr*norm2 << " ";
	  QQ+=TotalCorrQ[k*Nc+a];
#if (withEnergyPS==1)
	  DD+=TotalCorrD[k*Nc+a];
#endif
	}
	FF/=Nc;
	KK/=Nc;
	QQ/=Nc;
#if (withEnergyPS==1)
	DD/=Nc;
#endif
	
	pseudo_double omg=sqrt(KK/FF);
	specOut << omg<<' ';
	specOut <<FF*NormCorr*norm2<<' '<<KK*NormCorr*norm2<<' '<<sqrt(FF*KK)*NormCorr*norm2 << " ";
//	specOut << sqrt(FF*KK - QQ*QQ)*NormCorr*norm2 << " "<< (1./(2.*omg))*(FF*omg*omg+KK)*NormCorr*norm2 <<' '<<sqrt(FF*KK)*NormCorr*norm2 << " ";
//	specOut <<FF*NormCorr*norm2<<' '<<KK*NormCorr*norm2<<' '<<QQ*NormCorr*norm2<<' ';
#if (expansion==1)
	specOut<< scale_factor << " ";
#endif
#if (withEnergyPS ==1)  
	specOut <<DD*NormCorr*norm2<<' ';
#endif  
	specOut <<FF*NormCorr*norm2*pabs*pabs*pabs / 2./acos(-1.)/acos(-1.)<<' ';
#if (withEnergyPS ==1) 
	specOut <<DD*NormCorr*norm2*pabs*pabs*pabs / 2./acos(-1.)/acos(-1.)<<' ';
#endif
	 
	if (FF*KK - QQ*QQ<0)
	  std::cout<<"negative at "<<k<<std::endl;
	    
      }
    }
    specOut.close();
 //     modesout<<std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

