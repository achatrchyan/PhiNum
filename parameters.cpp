#ifndef _PARAMETERS_CPP_
#define _PARAMETERS_CPP_

#include "cfile.c"


/*=======Parameters=========*/

//Simulation parameters
int runs=1;
int runID;

#define doubleprecision 0

#if (doubleprecision == -1)
#define pseudo_double float
#define MPI_PSEUDO_DOUBLE MPI_FLOAT
#define really_double double

#define p_fftw_complex fftwf_complex
#define p_fftw_plan fftwf_plan
#define p_fftw_mpi_init fftwf_mpi_init
#define p_fftw_mpi_local_size_3d fftwf_mpi_local_size_3d
#define p_fftw_alloc_real fftwf_alloc_real
#define p_fftw_alloc_complex fftwf_alloc_complex
#define p_fftw_mpi_plan_dft_r2c_3d fftwf_mpi_plan_dft_r2c_3d
#define p_fftw_mpi_plan_dft_c2r_3d fftwf_mpi_plan_dft_c2r_3d
#define p_fftw_destroy_plan fftwf_destroy_plan
#define p_fftw_execute fftwf_execute
#endif

#if (doubleprecision == 0)
#define pseudo_double double
#define MPI_PSEUDO_DOUBLE MPI_DOUBLE
#define really_double double

#define p_fftw_complex fftw_complex
#define p_fftw_plan fftw_plan
#define p_fftw_mpi_init fftw_mpi_init
#define p_fftw_mpi_local_size_3d fftw_mpi_local_size_3d
#define p_fftw_alloc_real fftw_alloc_real
#define p_fftw_alloc_complex fftw_alloc_complex
#define p_fftw_mpi_plan_dft_r2c_3d fftw_mpi_plan_dft_r2c_3d
#define p_fftw_mpi_plan_dft_c2r_3d fftw_mpi_plan_dft_c2r_3d
#define p_fftw_destroy_plan fftw_destroy_plan
#define p_fftw_execute fftw_execute
#endif

#if (doubleprecision == 1)
#define pseudo_double long double
#define MPI_PSEUDO_DOUBLE MPI_LONG_DOUBLE
#define p_fftw_complex fftwl_complex
#define p_fftw_plan fftwl_plan
#define p_fftw_mpi_init fftwl_mpi_init
#define p_fftw_mpi_local_size_3d fftwl_mpi_local_size_3d
#define p_fftw_alloc_real fftwl_alloc_real
#define p_fftw_alloc_complex fftwl_alloc_complex
#define p_fftw_mpi_plan_dft_r2c_3d fftwl_mpi_plan_dft_r2c_3d
#define p_fftw_mpi_plan_dft_c2r_3d fftwl_mpi_plan_dft_c2r_3d
#define p_fftw_destroy_plan fftwl_destroy_plan
#define p_fftw_execute fftwl_execute

#define really_double double
#endif

//General Parameters
#define hardscale2 0	// 0 without hard scale 2, i.e. \int p_{t/z}^2 p^4 f,	1 with hard scale 2. NOTE: if hardscale2 == 0 then HSEta[2]==0 etc. in output
#define EXPANSION 0	// 0 for static (nonexpanding), 1 for anisotropic expansion with metric $g_{\mu\nu} = diag(1,-1,-1,-\tau^2)$
#define withCorr 1	// 0 without inclusion of Correlation functions in output, 1 with them
#define withBins 1	// 0 without binning the distribution function in momentum space, 1 with binning ;  for static box automatically binning on
#define LAMBDA_MODE 1	// 0 for Lambda == 0 (free Streaming), else Lambda != 0
#define withcharge 0
#define expansion 0
#define finebin 0
#define withEnergyPS 0
#define islanddetect 0
#define EnergybasedID 0 //to tag the islands based on their energy density or not (relevant only if withEnergyPS==1)
#define relicpockets 1

// MPI
int nodes;		// all CPUs used
int ID;			// the ID of the CPU

// Theory parameters -------------------------------------------------------------
#define theory 4

#if (theory==0) //N-component phi4phi6 theory

int Nc=1;					// number of scalar field components

//mass
int Mass_Mode=1;			// 1 - positive, -1 - negative, 0 - 0 mass
pseudo_double M_sqr=1.;	// value of renormalized mass

//quartic interaction
int Lambda_Mode=-1;			// 1 - positive, -1 - negative, 0 - 0 coupling
pseudo_double Lambda=100.;	// value of Lambda

//sextic interaction
pseudo_double Kappa=2.;	// value of sextic coupling = (Lambda*Lambda/Kappa*Kappa)


#endif

#if (theory==1) //N-component phi4 theory

int Nc=1;					// number of scalar field components

//mass
int Mass_Mode=-1;			// 1 - positive, -1 - negative, 0 - 0 mass
pseudo_double M_sqr=1.;		// value of renormalized mass

//quartic interaction
int Lambda_Mode=1;			// 1 - positive, -1 - negative, 0 - 0 coupling
pseudo_double Lambda=5.;	// value of Lambda

//linear term
pseudo_double h = -0.8; //value of the linear coupling = h/sqrt(Lambda)


#endif

#if (theory==2) //N-component sine-gordon theory (cos- potential)

int Nc=1;

pseudo_double M_sqr=0.49;                       //mass    
pseudo_double F_SSB=1e+3*0.7;                   //decay constant

pseudo_double Lambda = (M_sqr / (F_SSB*F_SSB)); //quartic coupling

#endif

#if (theory==3) //N-component massive sine-godron theory (monodromy potential)

int Nc = 1; //better stick to N=1 always

pseudo_double M_sqr = 9.;                  //instanton mass  
pseudo_double F_SSB = sqrt(M_sqr)*1e+24;                   //decay constant
pseudo_double alpha = 0;                        //phase of the cosine
pseudo_double M_mono_sqr = 1.;                  //monodromy mass

pseudo_double Lambda = (M_sqr / (F_SSB*F_SSB)); //quartic coupling

#endif

#if (theory==4) //second minimum

int Nc = 2;				// number of scalar field components

int Mass_Mode = 1;			// 1 - positive, -1 - negative, 0 - 0 mass
pseudo_double M_sqr = 1.;		// value of renormalized mass
					//quartic interaction
int Lambda_Mode = 1;			// 1 - positive, -1 - negative, 0 - 0 coupling
pseudo_double Lambda = 1;		// value of Lambda

					//cubic interaction
pseudo_double g = 2.;

#if (relicpockets == 1)
pseudo_double a_xi = 0.5;
pseudo_double phimax = 2.7;
pseudo_double Deltaphi = 0.2;
pseudo_double M_axion_sqr = 1000.;
#endif


#endif

#if (theory==5) //Monodromy potential

int Nc = 1; //better stick to N=1 always

pseudo_double M_sqr = 1.;                  //instanton mass  
pseudo_double F_SSB = 1.;                   //decay constant
pseudo_double p_mono = -0.5;                //monodromy mass

pseudo_double Lambda = (M_sqr / (F_SSB*F_SSB)); //quartic coupling

#endif
//Tad-pole renormalization (just for static box) -> cancels the quadratic divergence from the tadpole in the equation of motion and in the stress-energy tensor components (energy and pressure)
int renorm_tadpole = 0;			// 0 do not renormalize, 1 renormalize tad pole
int max_renSteps = 2000;		// maximal number of renormalization steps
pseudo_double tad_precision = 1e-6;		// renormalization precision

#if (expansion==1)
pseudo_double H_0=0.0171156;
pseudo_double H_I=F_SSB*1e-5; //Hubble scale of inlfation
#endif

// Specifying Lattice --------------------------------------------------------------

//Spatial Lattice parameters
int N_t = 256;
pseudo_double a_t=0.1;		// spatial lattice spacing in the transverse direction

pseudo_double a_eta;
int N_eta; 		

//Time parameters
pseudo_double T0=0;				// initial time
pseudo_double TMax=1000;			// the maximal time of the simulation after that the simulation ends
pseudo_double TSO=1;			// the time difference between output for energy, pressure, etc. (T.out) and the background field (Phi.out), i.e. of all observables computed in direct space. 
pseudo_double TMO=10;			// the time difference between the output of the distribution function, the correlations and possibly other observables generated in Fourier space. 
pseudo_double dtBaseTrans=0.1;		// usual time base (see lattice.cpp), dtau = dtBaseTrans*a_t
pseudo_double dtBaseLong=0.1;			// time base for the first part of dynamics (see lattice.cpp) when an adaptive time step is used (i.e. when adjustDT == 1), then one has dtau = dtBaseLong*tau*a_eta (see lattice.cpp)
//NOTE: The time variable used in the code is called tau. Its time step dtau is computed by specified ratios dtBase... that we will call time bases, see below.


// NOTE: if adjustDT == 0 (or static box) then dtBaseLong = dtBaseTrans is automatically set, and hence, dtau = dtBaseTrans*a_t is always valid.
#define adjustDT 0	// if 1 then time step is adjusted, otherwise dtau = dtBaseTrans*a_t

// SPECIFYING INITIAL CONDITIONS-----------------------------------------------------------

#define initialize 3    // 0 for ParametricResonance, 1 for SetBox, 2 for large Lambda IV (rather universal), 10 for WhiteNoise (somehow)

#if (initialize==0)

//  DESCRIPTION:  large background field, just vacuum fluctuations initially:
//
//	phi(t=0) = phi0/sqrt(Lambda)
//	pi(t=0) = pi0/sqrt(Lambda)
//	
//	f(t=0,p) = FAmp       NOTE: measured in orders of 1
//

pseudo_double FAmp=0.5;	// the vacuum 1/2
pseudo_double Temperature=1.;	// a possible parameter but leave Temperature==1
pseudo_double phi0=9.6; 	// the initial background field rescaled by sqrt(Lambda)
pseudo_double pi0=0.;		// the derivative of initial background field rescaled by sqrt(Lambda)

#endif

#if (initialize==10)  
// testing case: white noise (somehow implemented, see iv.cpp)
pseudo_double Noise=1.;
pseudo_double phi0=10.;
pseudo_double pi0=0.;

#endif

#if (initialize==1)

//  DESCRIPTION:  any homogeneous background field, large fluctuations Amplitude/Lambda until QS, vacuum fluctuations until 2 QS with amplitude Noise:
//
//  NOTE: please fo not use these initial conditions but instead the more involved initialize = 2 !!!
//
//	phi(t=0) = phi0/sqrt(Lambda)
//	pi(t=0) = pi0/sqrt(Lambda)
//	
//	f(t=0,p) = (sqrt(pchisqr)<QS ? Amplitude/Lambda : ( sqrt(pchisqr)<2QS ? Noise : 0 ) )
//
//	we have pchisqr = pt^2 + chi^2 pz^2 i.e. chi is the initial anisotropy parameter

pseudo_double QS=1.;   	// the initial momentum scale, NOTE: has to be specified!!!
pseudo_double Amplitude=25.;	// NOTE: Amplitude is measured in orders of 1/Lambda if Lambda!=0
pseudo_double chi=1.;		// initial anisotropy parameter
pseudo_double Noise=1.;	// NOTE: measured NOT in orders of 1/lambda but in O(1)!
pseudo_double phi0=9.6; 	// the initial background field rescaled by sqrt(Lambda)
pseudo_double pi0=0.;		// the derivative of initial background field rescaled by sqrt(Lambda)

#endif

#if (initialize==2)

//  DESCRIPTION:  any homogeneous background field, large initial fluctuations until QS and adjustable vacuum. 
//
//	phi(0) = phi0/sqrt(Lambda)
//	pi(0) = pi0/sqrt(Lambda)
//	f(0,p) = (Amplitude/Lambda)*Theta(QS - sqrt(pchisqr)) + Noise \Theta(scale2 - sqrt(pchisqr))
//
//	we have pchisqr = pt^2 + chi^2 pz^2, i.e. chi is the initial anisotropy parameter

int justnoise=0;	// 0 usual initialization for f(t=0,p) as above, 1 for just the noise until scale2 (Amplitude==0) -> simply leave justnoise=0

pseudo_double QS=1;   	// initial momentum scale - NOTE: has to be specified!!!
pseudo_double scale2=0.;	// second scale for vacuum fluctuations (proper renormalization)
// NOTE: if scale2 == 0, then scale2 is the momentum cutoff and the vacuum is occupied up to the cutoff

pseudo_double Amplitude=120.;	// amplitude parameter
pseudo_double chi=1.;		// anisotropy parameter
pseudo_double Noise=0.5;	// vacuum amplitude,  NOTE: measured NOT in orders of 1/lambda but in O(1)!
pseudo_double phi0=0.; 	// the initial background field rescaled by sqrt(Lambda)
pseudo_double pi0=0.;		// the derivative of initial background field rescaled by sqrt(Lambda)

#endif


#if (initialize==3)

//  DESCRIPTION:  Designed for N=2 Case. Initial conditions with large O(2) charge
//
//	phi(0) = phi0/sqrt(Lambda) 
//	pi(0) = pi0/sqrt(Lambda)
//	f(0,p) = (Amplitude/Lambda)*Theta(QS - sqrt(pchisqr)) + Noise \Theta(scale2 - sqrt(pchisqr))
//
//	we have pchisqr = pt^2 + chi^2 pz^2, i.e. chi is the initial anisotropy parameter

int justnoise=0;	// 0 usual initialization for f(t=0,p) as above, 1 for just the noise until scale2 (Amplitude==0) -> simply leave justnoise=0

pseudo_double QS=3.;   	// initial momentum scale - NOTE: has to be specified!!!
pseudo_double scale2=0.;	// second scale for vacuum fluctuations (proper renormalization)
// NOTE: if scale2 == 0, then scale2 is the momentum cutoff and the vacuum is occupied up to the cutoff

pseudo_double Amplitude=0.;	// amplitude parameter
#if (relicpockets==1)
pseudo_double Amplitude2=1.;	// amplitude parameter
#endif
pseudo_double chi=1.;		// anisotropy parameter
pseudo_double Noise=0.0;	// vacuum amplitude,  NOTE: measured NOT in orders of 1/lambda but in O(1)!
pseudo_double phi0_0=0, phi0_1=0.; 	// the initial background field rescaled by sqrt(Lambda), 0.88 = 3.6*sqrt(0.06)
pseudo_double pi0_0=0., pi0_1=0.;	// the derivative of initial background field rescaled by sqrt(Lambda)
pseudo_double M2_long=0., M2_trans=0.; 	//effective mass at initial time, to be calculated

#endif

//-----------------------------------------------------------------------------------------------------------------------------


void GetParameters(int argc, char ** argv)
{
  Konfig arguments(argc,argv);
  
  arguments.Getval("runs", runs);

#if (theory==0)  
   
  arguments.Getval("Nc", Nc); 
  arguments.Getval("M", M_sqr);
  arguments.Getval ("Mmode",Mass_Mode);
  arguments.Getval("L", Lambda);
  arguments.Getval("Lmode", Lambda_Mode);
  arguments.Getval("Kappa", Kappa);
#endif
  
#if (theory==1)
  arguments.Getval("Nc", Nc); 
  arguments.Getval("M", M_sqr);
  arguments.Getval ("Mmode",Mass_Mode);
  arguments.Getval("L", Lambda);
  arguments.Getval("Lmode", Lambda_Mode);
  arguments.Getval("h", h);
#endif
  
#if (theory==2)
  arguments.Getval("Nc", Nc); 
  arguments.Getval("M", M_sqr);
  arguments.Getval("L", Lambda);
  arguments.Getval("f", F_SSB);
  
#endif

#if (theory==3)
  arguments.Getval("Nc", Nc);
  arguments.Getval("M", M_sqr);
  arguments.Getval("L", Lambda);
  arguments.Getval("alpha", alpha);
  arguments.Getval("M_mono", M_mono_sqr);
  arguments.Getval("f", F_SSB);
#endif

#if (theory==4)
  arguments.Getval("Nc", Nc);
  arguments.Getval("M", M_sqr);
  arguments.Getval("Mmode", Mass_Mode);
  arguments.Getval("L", Lambda);
  arguments.Getval("Lmode", Lambda_Mode);
#endif

  arguments.Getval("tadpole",renorm_tadpole);
  arguments.Getval("tadprec",tad_precision);
  arguments.Getval("renSteps",max_renSteps);

  
  arguments.Getval("Nt", N_t);
  arguments.Getval("at", a_t);
  arguments.Getval("T0", T0);
  arguments.Getval("TMax",TMax);
  arguments.Getval("dt",dtBaseLong);
  arguments.Getval("TSO",TSO);
  arguments.Getval("TMO",TMO);

#if (initialize==2)
  arguments.Getval("jn",justnoise);
  arguments.Getval("Q2",scale2);
  arguments.Getval("Ampl",Amplitude);
  arguments.Getval("QS",QS);
  arguments.Getval("chi",chi);
  arguments.Getval("Noise",Noise);
  arguments.Getval("phi0",phi0);
  arguments.Getval("pi0",pi0);
#endif 
 
#if (initialize==3)
  arguments.Getval("jn",justnoise);
  arguments.Getval("Q2",scale2);
  arguments.Getval("Ampl",Amplitude);
  arguments.Getval("QS",QS);
  arguments.Getval("chi",chi);
  arguments.Getval("Noise",Noise);
  arguments.Getval("phi0_0",phi0_0);
  arguments.Getval("phi0_1",phi0_1);
  arguments.Getval("pi0_0",pi0_0);
  arguments.Getval("pi0_1",pi0_1);
#endif 
  
}


void CompleteParameters()
{
  a_eta = a_t;
  N_eta = N_t;
  dtBaseTrans=dtBaseLong;
  
#if (theory==2 || theory==3 || theory==5)
  Lambda = (M_sqr / (F_SSB*F_SSB));
#endif
    
#if (initialize==1 || initialize==2 || initialize==3)
  Noise*=Lambda;
#endif
}


int LASTRUN=0;



#endif //definition of the file
