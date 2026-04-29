#include <cmath>

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "util.cpp"
#include "comm.cpp"
#include "lattice.cpp"
#include "clusterfft.cpp"

#include "seedbubbles.cpp"



////// INITIALIZE MODES IN MOMENTUM SPACE ////////
void SetInitialConditions()
{
	// NOTE: internal structure: py until N_t/2+1, px until N_t, pz until LocalLattice::N_eta_loc (0 <= pz < ...N_eta_loc)

	//INITIALIZE ARRAYS IN MOMENTUM AND COORDINATE SPACE
	std::complex<pseudo_double> **AFieldP=new std::complex<pseudo_double>*[Nc];
	std::complex<pseudo_double> **EFieldP=new std::complex<pseudo_double>*[Nc];

	pseudo_double **AFieldX=new pseudo_double*[Nc];
	pseudo_double **EFieldX=new pseudo_double*[Nc];

	for(int a=0;a<Nc;a++)
	{
		AFieldP[a]=new std::complex<pseudo_double>[(N_t/2+1)*N_t*LocalLattice::N_eta_loc];
		EFieldP[a]=new std::complex<pseudo_double>[(N_t/2+1)*N_t*LocalLattice::N_eta_loc];

		AFieldX[a]=new pseudo_double[N_t*N_t*LocalLattice::N_eta_loc];
		EFieldX[a]=new pseudo_double[N_t*N_t*LocalLattice::N_eta_loc];
	}
	
	//INITIALIZE LOCAL VARIABLES
	pseudo_double ppx,ppy,nu,pt,pchisqr,psqr,Omega;

	int pzGlobal, posMom;
	
	pseudo_double Occupation;
	
	pseudo_double gauss1,gauss2;
	pseudo_double phase1,phase2;
	
	std::complex<pseudo_double> F1;
	std::complex<pseudo_double> F2;
	
	const pseudo_double timetau=taufunc();
	
	randGen::init();
		
	//INITIALIZE VALUES WITHOUT RESPECTING REALITY
	for(int pz=0;pz<LocalLattice::N_eta_loc;pz++){

		pzGlobal=pz+ID*LocalLattice::N_eta_loc;

		//setting nu
		if (pzGlobal==0) nu=0;
		else
		{	//nu!=0	
			if (pzGlobal<N_eta/2+1) nu=2.0/(a_eta)*sin((pzGlobal)*acos(-1.)/N_eta);
			else nu=-2.0/(a_eta)*sin((pzGlobal)*acos(-1.)/N_eta);
		}



		for(int px=0;px<N_t;px++)
		{

			//setting ppx
			if (px<N_t/2+1) ppx=2.0/(a_t)*sin(px*acos(-1.)/N_t);
			else ppx=-2.0/(a_t)*sin(px*acos(-1.)/N_t);



			for(int py=0;py<N_t/2+1; py++)
			{

				
				//setting ppy
				ppy=2.0/(a_t)*sin(py*acos(-1.)/N_t);

								
				//pt
				pt = sqrt(sqr(ppx)+sqr(ppy));

				//psqr
				psqr = sqr(pt)+sqr(nu);
			

				//position of fields in momentum space
				posMom = py+px*(N_t/2+1)+pz*N_t*(N_t/2+1);
				

				//SET MODES TO ZERO TO AVOID DOUBLE INITIALIZATION OF SELF-INVERSE MODES
				for(int a=0;a<Nc;a++)
				{
					AFieldP[a][posMom]=0.0;
					EFieldP[a][posMom]=0.0;
				}

#if (withcharge==1)
				//INITIALIZE MODES FOR ALL COLORS
				for(int a=0;a<Nc;a++){
				  
					if (Nc==2)
					{
					  if (a==1)
					  {
					//Complex phase separately for lambda=1,2 modes
					  phase2=phase1; 
					  phase1=2.0*acos(-1.)*randGen::rng();
					
					//Gaussian ampltiude for lambda=1,2 modes
					  gauss2=gauss1;
					  gauss1=fabs(randGen::Gauss());
					  }
					    
					}
					else
					{
					//Complex phase separately for lambda=1,2 modes
					  phase1=2.0*acos(-1.)*randGen::rng(); phase2=2.0*acos(-1.)*randGen::rng();
					
					//Gaussian ampltiude for lambda=1,2 modes
					  gauss1=fabs(randGen::Gauss()); gauss2=fabs(randGen::Gauss());
					}


					////// SETTING THE OCCUPATION //////


					// PARAMETRIC RESONANCE
					#if (initialize==0)

						if(a==0){
							Omega=sqrt(psqr+fabs(M_sqr)-3.*sqr(Temperature)*sqr(phi0)/(6.*Nc));
						}
						else{
							Omega=sqrt(psqr+fabs(M_sqr)-sqr(Temperature)*sqr(phi0)/(6.*Nc));
						}

						Occupation = FAmp*Lambda;

					#endif	// end initialize 0



					// BOX INITIALIZATION
					#if (initialize==1)

						//pchisqr
							// pchisqr (= pt² + chi²*(nu/t)² )
							pchisqr = sqr(pt)+sqr(chi)*sqr(nu);
						
						
						//Amplitude of Modes
						if(sqrt(pchisqr)<QS){Occupation=Amplitude;}
						
						//Chopped Plateau
						else{
							if(sqrt(pchisqr)<2.0*QS){Occupation=Noise;}
							else {Occupation=0.0;}
						}


						//Omega
						if(a==0){
							Omega=sqrt(psqr+fabs(M_sqr)-3.*sqr(phi0)/(6.*Nc));
						}
						else{
							Omega=sqrt(psqr+fabs(M_sqr)-sqr(phi0)/(6.*Nc));
						}


					#endif	// end initialize 1




					// SMALL N LARGE LAMBDA BOX INITIALIZATION
					#if (initialize==2)

						pchisqr = sqr(pt)+sqr(nu);

						if (justnoise==0)
						{
						   
							//Amplitude of Modes (+ Noise)
							if(sqrt(pchisqr)<QS)
								Occupation=Amplitude+Noise;
							
							//Chopped Plateau
							else{
								if (scale2==0)
									Occupation=Noise;
								else
								{
									if(sqrt(pchisqr)<scale2)
										Occupation=Noise;
									else 
										Occupation=0.0;
								}
							}
						}
						else 
						{ // just the noise
							if (scale2==0)
								Occupation=Noise;
							else
							{
								if(sqrt(pchisqr)<scale2)
									Occupation=Noise;
								else 
									Occupation=0.0;
							}
						}


						//Omega
						if(a==0)
							Omega=sqrt(psqr+fabs(M_sqr)-3.*sqr(phi0)/(6.*Nc));
						else
							Omega=sqrt(psqr+fabs(M_sqr)-sqr(phi0)/(6.*Nc));


					#endif	// end initialize 2

						
					// Possibility for nonzero charge
					#if (initialize==3)

						pchisqr = sqr(pt)+sqr(nu);

						if (justnoise==0)
						{
						   
							//Amplitude of Modes (+ Noise)
							if(sqrt(pchisqr)<QS)
								Occupation=Amplitude+Noise;
							
							//Chopped Plateau
							else{
								if (scale2==0)
									Occupation=Noise;
								else
								{
									if(sqrt(pchisqr)<scale2)
										Occupation=Noise;
									else 
										Occupation=0.0;
								}
							}
						}
						else 
						{ // just the noise
							if (scale2==0)
								Occupation=Noise;
							else
							{
								if(sqrt(pchisqr)<scale2)
									Occupation=Noise;
								else 
									Occupation=0.0;
							}
						}


						//Omega
						if(a==0)
							Omega=sqrt(psqr+M2_long);
						else
							Omega=sqrt(psqr+M2_trans);


					#endif	// end initialize 3
						


					//TODO
					//Complex Amplitude including occupation number
					F1=(pseudo_double)(sqrt(Occupation)*gauss1)*std::complex<pseudo_double>(cos(phase1),sin(phase1));
					F2=(pseudo_double)(sqrt(Occupation)*gauss2)*std::complex<pseudo_double>(cos(phase2),sin(phase2));
			
					  
				

					//INITIALIZE MODES IN MOMENTUM SPACE


					
						//USE CONTINUUM EXPRESSION FOR A AND E FIELDs
						if(fabs(Omega) < 1e-8){	// approximatelay Omega==0
							AFieldP[a][posMom]=0;
							EFieldP[a][posMom]=0;
							//AFieldP[a][posMom]=F1;
							//EFieldP[a][posMom]=F2;

							std::cerr << "Omega is zero!" << std::endl;
						}
						else{
						

							AFieldP[a][posMom]=F1/(pseudo_double)(sqrt(Omega));
							EFieldP[a][posMom]=F2*(pseudo_double)(sqrt(Omega));						
						}
				}
#endif
#if (withcharge==0)
				//INITIALIZE MODES FOR ALL COLORS
				for(int a=0;a<Nc;a++)
				{
				  
					
					//Complex phase separately for lambda=1,2 modes
					phase1=2.0*acos(-1.)*randGen::rng(); phase2=2.0*acos(-1.)*randGen::rng();
					
					//Gaussian ampltiude for lambda=1,2 modes
					gauss1=fabs(randGen::Gauss()); gauss2=fabs(randGen::Gauss());
					


					////// SETTING THE OCCUPATION //////


					// PARAMETRIC RESONANCE
#if (initialize==0)

						if(a==0){
							Omega=sqrt(psqr+fabs(M_sqr)-3.*sqr(Temperature)*sqr(phi0)/(6.*Nc));
						}
						else{
							Omega=sqrt(psqr+fabs(M_sqr)-sqr(Temperature)*sqr(phi0)/(6.*Nc));
						}

						Occupation = FAmp*Lambda;

#endif	// end initialize 0



					// BOX INITIALIZATION
#if (initialize==1)

						//pchisqr
							// pchisqr (= pt² + chi²*(nu/t)² )
							pchisqr = sqr(pt)+sqr(chi)*sqr(nu);
						
						
						//Amplitude of Modes
						if(sqrt(pchisqr)<QS){Occupation=Amplitude;}
						
						//Chopped Plateau
						else{
							if(sqrt(pchisqr)<2.0*QS){Occupation=Noise;}
							else {Occupation=0.0;}
						}


						//Omega
						if(a==0){
							Omega=sqrt(psqr+fabs(M_sqr)-3.*sqr(phi0)/(6.*Nc));
						}
						else{
							Omega=sqrt(psqr+fabs(M_sqr)-sqr(phi0)/(6.*Nc));
						}


#endif	// end initialize 1




					// SMALL N LARGE LAMBDA BOX INITIALIZATION
#if (initialize==2)

						pchisqr = sqr(pt)+sqr(nu);

						if (justnoise==0)
						{
						   
							//Amplitude of Modes (+ Noise)
							if(sqrt(pchisqr)<QS)
								Occupation=Amplitude+Noise;
							
							//Chopped Plateau
							else{
								if (scale2==0)
									Occupation=Noise;
								else
								{
									if(sqrt(pchisqr)<scale2)
										Occupation=Noise;
									else 
										Occupation=0.0;
								}
							}
						}
						else 
						{ // just the noise
							if (scale2==0)
								Occupation=Noise;
							else
							{
								if(sqrt(pchisqr)<scale2)
									Occupation=Noise;
								else 
									Occupation=0.0;
							}
						}


						//Omega
						if(a==0)
							Omega=sqrt(psqr+fabs(M_sqr)-3.*sqr(phi0)/(6.*Nc));
						else
							Omega=sqrt(psqr+fabs(M_sqr)-sqr(phi0)/(6.*Nc));


#endif	// end initialize 2

						
#if (initialize==3)

					pchisqr = sqr(pt)+sqr(nu);
					if (justnoise==0)
					{
//Amplitude of Modes (+ Noise)
#if (expansion==1)
						if (psqr>0. && psqr<(QS*QS))
						{
							pseudo_double omegap=sqrt(psqr+(1.+2.*H_0*T0)*M2_long);
							Occupation = 2.*acos(-1.)*acos(-1.)*2.1e-9*phi0_0*phi0_0/psqr/sqrt(psqr)*omegap;
						}
#endif
#if (expansion==0)
						if(pchisqr<QS*QS)
						{
		   					Occupation=Amplitude+Noise;
#if (relicpockets==1)
							if (a==1)
								Occupation= Amplitude2 + Noise; 
#endif
						}
#endif
							//Chopped Plateau
						else
						{
							if (scale2==0)
								Occupation=Noise;
							else
							{
								if(sqrt(pchisqr)<scale2)
									Occupation=Noise;
								else 
									Occupation=0.0;
							}
						}
					}
					else 
					{ // just the noise
						if (scale2==0)
							Occupation=Noise;
						else
						{
							if(sqrt(pchisqr)<scale2)
								Occupation=Noise;
						else
							Occupation=0.0;
						}
					}

#if (expansion==0)
					//Omega
					if(a==0)
					  Omega=sqrt(psqr+M2_long);
					else
					  Omega=sqrt(psqr+M2_trans);
#endif
#if (expansion==1)
					//Omega
					if(a==0)
					  Omega=sqrt(psqr+(1.+2.*H_0*T0)*M2_long);
					else
					  Omega=sqrt(psqr+(1.+2.*H_0*T0)*M2_trans);
#endif

#endif	// end initialize 3
                
#if (initialize==4)

					pchisqr = sqr(pt)+sqr(nu);
					if (justnoise==0)
					{
//Amplitude of Modes (+ Noise)

						if(pchisqr<QS*QS)
		   					Occupation=Temperature/sqrt(psqr+M2_long) + Noise;
							//Chopped Plateau
						else
						{
							if (scale2==0)
								Occupation=Noise;
							else
							{
								if(sqrt(pchisqr)<scale2)
									Occupation=Noise;
								else 
									Occupation=0.0;
							}
						}
					}
					else 
					{ // just the noise
						if (scale2==0)
							Occupation=Noise;
						else
						{
							if(sqrt(pchisqr)<scale2)
								Occupation=Noise;
							else 
								Occupation=0.0;
						}
					}
					//Omega
					if(a==0)
					  Omega=sqrt(psqr+M2_long);
					else
					  Omega=sqrt(psqr+M2_trans);

#endif	// end initialize 4
		
    
						


					//TODO
					//Complex Amplitude including occupation number
					F1=(pseudo_double)(sqrt(Occupation)*gauss1)*std::complex<pseudo_double>(cos(phase1),sin(phase1));
					F2=(pseudo_double)(sqrt(Occupation)*gauss2)*std::complex<pseudo_double>(cos(phase2),sin(phase2));
			
					//INITIALIZE MODES IN MOMENTUM SPACE
					
						//USE CONTINUUM EXPRESSION FOR A AND E FIELDs
					if(fabs(Omega) < 1e-8)
					{	// approximatelay Omega==0
						AFieldP[a][posMom]=0;
						EFieldP[a][posMom]=0;
						//AFieldP[a][posMom]=F1;
						//EFieldP[a][posMom]=F2;

						std::cerr << "Omega is zero!" << std::endl;
					}
					else
					{
						AFieldP[a][posMom]=F1/(pseudo_double)(sqrt(Omega));
#if (expansion==1)
						EFieldP[a][posMom]=AFieldP[a][posMom]*H_0;
#endif
#if (expansion==0)
						EFieldP[a][posMom]=F2*(pseudo_double)(sqrt(Omega));
#endif	
					}
				}

#endif

			}
		}
	}
	
	//SYNCHRONIZE
	MPI_Barrier(MPI_COMM_WORLD);
	
	//NOTE: Do not forget that 0 <= pz < Lo....::N_eta_loc!!!	


	//MAKE SURE THE RESULT IS REAL FOR PHI-FIELDS
	for(int a=0;a<Nc;a++)
	{
			FFT::ensureRealResult(AFieldP[a]);
			MPI_Barrier(MPI_COMM_WORLD);
	}
	
	//MAKE SURE THE RESULT IS REAL FOR acos(-1.)-FIELDS
	for(int a=0;a<Nc;a++)
	{
			FFT::ensureRealResult(EFieldP[a]);
			MPI_Barrier(MPI_COMM_WORLD);
	}
	
	//SYNCHRONIZE
	MPI_Barrier(MPI_COMM_WORLD);
	



	//FOURIER TRANSFORM PHI-FIELDS
	for(int a=0;a<Nc;a++)
	{
			FFT::fftPtoX(AFieldP[a],AFieldX[a]);
			MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//FOURIER TRANSFORM acos(-1.)-FIELDS
	for(int a=0;a<Nc;a++)
	{
			FFT::fftPtoX(EFieldP[a],EFieldX[a]);
			MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
	//SET VALUES ON THE LATTICE
	pseudo_double NormFactor=sqrt(sqr(N_t*a_t)*N_eta*a_eta);
	
	int posLOC, posBUFF;

	//NOTE: HERE Z FROM 0 BUT POSITION IN Z+1
	for(int z=0;z<LocalLattice::N_eta_loc;z++)
	{
		for(int y=0;y<N_t;y++)
		{
			for(int x=0;x<N_t;x++)
			{
				

					//Field position and position in the buffers	
					posLOC=pos(x,y,z+1);
					posBUFF=x+y*N_t+z*N_t*N_t;
#if (initialize==3 || initialize==4)
					LocalLattice::Pi[posLOC]=EFieldX[0][posBUFF]*NormFactor+pi0_0;
#if (expansion==1)
					LocalLattice::Pi[posLOC]=EFieldX[0][posBUFF]*NormFactor+pi0_0 + phi0_0*H_0;
#endif
					LocalLattice::Phi[posLOC]=AFieldX[0][posBUFF]*NormFactor+phi0_0;
					if (Nc>1)
					{
					  LocalLattice::Pi[posLOC+1]=EFieldX[1][posBUFF]*NormFactor+pi0_1;
					  LocalLattice::Phi[posLOC+1]=AFieldX[1][posBUFF]*NormFactor+phi0_1;
					}
					for(int a=2;a<Nc;a++)
					{

						//SET E-FIELDS
						LocalLattice::Pi[a+posLOC]=EFieldX[a][posBUFF]*NormFactor;
					
						//SET A-FIELDS
						LocalLattice::Phi[a+posLOC]=AFieldX[a][posBUFF]*NormFactor;
						
						
					}
#else
					//SET E-FIELDS
					LocalLattice::Pi[posLOC]=EFieldX[0][posBUFF]*NormFactor+pi0;
				
					//SET A-FIELDS
					LocalLattice::Phi[posLOC]=AFieldX[0][posBUFF]*NormFactor+phi0;
					
					for(int a=1;a<Nc;a++)
					{

						//SET E-FIELDS
						LocalLattice::Pi[a+posLOC]=EFieldX[a][posBUFF]*NormFactor;
					
						//SET A-FIELDS
						LocalLattice::Phi[a+posLOC]=AFieldX[a][posBUFF]*NormFactor;
						
						
					}	
#endif
				
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);



    if (ID==0) std::cout << "Add some bubbles!\n";


	// >>> LOAD BOUNCE PROFILE HERE <<<
	BounceProfile bp = load_bounce_profile("bounce.txt");

	double phi_mid = 0.5 * (bp.phi[0] + bp.phi[bp.N - 1]);
	double wall_r = 0.0;
	for (int j = 0; j < bp.N; j++)
	   	if (bp.phi[j] < phi_mid)
    	{
        	wall_r = bp.r[j];
        	break;
    	}
	
	int N_bubbles = 1;                 // your choice
	double L = N_t * a_t;
	double d_min = 3.*wall_r;
#if (seedpockets==1)
	d_min = 3.* pocket_wall_r;
#endif

	Point *points = new Point[N_bubbles];

	generate_and_broadcast_points(points, N_bubbles, L, d_min, ID, MPI_COMM_WORLD);

	if (ID==0)
		cout << "Initializing bubbles...\n";

#if (seedpockets==1)
	initialize_bubbles(points, N_bubbles, &bp, wall_r, pocket_wall_r);
#else
	initialize_bubbles(points, N_bubbles, &bp, 0, 0);
#endif
	MPI_Barrier(MPI_COMM_WORLD);
	delete[] points;
	free_bounce_profile(&bp);
	MPI_Barrier(MPI_COMM_WORLD);
		
	//Exchange boundaries
	Communication::exchange();
	MPI_Barrier(MPI_COMM_WORLD);
	
	for(int a=0;a<Nc;a++)
	{
		delete[] AFieldP[a];
		delete[] EFieldP[a];
		
		delete[] AFieldX[a];
		delete[] EFieldX[a];
	}
	
	//if (ID==0) std::cerr << "also deleted!\n";

	delete[] AFieldP;
	delete[] EFieldP;
	
	delete[] AFieldX;
	delete[] EFieldX;


} //end of setIV()

 
