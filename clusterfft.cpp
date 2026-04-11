#ifndef _FFT_
#define _FFT_

#include <mpi.h>
#include <fftw3-mpi.h>
#include <complex>

#include "parameters.cpp"

int sgn(pseudo_double x)
{
  if (x>0)
    return 1;
  else
    return -1;
}

namespace FFT 
{
//NOTE: in x-space we have the structure y,x,z whereas in p-space,
//	the coordinate order is py,px,pz 
// 0 <= py < N_t/2+1
// 0 <= px < N_t
// 0 <= pz < LocalLattice::N_eta_loc
	
  ptrdiff_t xLow,locNz;
  ptrdiff_t locSize;
	
  pseudo_double *XData;
  fftw_complex *PData;
	
  fftw_plan XtoP,PtoX;

  
  
  void init()
  {
    fftw_mpi_init();
//Split over z Direction, real in x Direction
    locSize=fftw_mpi_local_size_3d(N_eta,N_t,N_t/2+1,MPI_COMM_WORLD,&locNz,&xLow);

    XData=fftw_alloc_real(2*locSize);
    PData=fftw_alloc_complex(locSize);
		
    XtoP=fftw_mpi_plan_dft_r2c_3d(N_eta,N_t,N_t,XData,PData,MPI_COMM_WORLD,FFTW_MEASURE);
    PtoX=fftw_mpi_plan_dft_c2r_3d(N_eta,N_t,N_t,PData,XData,MPI_COMM_WORLD,FFTW_MEASURE);
  }

  void fin()
  {
    fftw_destroy_plan(XtoP);
    fftw_destroy_plan(PtoX);
  }	

  void ensureRealResult(std::complex<pseudo_double> * FieldP)
  {
    MPI_Status status;

    MPI_Request *TransferIV=new MPI_Request[LocalLattice::N_eta_loc*N_t*(N_t/2+1)];
    for(int i=0;i<LocalLattice::N_eta_loc*N_t*(N_t/2+1);i++)
      TransferIV[i]=MPI_REQUEST_NULL;
	
    int globalPz,globalMPz;
    int localMPz,partnerID;
    int globalPx,globalMPx,localMPx;
    int globalPy,globalMPy,localMPy;

    //TRANSFER VALUES TO OTHER NODES AND SET VALUES ON OWN NODE
    for(int localPz=0;localPz<LocalLattice::N_eta_loc;localPz++)
      for(int localPx=0;localPx<N_t;localPx++)
	for(int localPy=0;localPy<N_t/2+1;localPy++)
	{
	  globalPz=localPz+ID*LocalLattice::N_eta_loc;
	  globalMPz=mod(N_eta-globalPz,N_eta);
					
	  partnerID=globalMPz/LocalLattice::N_eta_loc;
	  localMPz=globalMPz-partnerID*LocalLattice::N_eta_loc;
					
	  globalPx=localPx;
	  globalMPx=mod(N_t-globalPx,N_t);
	  localMPx=globalMPx;
					
	  globalPy=localPy;
	  globalMPy=mod(N_t-globalPy,N_t);
	  localMPy=globalMPy;

//Modes for which reality must be enforce by hand
	  if(globalMPy==globalPy)
	  {
//If the p=-p mode is on the same lattice handle locally
	    if(partnerID==ID)
	    {
//If mode has p=-p it has to be purely real => Set to absolute value
	      if(globalPx==globalMPx && globalPz==globalMPz)
	      	FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]=abs(FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)])*sgn(real(FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]));
//Set the -p mode to the conugate value of the p mode
	      else 
		FieldP[localMPy+localMPx*(N_t/2+1)+localMPz*N_t*(N_t/2+1)]=conj(FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]);
	    }
					
//If p=-p mode is somewhere else communicate
	    else
	    {
            
#if (longpseudo_double == 1)            
          pseudo_double REAL1[2]={real(FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]),imag(FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)])};
//Send from left to right
	      if(ID<partnerID)
	      {
            MPI_Isend(&REAL1, 2, MPI_DOUBLE, partnerID, localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1), MPI_COMM_WORLD, &TransferIV[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]);
	      }
//Receive from left to right (STILL HAS TO BE CONJUGATED)
	      if(ID>partnerID)
	      {
            MPI_Irecv(&REAL1, 2, MPI_DOUBLE, partnerID, localMPy+localMPx*(N_t/2+1)+localMPz*N_t*(N_t/2+1), MPI_COMM_WORLD, &TransferIV[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]);
            FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]= (REAL1[0], REAL1[1]);
	      }
#endif

#if (longpseudo_double == 0)
//Send from left to right
	      if(ID<partnerID)
	      	MPI_Isend(&FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)], 1, MPI_DOUBLE_COMPLEX, partnerID, localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1), MPI_COMM_WORLD, &TransferIV[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]);
//Receive from left to right (STILL HAS TO BE CONJUGATED)
	      if(ID>partnerID)
            MPI_Irecv(&FieldP[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)], 1, MPI_DOUBLE_COMPLEX, partnerID, localMPy+localMPx*(N_t/2+1)+localMPz*N_t*(N_t/2+1), MPI_COMM_WORLD, &TransferIV[localPy+localPx*(N_t/2+1)+localPz*N_t*(N_t/2+1)]);          
#endif
	    }
	  }
	  
	}

//SYNCHRONIZE , COMPLETE TRANSFER AND COMPLEX CONJUGATE
    MPI_Barrier(MPI_COMM_WORLD);
		
    for(int i=0;i<LocalLattice::N_eta_loc*N_t*(N_t/2+1);i++)
      if(TransferIV[i]!=MPI_REQUEST_NULL)
      {
	MPI_Wait(&TransferIV[i],&status);
	if(ID>(nodes-1)/2)
	  FieldP[i]=conj(FieldP[i]);
      }
	
//SYNCHRONIZE AND CLEANUP
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] TransferIV;
  }

//FFTW from momentum to spatial: In has structure py,px,pz, Out has structure x y z
  void fftPtoX(std::complex<pseudo_double> *In,pseudo_double *Out)
  {
    pseudo_double Norm=(sqr(a_t*N_t)*N_eta*a_eta);
    int posMom;
		
//Set Input
    ptrdiff_t px,py,pz;
    for(pz=0;pz<locNz;pz++)
      for(int px=0;px<N_t;px++)
	for(int py=0;py<(N_t/2+1);py++)
	{
	  posMom=py+px*(N_t/2+1)+pz*N_t*(N_t/2+1);
	  PData[posMom][0]=real(In[posMom]);
	  PData[posMom][1]=imag(In[posMom]);
	}
    for (int k=0; k < locNz*N_t*2*(N_t/2+1); k++)
      XData[k]=0;
    MPI_Barrier(MPI_COMM_WORLD);
	
    
    
//Transformation
    fftw_execute(PtoX);
    MPI_Barrier(MPI_COMM_WORLD);

    
//Set Output
    ptrdiff_t x,y,z;
    for (z=0;z<locNz;z++)
      for (x=0;x<N_t;x++)
	for (y=0;y<N_t;y++)
	  Out[x+y*N_t+z*N_t*N_t]= XData[y+x*(2*(N_t/2+1))+z*N_t*(2*(N_t/2+1))] / Norm;
    MPI_Barrier(MPI_COMM_WORLD);	
    
  }
	
// FFTW from spatial to momentum space: In has structure of Field (see lattice, there the function pos(int) ), Out has structure py,px,pz
  void fftXtoP_direct(int a, pseudo_double * const & In,std::complex<pseudo_double> * &Out)
  {
    pseudo_double Factor=(sqr(a_t)*a_eta);
    
    int posMom,posLOC;
		
//Set Input
    ptrdiff_t x,y,z;
    for (z=0;z<locNz;z++)
      for (x=0;x<N_t;x++)
	for (y=0;y<N_t;y++)
	{
	  posLOC=pos(x,y,z+1);
	  XData[y+x*(2*(N_t/2+1))+z*N_t*(2*(N_t/2+1))]=In[a+posLOC];
	}
    for (int k=0; k < locNz*N_t*(N_t/2+1); k++)
    {
      PData[k][0]=0;
      PData[k][1]=0;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    
    
//Transformation
    fftw_execute(XtoP);
    MPI_Barrier(MPI_COMM_WORLD);
    
    
//Set Output
    ptrdiff_t px,py,pz;
    for(pz=0;pz<locNz;pz++)
      for(int px=0;px<N_t;px++)
	for(int py=0;py<(N_t/2+1);py++)
	{
	  posMom=py+px*(N_t/2+1)+pz*N_t*(N_t/2+1);
	  Out[posMom]=Factor*std::complex<pseudo_double>(PData[posMom][0],PData[posMom][1]);
	}
    MPI_Barrier(MPI_COMM_WORLD);		
  }

//FFTW from momentum to spatial: In has structure py,px,pz, Out has structure x y z
  void fftPtoX_direct(int a, std::complex<pseudo_double> *In, pseudo_double *Out)
  {
	  pseudo_double Norm = (sqr(a_t*N_t)*N_eta*a_eta);
	  int posMom, posLOC;

	  //Set Input
	  ptrdiff_t px, py, pz;
	  for (pz = 0;pz<locNz;pz++)
		  for (int px = 0;px<N_t;px++)
			  for (int py = 0;py<(N_t / 2 + 1);py++)
			  {
				  posMom = py + px*(N_t / 2 + 1) + pz*N_t*(N_t / 2 + 1);
				  PData[posMom][0] = real(In[posMom]);
				  PData[posMom][1] = imag(In[posMom]);
			  }
	  for (int k = 0; k < locNz*N_t * 2 * (N_t / 2 + 1); k++)
		  XData[k] = 0;
	  MPI_Barrier(MPI_COMM_WORLD);

	  //Transformation
	  fftw_execute(PtoX);
	  MPI_Barrier(MPI_COMM_WORLD);

	  //Set Output
	  ptrdiff_t x, y, z;
	  for (z = 0;z < locNz;z++)
		  for (x = 0;x < N_t;x++)
			  for (y = 0;y < N_t;y++)
			  {
				  posLOC = pos(x, y, z + 1);
				  Out[a + posLOC] = XData[y + x*(2 * (N_t / 2 + 1)) + z*N_t*(2 * (N_t / 2 + 1))] / Norm;
			  }
	  MPI_Barrier(MPI_COMM_WORLD);
  }
	
//FFTW from spatial to momentum: In has structure x y z, Out has structure py,px,pz
	void fftXtoP(pseudo_double *In,std::complex<pseudo_double> *Out)
	{
		
		pseudo_double Factor=(sqr(a_t)*a_eta);

		int posMom;
		
		//Set Input
		ptrdiff_t x,y,z;
		for (z=0;z<locNz;z++){
           		for (x=0;x<N_t;x++){
				for (y=0;y<N_t;y++){

					XData[y+x*(2*(N_t/2+1))+z*N_t*(2*(N_t/2+1))]=In[x+y*N_t+z*N_t*N_t];
				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);


		//TODO
		for (int k=0; k < locNz*N_t*(N_t/2+1); k++){
			PData[k][0]=0;
			PData[k][1]=0;
		}
		

		//Perform Transformation
		fftw_execute(XtoP);
		
		MPI_Barrier(MPI_COMM_WORLD);

		
		//Set Output
		ptrdiff_t px,py,pz;
		for(pz=0;pz<locNz;pz++){
			for(int px=0;px<N_t;px++){
				for(int py=0;py<(N_t/2+1);py++){

					posMom=py+px*(N_t/2+1)+pz*N_t*(N_t/2+1);

					Out[posMom]=Factor*std::complex<pseudo_double>(PData[posMom][0],PData[posMom][1]);

				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);		
	}


}
#endif
