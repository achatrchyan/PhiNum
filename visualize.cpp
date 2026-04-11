#ifndef _VISUALIZECPP_
#define _VISUALIZECPP_

#include <iostream>
#include <fstream> 
#include <iomanip>
#include <mpi.h>

#include "parameters.cpp"
#include "lattice.cpp"
#include "dynamics.cpp"



int jJj=0;



//Printing 3d regions
void Print3D(int x,int y,int z,int l, int delta)
{
  if (l==0)
    return;
  
  std::ofstream fout;
  char *fname=new char[256];
  snprintf(fname,256,"./output/%i/3d-%i-%i-%i-%i-%i.txt",LASTRUN,int(LocalLattice::tau),x,y,z,l);
      
  int rID=z/LocalLattice::N_eta_loc;
    
  if (ID==rID)
    fout.open(fname,std::fstream::app);
    
  for (int k=z,kloc=z%LocalLattice::N_eta_loc;k<z+l;k+=delta,kloc+=delta)
  {
    if ( (k%N_t)/LocalLattice::N_eta_loc!=rID)
    {
      if (ID==rID)
        fout.close();
    
      MPI_Barrier(MPI_COMM_WORLD);
      
      rID=(k%N_t)/LocalLattice::N_eta_loc;
      kloc=0;
      
      if (ID==rID)
        fout.open(fname,std::fstream::app);
    }
    
    if (ID==rID)
    {
      for (int i=x;i<x+l;i+=delta)
        for (int j=y;j<y+l;j+=delta)
        {
            pseudo_double FieldHere=0., MomentumHere=0., FieldAmplitudeHere=0., MomentumAmplitudeHere=0., InteractionHere=0., EnergyHere=0., GradientHere=0., ChargeHere=0., DeltaHere=0.;
	  
            int Position= pos( (i%N_t), (j%N_t) , kloc+1 ) ;
            int Kuoo=pos(mod(i+1,N_t), j%N_t, kloc+1);
            int Kouo=pos(i%N_t, mod(j+1,N_t), kloc+1);
            int Koou=pos(i%N_t, j%N_t, kloc+2);
	  
#if (withcharge==1)
            if (Nc==2)
                ChargeHere+=LocalLattice::Phi[Position]*LocalLattice::Pi[Position+1]-LocalLattice::Phi[Position+1]*LocalLattice::Pi[Position];
#endif
	  
            FieldHere=LocalLattice::Phi[Position];	//The first component
            MomentumHere=LocalLattice::Pi[Position];
	  
            for (int a=0;a<Nc;a++)
            {
                FieldAmplitudeHere+=(LocalLattice::Phi[Position]*LocalLattice::Phi[Position]);
                MomentumAmplitudeHere+=(LocalLattice::Pi[Position]*LocalLattice::Pi[Position]);
	    
                GradientHere+= ( sqr( LocalLattice::Phi[Kuoo]-LocalLattice::Phi[Position] )/sqr(a_t) + sqr( LocalLattice::Phi[Kouo]-LocalLattice::Phi[Position] )/sqr(a_t) + sqr( LocalLattice::Phi[Koou]-LocalLattice::Phi[Position] )/sqr(a_t) );

#if (theory==0)
                if (a==0)
                    InteractionHere+= (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_long*sqr(LocalLattice::Phi[Position]) ;
                else
                    InteractionHere+=(pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_trans*sqr(LocalLattice::Phi[Position]);	            
#endif
#if (theory==1)
                if (a==0)
                    InteractionHere+=(  h*LocalLattice::Phi[Position] + (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_long*sqr(LocalLattice::Phi[Position]) );
                else
                    InteractionHere+=(pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_trans*sqr(LocalLattice::Phi[Position]);	  
#endif        
#if (theory==3)
                if (a==0)
                    InteractionHere+=0.5*LocalLattice::BareM_mono_sqr_long*sqr(LocalLattice::Phi[Position]);
                else
                    InteractionHere+=0.5*LocalLattice::BareM_mono_sqr_trans*sqr(LocalLattice::Phi[Position]);	
#endif	  
#if (theory==4)
				if (a == 0)
					InteractionHere += (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_long*sqr(LocalLattice::Phi[Position]);
				else
					InteractionHere += (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_trans*sqr(LocalLattice::Phi[Position]);
#endif
	    
                Position++;
                Kuoo++; Kouo++; Koou++;
            }
	  	  
#if (theory==0)
            InteractionHere+=(pseudo_double)(Lambda_Mode)*(sqr(FieldAmplitudeHere)/(24.0*Nc))+(sqr(FieldAmplitudeHere)*FieldAmplitudeHere/(720.0*Nc*Nc*Kappa*Kappa));
#endif
#if (theory==1)
            InteractionHere+=(pseudo_double)(Lambda_Mode)*(sqr(FieldAmplitudeHere)/(24.0*Nc));
#endif
#if (theory==2)
            InteractionHere+= sqr(LocalLattice::BareM_sqr_long)*(1.-cos( sqrt(FieldAmplitudeHere/LocalLattice::BareM_sqr_long) ));
#endif
#if (theory==3)
            InteractionHere+= sqr(LocalLattice::BareM_sqr_long)*(1.-cos( sqrt(FieldAmplitudeHere/LocalLattice::BareM_sqr_long) ));
#endif
#if (theory==4)
	    InteractionHere += ( ((-g*FieldHere*FieldAmplitudeHere)/(6.0)) +  (pseudo_double) (Lambda_Mode)*(sqr(FieldAmplitudeHere) / (24.0*Nc)) );
#endif
	 
            EnergyHere=0.5*MomentumAmplitudeHere + 0.5*GradientHere + InteractionHere;
            FieldAmplitudeHere=sqrt(FieldAmplitudeHere);
            MomentumAmplitudeHere=sqrt(MomentumAmplitudeHere);
	  
            fout<<i<<' '<<j<<' '<<k<<' '<<FieldHere<<' '<<MomentumHere<<' '<<EnergyHere<<' '<<GradientHere<<' '<<InteractionHere<<' '<<FieldAmplitudeHere<<' '<<MomentumAmplitudeHere;
       
        
#if (withcharge==1)
            if (Nc==2)
                fout<<' '<<ChargeHere;
#endif
            fout<<std::endl;    
	
	  }
    }    
  }
  
  delete[] fname;
  if (ID==rID)
    fout.close();
  MPI_Barrier(MPI_COMM_WORLD);
}

void Print3D_Full(int delta)
{
  Print3D(0,0,0,N_t,delta);
}

void Print3D_Center(int px, int py, int pz,int r, int delta)
{
  int x=px-r;
  int y=py-r;
  int z=pz-r;
  if (x<0)
    x+=N_t;
  if (y<0)
    y+=N_t;
  if (z<0)
    z+=N_t;
  Print3D(x,y,z,2*r, delta);
}

void Print2Dz(point P,int x,int y, int l)
{
  if (ID==P.id)
  {
    if (l==0)
      return;
  
    std::ofstream fout;
    char *fname=new char[256];
    snprintf(fname,256,"./output/%i/2dz-%i.txt",LASTRUN,jJj++);
    fout.open(fname);
    delete[] fname;
    
    for (int i=x;i<x+l;i++)
      for (int j=y;j<y+l;j++)
      {
        pseudo_double FieldHere=0., MomentumHere=0., FieldAmplitudeHere=0., MomentumAmplitudeHere=0., InteractionHere=0., EnergyHere=0., GradientHere=0., ChargeHere=0., DeltaHere=0.;
	  
        int Position= pos( (i%N_t), (j%N_t) , (P.z()%LocalLattice::N_eta_loc)+1 ) ;
        int Kuoo=pos(mod(i+1,N_t), j%N_t, (P.z()%LocalLattice::N_eta_loc)+1);
        int Kouo=pos(i%N_t, mod(j+1,N_t), (P.z()%LocalLattice::N_eta_loc)+1);
        int Koou=pos(i%N_t, j%N_t,(P.z()%LocalLattice::N_eta_loc)+2);

#if (withcharge==1)
        if (Nc==2)
            ChargeHere+=LocalLattice::Phi[Position]*LocalLattice::Pi[Position+1]-LocalLattice::Phi[Position+1]*LocalLattice::Pi[Position];
#endif  
        FieldHere=LocalLattice::Phi[Position];	//The first component
        MomentumHere=LocalLattice::Pi[Position];
	  
#if (withEnergyPS==1)	
        EnergyHere = LocalLattice::Delta[Position];  
#endif
#if (withEnergyPS==0) //remove?
        for (int a=0;a<Nc;a++)
        {
            FieldAmplitudeHere+=(LocalLattice::Phi[Position]*LocalLattice::Phi[Position]);
            MomentumAmplitudeHere+=(LocalLattice::Pi[Position]*LocalLattice::Pi[Position]);
	    
            GradientHere+= ( sqr( LocalLattice::Phi[Kuoo]-LocalLattice::Phi[Position] )/sqr(a_t) + sqr( LocalLattice::Phi[Kouo]-LocalLattice::Phi[Position] )/sqr(a_t) + sqr( LocalLattice::Phi[Koou]-LocalLattice::Phi[Position] )/sqr(a_t) );

#if (theory==0)
            if (a==0)
                InteractionHere+=(pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_long*sqr(LocalLattice::Phi[Position]) ;
            else
                InteractionHere+=(pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_trans*sqr(LocalLattice::Phi[Position]);	            
#endif
#if (theory==1)
            if (a==0)
                InteractionHere+=(  h*LocalLattice::Phi[Position] + (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_long*sqr(LocalLattice::Phi[Position]) );
            else
                InteractionHere+=(pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_trans*sqr(LocalLattice::Phi[Position]);	  
#endif
#if (theory==3)
#if (expansion==0)
            if (a==0)
                InteractionHere+=0.5*LocalLattice::BareM_mono_sqr_long*sqr(LocalLattice::Phi[Position]);
            else
                InteractionHere+=0.5*LocalLattice::BareM_mono_sqr_trans*sqr(LocalLattice::Phi[Position]);
#endif
#if (theory==4)
			if (a == 0)
				InteractionHere += (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_long*sqr(LocalLattice::Phi[Position]);
			else
				InteractionHere += (pseudo_double)(Mass_Mode)*0.5*LocalLattice::BareM_sqr_trans*sqr(LocalLattice::Phi[Position]);
#endif
#if (expansion==1)
            if (a==0)
                InteractionHere+=(0.5*LocalLattice::BareM_mono_sqr_long*( sqr(LocalLattice::scale_factor)+sqr(H_0/LocalLattice::scale_factor)  )*sqr(LocalLattice::Phi[Position])   -   LocalLattice::Phi[Position]*LocalLattice::Pi[Position]*H_0/LocalLattice::scale_factor    );
            else
                InteractionHere+=(0.5*LocalLattice::BareM_mono_sqr_trans*( sqr(LocalLattice::scale_factor)+sqr(H_0/LocalLattice::scale_factor)  )*sqr(LocalLattice::Phi[Position])    -   LocalLattice::Phi[Position]*LocalLattice::Pi[Position]*H_0/LocalLattice::scale_factor     );	
                InteractionHere+=sqr(LocalLattice::scale_factor*LocalLattice::scale_factor*LocalLattice::BareM_sqr_long)*(1.-cos( LocalLattice::Phi[Position]/(LocalLattice::scale_factor*sqrt(LocalLattice::BareM_sqr_long)) ));
#endif
#endif	    

            Position++;
            Kuoo++; Kouo++; Koou++;
        }

#if (theory==0)
        InteractionHere+=(pseudo_double)(Lambda_Mode)*(sqr(FieldAmplitudeHere)/(24.0*Nc))+(sqr(FieldAmplitudeHere)*FieldAmplitudeHere/(720.0*Nc*Nc*Kappa*Kappa));
#endif
#if (theory==1)
        InteractionHere+=(pseudo_double)(Lambda_Mode)*(sqr(FieldAmplitudeHere)/(24.0*Nc));
#endif
#if (theory==2)
        InteractionHere+= sqr(LocalLattice::BareM_sqr_long)*(1.-cos( sqrt(FieldAmplitudeHere/LocalLattice::BareM_sqr_long) ));
#endif
#if (theory==3)
        InteractionHere+= sqr(LocalLattice::BareM_sqr_long)*(1.-cos( sqrt(FieldAmplitudeHere/LocalLattice::BareM_sqr_long) ));
#endif
#if (theory==4)
		InteractionHere += (((-g*FieldHere*FieldAmplitudeHere) / (6.0)) + (pseudo_double)(Lambda_Mode)*(sqr(FieldAmplitudeHere) / (24.0*Nc)));
#endif
	 
        EnergyHere=0.5*MomentumAmplitudeHere + 0.5*GradientHere + InteractionHere;
        FieldAmplitudeHere=sqrt(FieldAmplitudeHere);
        MomentumAmplitudeHere=sqrt(MomentumAmplitudeHere);
#endif
	  
        fout<<i*a_t<<' '<<j*a_t<<' ';
#if (expansion==0)       
        fout<<LocalLattice::tau<<' ';
#endif
#if (expansion==1)
        fout<<LocalLattice::scale_factor<<' ';
#endif
        fout<<FieldHere<<' '<<MomentumHere;
#if (withEnergyPS==1)
        fout<<' '<<EnergyHere;
#endif
        
#if (withcharge==1)
        if (Nc==2)
            fout<<' '<<ChargeHere;
#endif

#if (islanddetect==1)
        fout<<' '<<island_analysis::b[Position];
        if (island_analysis::b[Position]<0 || island_analysis::b[Position]>=island_analysis::found_islands.size())
          fout<<' '<<-1;
        else
          fout<<' '<<island_analysis::found_islands[island_analysis::b[Position]].status;
#endif

        fout<<std::endl;    
	
    }
  fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void Print2Dz_Full(point P)
{
  Print2Dz(P,0,0,N_t);
}

void Print2Dz_Center(point P, int px, int py, int r)
{
  int x=px-r;
  int y=py-r;
  if (x<0)
    x+=N_t;
  if (y<0)
    y+=N_t;
  Print2Dz(P,x,y,2*r);
}

//Printing 2d slices in 3 directions
void Print2Dx3(point P,int x,int y,int z,int l)
{
  if (l==0)
    return;
  
  Print2Dz(P,x,y,l);
  
  std::ofstream fxout,fyout;
  char *fnamex=new char[256], *fnamey=new char[256];
  snprintf(fnamex,256,"./output/%i/2dx-%i-%i-%i-%i-%i.txt",LASTRUN,int(LocalLattice::tau),P.x(),y,z,l);
  snprintf(fnamey,256,"./output/%i/2dy-%i-%i-%i-%i-%i.txt",LASTRUN,int(LocalLattice::tau),x,P.y(),z,l);
      
  int rID=z/LocalLattice::N_eta_loc;
    
  if (ID==rID)
  {
    fxout.open(fnamex,std::fstream::app);
    fyout.open(fnamey, std::fstream::app);
  }
  
  for (int k=z,kloc=z%LocalLattice::N_eta_loc;k<z+l;k++,kloc++)
  {
    if ( (k%N_t)/LocalLattice::N_eta_loc!=rID)
    {
      if (ID==rID)
      {
        fxout.close();
        fyout.close();
      }
    
      MPI_Barrier(MPI_COMM_WORLD);
      
      rID=(k%N_t)/LocalLattice::N_eta_loc;
      kloc=0;
      
      if (ID==rID)
      {
        fxout.open(fnamex,std::fstream::app);
        fyout.open(fnamey,std::fstream::app);
      }
    }
    
    if (ID==rID)
    {
      for (int i=x;i<x+l;i++)
      {
        pseudo_double FieldHere=0., MomentumHere=0., FieldAmplitudeHere=0., MomentumAmplitudeHere=0., InteractionHere=0., EnergyHere=0., GradientHere=0., ChargeHere=0., DeltaHere=0.;
        int Position= pos( (i%N_t), (P.y()%N_t) , kloc+1 ) ;
        
        FieldHere=LocalLattice::Phi[Position];	//The first component
        MomentumHere=LocalLattice::Pi[Position];
	  
        for (int a=0;a<Nc;a++)
        {
            FieldAmplitudeHere+=(LocalLattice::Phi[Position+a]*LocalLattice::Phi[Position+a]);
            MomentumAmplitudeHere+=(LocalLattice::Pi[Position+a]*LocalLattice::Pi[Position+a]);	    
        }
        FieldAmplitudeHere=sqrt(FieldAmplitudeHere);
        MomentumAmplitudeHere=sqrt(MomentumAmplitudeHere);
        
        fyout<<i<<' '<<k<<' '<<FieldHere<<' '<<MomentumHere<<' '<<FieldAmplitudeHere<<' '<<MomentumAmplitudeHere;
#if (withcharge==1) 
        if (Nc==2)
        {
            pseudo_double ChargeHere=LocalLattice::Phi[Position]*LocalLattice::Pi[Position+1]-LocalLattice::Pi[Position]*LocalLattice::Phi[Position+1];
            fyout<<' '<<ChargeHere;
        }  
#endif
        fyout<<std::endl;
      }
      for (int j=y;j<y+l;j++)
      {
        pseudo_double FieldHere=0., MomentumHere=0., FieldAmplitudeHere=0., MomentumAmplitudeHere=0., InteractionHere=0., EnergyHere=0., GradientHere=0., ChargeHere=0., DeltaHere=0.;
        int Position= pos( (P.x()%N_t), (j%N_t) , kloc+1 ) ;
        
        FieldHere=LocalLattice::Phi[Position];	//The first component
        MomentumHere=LocalLattice::Pi[Position];
	  
        for (int a=0;a<Nc;a++)
        {
            FieldAmplitudeHere+=(LocalLattice::Phi[Position+a]*LocalLattice::Phi[Position+a]);
            MomentumAmplitudeHere+=(LocalLattice::Pi[Position+a]*LocalLattice::Pi[Position+a]);	    
        }
        FieldAmplitudeHere=sqrt(FieldAmplitudeHere);
        MomentumAmplitudeHere=sqrt(MomentumAmplitudeHere);
        
        fxout<<j<<' '<<k<<' '<<FieldHere<<' '<<MomentumHere<<' '<<FieldAmplitudeHere<<' '<<MomentumAmplitudeHere;
#if (withcharge==1) 
        if (Nc==2)
        {
            pseudo_double ChargeHere=LocalLattice::Phi[Position]*LocalLattice::Pi[Position+1]-LocalLattice::Pi[Position]*LocalLattice::Phi[Position+1];
            fxout<<' '<<ChargeHere;
        }  
#endif
        fyout<<std::endl;
      }
    }    
  }
  
  delete[] fnamex;
  delete[] fnamey;
  if (ID==rID)
  {
    fxout.close();
    fyout.close();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
}

void Print2Dx3_Full(point P)
{
  Print2Dx3(P,0,0,0,N_t);
}

void Print2Dx3_Center(point P,int r)
{
  int x=P.x()-r;
  int y=P.y()-r;
  int z=P.z()-r;
  if (x<0)
    x+=N_t;
  if (y<0)
    y+=N_t;
  if (z<0)
    z+=N_t;
  Print2Dx3(P,x,y,z,2*r);
}



//-----------------------

void Print_snapshot_of_islands()
{
  using namespace island_analysis;
    
  char * fname = new char[256];
  snprintf(fname, 256, "./output/%i/islands-%i.txt", LASTRUN, jJj++);

  for (int j=0;j<nodes;j++)
  {
    if (ID!=j)
      continue;
    ofstream fout (fname);
    for (int i=0;i<found_islands.size();i++)
    {
      if (found_islands[i].disp)
        fout<<found_islands[i].index << ' '<<found_islands[i].status<<' '<<pow(found_islands[i].V.back()*3./4./acos(-1.), 1./3.)<<' '<<found_islands[i].Peak_value.back()<<' ' << found_islands[i].at_time.back() - found_islands[i].at_time.front()<<' ' << found_islands[i].PotEn.back()+ found_islands[i].KinEn.back()+ found_islands[i].GradEn.back()<<' '<<pow(found_islands[i].V.back()*3./4./acos(-1.), 1./3.)/sqrt( log( found_islands[i].Peak_value[0]* found_islands[i].Peak_value[0] / thres / thres) )<< std::endl;
    }
    fout.close();
  }	
  delete[] fname;
}


#endif
