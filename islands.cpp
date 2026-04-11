#ifndef _ISLANDCPP_
#define _ISLANDCPP_
#include <fstream>
#include <cmath>
#include <queue>
#include <vector>
#include <algorithm>
using namespace std;

#include "util.cpp"
#include "fieldanalysis.cpp"

struct island
{
  int index;	//this is a label attached to a given island, to make it identifiable

  bool disp;	//if the island got split by the nodes or boundaries, disp=1 corresponds to the one where the "center" of the bubble is.
  vector<point> center;  //location of the center
  vector<pseudo_double> V, at_time;

  vector<pseudo_double> PotEn, KinEn, GradEn;
  vector<pseudo_double> avPhi;
  vector<pseudo_double> Peak_value;

  vector<int> addr, addl;  

  int status;
};

//compare them according to their volume?
bool operator < (island& B1, island& B2)
{
  return (B1.V.back()>B2.V.back());
}

namespace island_analysis
{
  pseudo_double thres = 1.125; //the threshold value for island identification
  long long number_of_islands; //a counter for assigning
 
  vector< island > found_islands;
	
  int *b; //store the indices of the islands on the grid		
  queue<int> q;
    
  void init()
  {
    b=new int[Nc*N_t*N_t*(LocalLattice::N_eta_loc+2)];
	  number_of_islands = 0;
  }
  
  void fin()
  {
	  found_islands.clear();
    delete[] b;
  }

  std::ofstream Safe_islands, Safe_islands2, Expanding_bubbles, Init_islands;
}

bool is_bubble(int i)
{
  return ( island_analysis::found_islands[i].GradEn.back() + island_analysis::found_islands[i].PotEn.back() < 0.) ; //Could also be 0, but to be safe.
//	return ((island_analysis::found_islands[i].V.back() > 4.*acos(-1.)*10.*10.*10. / 3.));//This one is the old one, but no so good.
}

bool is_oscillon(int i)
{
  return ( island_analysis::found_islands[i].at_time.back() - island_analysis::found_islands[i].at_time.front()> 40.);
}

bool is_island(int K)
{
	// return (LocalLattice::Phi[K]>island_analysis::thres);  //based on the field value (only positive), best at the initial time
#if (withEnergyPS==0)	
  return (LocalLattice::Phi[K] * LocalLattice::Phi[K] * 0.5 + LocalLattice::Pi[K] * LocalLattice::Pi[K] * 0.5 > island_analysis::thres*island_analysis::thres*0.5); //including the kinetic term
#endif
#if (withEnergyPS==1)
#if (EnergybasedID==0)
  return (LocalLattice::Phi[K] * LocalLattice::Phi[K] * 0.5 + LocalLattice::Pi[K] * LocalLattice::Pi[K] * 0.5 > island_analysis::thres*island_analysis::thres*0.5); //including the kinetic term
#endif
#if (EnergybasedID==1)
  return (LocalLattice::Delta[K]>island_analysis::thres);//oscillons
#endif
#endif	
  // return ((LocalLattice::Phi[K]-LocalLattice::avgphi[0]) * (LocalLattice::Phi[K]-LocalLattice::avgphi[0]) * 0.5 + (LocalLattice::Pi[K]-LocalLattice::avgpi[0]) * (LocalLattice::Pi[K]-LocalLattice::avgpi[0]) * 0.5> island_analysis::thres*island_analysis::thres*0.5); //strange one
}
  
// this one does a full search for the extent of the island
void BFS(int P, int fflag, int flag, island *B) //fflag is the index that will be assigned in the end, flag is the array index in found_bubbles
{
  using namespace LocalLattice;
  using namespace island_analysis;
  
  int NN=0;
  pseudo_double cPotEn=0.,cGradEn=0.,cKinEn=0.,cavPhi=0.;
  
  int cnt=P;
  int x=0, y=0, z=0;
  
  b[P]=flag;
  q.push(P);
  
  while (!q.empty())
  {
    int Kooo=q.front();
    q.pop();

#if (EnergybasedID==0 || withEnergyPS==0)  
	  if (Phi[Kooo] * Phi[Kooo] >Phi[cnt] * Phi[cnt]) //find the peak
      cnt=Kooo;
#endif
#if (EnergybasedID==1 && withEnergyPS==1)  
	  if (Delta[Kooo]>Delta[cnt]) //find the peak
      cnt=Kooo;
#endif

    z=Kooo/(Nc*N_t*N_t);
    y=(Kooo - Nc*N_t*N_t*z)/(Nc*N_t);
    x=(Kooo-Nc*N_t*N_t*z-Nc*N_t*y)/Nc;
         
    int Kuoo=pos(mod(x+1,N_t), y, z);
    int Kdoo=pos(mod(x-1,N_t), y, z);
    int Kouo=pos(x, mod(y+1,N_t), z);
    int Kodo=pos(x, mod(y-1,N_t), z);
    int Koou=pos(x, y, z+1);
    int Kood=pos(x, y, z-1);
        
    NN++;

    cGradEn+= 0.5*( sqr(Phi[Kuoo]-Phi[Kooo]) + sqr(Phi[Kouo]-Phi[Kooo]) + sqr(Phi[Koou]-Phi[Kooo]) ) / sqr(a_t);

#if (expansion==0)
    cKinEn+= 0.5*sqr(Pi[Kooo]);
#endif 
#if (expansion==1)
    cKinEn+= ( 0.5*sqr(Pi[Kooo]) + 0.5*sqr(H_0 / scale_factor)*sqr(Phi[Kooo]) - Phi[Kooo] * Pi[Kooo] * H_0 / scale_factor);
#endif

#if (theory==0)
#if (expansion==0)
    cPotEn+=(pseudo_double)(Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo])+(pseudo_double)(Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc))+(sqr(sqr(Phi[Kooo]))*sqr(Phi[Kooo])/(720.0*Nc*Nc*Kappa*Kappa)) ;
#endif
#if (expansion==1)
    cPotEn+=(Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
    cPotEn+=(Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc))+(sqr(sqr(Phi[Kooo]))*sqr(Phi[Kooo])/(720.0*Nc*Nc*Kappa*Kappa*scale_factor*scale_factor)) ;
#endif
#endif

#if (theory==1)
#if (expansion==0)
    cPotEn+=(  h*Phi[Kooo] + (pseudo_double)(Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo])+(pseudo_double)(Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc)) );
#endif
#if (expansion==1)
    cPotEn+=(Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
    cPotEn+=(  h*Phi[Kooo]*scale_factor*scale_factor*scale_factor + (Lambda_Mode)*(sqr(sqr(Phi[Kooo]))/(24.0*Nc)) ) ;
#endif
#endif

#if (theory==2)
#if (expansion==0)
    cPotEn+=sqr(BareM_sqr_long)*(1.-cos( sqrt(sqr(Phi[Kooo])/BareM_sqr_long) ));
#endif
#if (expansion==1)
    cPotEn+=sqr(scale_factor*scale_factor*BareM_sqr_long)*(1.-cos(  Phi[Kooo] /(scale_factor*sqrt(BareM_sqr_long)) ));
#endif
#endif

#if (theory==3)
#if (expansion==0)
    cPotEn+=0.5*BareM_mono_sqr_long*sqr(Phi[Kooo])+sqr(BareM_sqr_long)*(1.-cos( sqrt(sqr(Phi[Kooo])/BareM_sqr_long) ));
#endif
#if (expansion==1)
    cPotEn += (0.5*BareM_mono_sqr_long*( sqr(scale_factor)  )*sqr(Phi[Kooo])   );
    cPotEn += sqr(scale_factor*scale_factor*BareM_sqr_long)*(1.-cos(  Phi[Kooo] /(scale_factor*sqrt(BareM_sqr_long)) ));
#endif
#endif

#if (theory==4)
#if (expansion==0)
	  cPotEn += (pseudo_double)(Mass_Mode)*0.5*BareM_sqr_long*sqr(Phi[Kooo]) - g*sqr(Phi[Kooo])*Phi[Kooo] / 6. + (pseudo_double)(Lambda_Mode)*(sqr(sqr(Phi[Kooo])) / (24.0*Nc));
#endif
#if (expansion==1)
	  cPotEn += (Mass_Mode)*0.5*BareM_sqr_long*scale_factor*scale_factor*sqr(Phi[Kooo]);
    cPotEn += (((-g*scale_factor*sqr(Phi[Kooo])*Phi[Kooo]) / (6.0)) + (Lambda_Mode)*(sqr(sqr(Phi[Kooo])) / (24.0*Nc)));
#endif
#endif

#if (theory==5)
#if (expansion==0)
    cPotEn += sqr(BareM_sqr_long)*0.5 / p_mono * (pow((1. + Phi[Kooo] * Phi[Kooo] / BareM_sqr_long), p_mono) - 1.);
#endif
#if (expansion==1)
    cPotEn += sqr(scale_factor*scale_factor*BareM_sqr_long)*0.5/p_mono * (pow((1. + Phi[Kooo] * Phi[Kooo]/BareM_sqr_long/scale_factor/scale_factor), p_mono) - 1.); 
#endif
#endif

    cavPhi+= Phi[Kooo];
    
    if (is_island(Kuoo) && b[Kuoo]==-1)
    {
      b[Kuoo]=flag;
      q.push(Kuoo);
    }
    
    if (is_island(Kdoo) && b[Kdoo]==-1)
    {
      b[Kdoo]=flag;
      q.push(Kdoo);
    }
    
    if (is_island(Kouo) && b[Kouo]==-1)
    {
      b[Kouo]=flag;
      q.push(Kouo);
    }
    
    if (is_island(Kodo) && b[Kodo]==-1)
    {
      b[Kodo]=flag;
      q.push(Kodo);
    }
    
    if (z<N_eta_loc)
      if (is_island(Koou) && b[Koou]==-1)
      {
	      b[Koou]=flag;
	      q.push(Koou);
      }
   
    if (z>1)
      if (is_island(Kood) && b[Kood]==-1)
      {
	      b[Kood]=flag;
	      q.push(Kood);
      }
    
  }
  
  cavPhi/=(pseudo_double)(NN);
  pseudo_double cV=NN*a_t*a_t*a_eta;

  cPotEn *= (a_t*a_t*a_eta);
  cKinEn *= (a_t*a_t*a_eta);
  cGradEn *= (a_t*a_t*a_eta);
  
  point Bp;
  Bp.pos=cnt; Bp.id=ID;
  B->center.push_back(Bp);
  
  B->PotEn.push_back(cPotEn);
  B->KinEn.push_back(cKinEn);
  B->GradEn.push_back(cGradEn);
  B->avPhi.push_back(cavPhi);
  B->V.push_back(cV);
#if (withEnergyPS==0 || EnergybasedID==0)
  B->Peak_value.push_back(Phi[cnt]);
#endif
#if (withEnergyPS==1 && EnergybasedID==1)
  B->Peak_value.push_back(Delta[cnt]);
#endif
  B->addr.push_back(-1); B->addl.push_back(-1);
  B->at_time.push_back(tau);
  B->index = fflag;
  B->status = 0;
  if (B->center.size()>5)
  { 
     B->PotEn.erase(B->PotEn.begin()+2);
     B->center.erase(B->center.begin()+2);
     B->KinEn.erase(B->KinEn.begin()+2);
     B->GradEn.erase(B->GradEn.begin()+2);
     B->avPhi.erase(B->avPhi.begin()+2);
     B->V.erase(B->V.begin()+2);
     B->Peak_value.erase(B->Peak_value.begin()+2);
     B->addr.erase(B->addr.begin()+2);
     B->addl.erase(B->addl.begin()+2);
     B->at_time.erase(B->at_time.begin()+2);
  }
}

//This function is a bit long, but it does communicate islands along the borders
void Analysis_Neighbours()
{
  using namespace LocalLattice;
  using namespace island_analysis;
  
//Communicate the borders
  MPI_Status bstatus;
  MPI_Request bRecvRL;
  MPI_Request bSendRL;
  MPI_Request bRecvLR;
  MPI_Request bSendLR;
  MPI_Irecv(&b[pos(0,0,N_eta_loc+1)], N_t*N_t*Nc, MPI_INT, MPI::right, 31, MPI_COMM_WORLD, &bRecvRL);
  MPI_Isend(&b[pos(0,0,1)], N_t*N_t*Nc, MPI_INT, MPI::left, 31, MPI_COMM_WORLD, &bSendRL);
  MPI_Wait(&bRecvRL,&bstatus);
  MPI_Wait(&bSendRL,&bstatus);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Irecv(&b[pos(0,0,0)], N_t*N_t*Nc, MPI_INT, MPI::left, 52, MPI_COMM_WORLD, &bRecvLR);
  MPI_Isend(&b[pos(0,0,N_eta_loc)], N_t*N_t*Nc, MPI_INT, MPI::right, 52, MPI_COMM_WORLD, &bSendLR);
  MPI_Wait(&bRecvLR,&bstatus);
  MPI_Wait(&bSendLR,&bstatus);	
  MPI_Barrier(MPI_COMM_WORLD);

//Communicate with number of islands
  MPI_Request ntr, nfr, ntl, nfl;
  int nr=0,nl=0, N=found_islands.size();
  MPI_Irecv(&nl,1 , MPI_INT, MPI::left, 114, MPI_COMM_WORLD, &nfl);
  MPI_Irecv(&nr,1 , MPI_INT, MPI::right, 115, MPI_COMM_WORLD, &nfr);
  MPI_Isend(&N, 1, MPI_INT, MPI::right, 114, MPI_COMM_WORLD, &ntr);
  MPI_Isend(&N, 1, MPI_INT, MPI::left, 115, MPI_COMM_WORLD, &ntl);
  MPI_Wait(&nfl,&bstatus);
  MPI_Wait(&nfr,&bstatus);
  MPI_Wait(&ntl,&bstatus);
  MPI_Wait(&ntr,&bstatus);	
  MPI_Barrier(MPI_COMM_WORLD);

//for BFS  
  bool **Dl=new bool* [nl+N];
  bool **Dr=new bool* [nr+N];
  int *adl=new int[nl];
  int *adr=new int[nr];
  for (int i=0;i<nl+N;i++)
  {
    Dl[i]=new bool[nl+N];
    for (int j=0;j<nl+N;j++)
      Dl[i][j]=0;
  }
  for (int i=0;i<nr+N;i++)
  {
    Dr[i]=new bool[nr+N];
    for (int j=0;j<nr+N;j++)
      Dr[i][j]=0;
  }
  for (int i=0;i<nl;i++)
    adl[i]=-1;  
  for (int i=0;i<nr;i++)
    adr[i]=-1;
  
  int addcountr=0,addcountl=0, cl=0,cr=0;
  MPI_Request RRRd,SSSd, SSSu, RRRu;
  
//count the number of items to send and to recieve  
  for (int i=0;i<N_t;i++)
    for (int j=0;j<N_t;j++)
    {
      if (b[pos(i,j,1)]!=-1 && b[pos(i,j,0)]!=-1)
      {
	Dl[ b[pos(i,j,1)] ][ N+b[pos(i,j,0)] ]=1;
	Dl[ N+b[pos(i,j,0)] ][ b[pos(i,j,1)] ]=1;
	
	if (found_islands[ b[pos(i,j,1)] ].addl.back()==-1)
	{
	  found_islands[ b[pos(i,j,1)] ].addl.back()=-2;
	  addcountl++;
	}
	if (adl[ b[pos(i,j,0)] ] ==-1)
	{
	  adl[ b[pos(i,j,0)] ]=-2;
	  cl++;
	}
      }
      if (b[pos(i,j,N_eta_loc)]!=-1 && b[pos(i,j,N_eta_loc+1)]!=-1)
      {
	Dr[ b[pos(i,j,N_eta_loc)] ][ N+b[pos(i,j,N_eta_loc+1)] ]=1;
	Dr[ N+b[pos(i,j,N_eta_loc+1)] ][ b[pos(i,j,N_eta_loc)] ]=1;
	
	if (found_islands[ b[pos(i,j,N_eta_loc)] ].addr.back()==-1)
	{
	//  std::cout<<ID<<' '<<b[pos(i,j,N_eta_loc)]<<' '<<ID+1<<' '<<b[pos(i,j,N_eta_loc+1)]<<std::endl;
	  found_islands[ b[pos(i,j,N_eta_loc)] ].addr.back()=-2;	
	  addcountr++;
	}
	if (adr[ b[pos(i,j,N_eta_loc+1)] ] ==-1)
	{
	  adr[ b[pos(i,j,N_eta_loc+1)] ]=-2;
	  cr++;
	}
      }
    }

// data to send and recieve    
  pseudo_double *tempdr, *tempdl, *tempur, *tempul;
  
 // if (addcountr!=0)
  //{
    tempur=new pseudo_double[6*addcountr];
    tempdr=new pseudo_double[6*cr];
    MPI_Irecv(tempdr, 6*cr, MPI_DOUBLE, MPI::right, 3, MPI_COMM_WORLD, &RRRd);
 // } 

  //if (addcountl!=0)
 // {
    tempdl=new pseudo_double[6*addcountl];
    tempul=new pseudo_double[6*cl];
    MPI_Irecv(tempul, 6*cl, MPI_DOUBLE, MPI::left, 6, MPI_COMM_WORLD, &SSSu);
 // }

      
  addcountr=0; 
  addcountl=0;
  
  for (int i=0;i<N_t;i++)
    for (int j=0;j<N_t;j++)
    {
      if (b[pos(i,j,1)]!=-1 && b[pos(i,j,0)]!=-1)
	if (found_islands[ b[pos(i,j,1)] ].addl.back()==-2)
	{
	  found_islands[ b[pos(i,j,1)] ].addl.back()=-3;
	  
	  tempdl[addcountl*6]=found_islands[ b[pos(i,j,1)] ].V.back();
	  tempdl[addcountl*6+1]=found_islands[ b[pos(i,j,1)] ].avPhi.back();
	  tempdl[addcountl*6+2]=found_islands[ b[pos(i,j,1)] ].GradEn.back();
	  tempdl[addcountl*6+3]=found_islands[ b[pos(i,j,1)] ].PotEn.back();
	  tempdl[addcountl*6+4]=found_islands[ b[pos(i,j,1)] ].KinEn.back();
#if (withEnergyPS==0 || EnergybasedID==0)
 	  tempdl[addcountl*6+5]=Phi[ found_islands[ b[pos(i,j,1)] ].center.back().pos ];
#endif
#if (withEnergyPS==1 && EnergybasedID==1)
 	  tempdl[addcountl*6+5]=Delta[ found_islands[ b[pos(i,j,1)] ].center.back().pos ];//oscillons
#endif

	  addcountl++;
	}
      if (b[pos(i,j,N_eta_loc)]!=-1 && b[pos(i,j,N_eta_loc+1)]!=-1)
	if (found_islands[ b[pos(i,j,N_eta_loc)] ].addr.back()==-2)
	{
	  found_islands[ b[pos(i,j,N_eta_loc)] ].addr.back()=-3;
	    
	  tempur[addcountr*6]=found_islands[ b[pos(i,j,N_eta_loc)] ].V.back();
	  tempur[addcountr*6+1]=found_islands[ b[pos(i,j,N_eta_loc)] ].avPhi.back();
	  tempur[addcountr*6+2]=found_islands[ b[pos(i,j,N_eta_loc)] ].GradEn.back();
	  tempur[addcountr*6+3]=found_islands[ b[pos(i,j,N_eta_loc)] ].PotEn.back();
	  tempur[addcountr*6+4]=found_islands[ b[pos(i,j,N_eta_loc)] ].KinEn.back();
#if (withEnergyPS==0 || EnergybasedID==0)
	  tempur[addcountr*6+5]=Phi[ found_islands[ b[pos(i,j,N_eta_loc)] ].center.back().pos ];
#endif
#if (withEnergyPS==1 && EnergybasedID==1)
	  tempur[addcountr*6+5]=Delta[ found_islands[ b[pos(i,j,N_eta_loc)] ].center.back().pos ];//oscillons
#endif
	  
	  addcountr++;
	}
    } 
         
//  if (addcountr!=0)
    MPI_Isend(tempur, 6*addcountr, MPI_DOUBLE, MPI::right, 6, MPI_COMM_WORLD, &RRRu);
//  if (addcountl!=0)
    MPI_Isend(tempdl, 6*addcountl, MPI_DOUBLE, MPI::left, 3, MPI_COMM_WORLD, &SSSd);

 // if (addcountr!=0)
    MPI_Wait(&RRRd,&bstatus);
 // if (addcountl!=0)
  //{
    MPI_Wait(&SSSu, &bstatus);
    MPI_Wait(&SSSd,&bstatus);
 // }
 // if (addcountr!=0)
    MPI_Wait(&RRRu, &bstatus);
  MPI_Barrier(MPI_COMM_WORLD);
 
  //BFS
  queue<int> QQQ;
  
  int n_of_conn_compl=0;
  int *coll=new int [N+nl];
  for (int i=0;i<N+nl;i++)
    coll[i]=-1;
  for (int i=0;i<N+nl;i++)
    if (coll[i]==-1)
    {
      coll[i]=n_of_conn_compl;
      n_of_conn_compl++;
      QQQ.push(i);
      while (!QQQ.empty())
      {
	int k=QQQ.front();
	QQQ.pop();
	for (int j=0;j<N+nl;j++)
	  if (Dl[k][j] && coll[j]==-1)
	  {
	    QQQ.push(j);
	    coll[j]=coll[i];
	  }
      }
    }
    
  int n_of_conn_compr=0;
  int *colr=new int [N+nr];
  for (int i=0;i<N+nr;i++)
    colr[i]=-1;
  for (int i=0;i<N+nr;i++)
    if (colr[i]==-1)
    {
      colr[i]=n_of_conn_compr;
      n_of_conn_compr++;
      QQQ.push(i);
      while (!QQQ.empty())
      {
	int k=QQQ.front();
	QQQ.pop();
	for (int j=0;j<N+nr;j++)
	  if (Dr[k][j] && colr[j]==-1)
	  {
	    QQQ.push(j);
	    colr[j]=colr[i];
	  }
      }
    }
    
  for (int i=0;i<N+nr;i++)
    delete [] Dr[i];
  for (int i=0;i<N+nl;i++)
    delete [] Dl[i];
  delete[] Dl;
  delete[] Dr;
   
//Calculater Parameters for each connected component    
  pseudo_double *VCCl=new pseudo_double[n_of_conn_compl];
  pseudo_double *avPhiCCl=new pseudo_double[n_of_conn_compl];
  pseudo_double *GradEnCCl=new pseudo_double[n_of_conn_compl];
  pseudo_double *PotEnCCl=new pseudo_double[n_of_conn_compl];
  pseudo_double *KinEnCCl=new pseudo_double[n_of_conn_compl];
  pseudo_double *Peak_valueCCl=new pseudo_double[n_of_conn_compl];
  
  pseudo_double *VCCr=new pseudo_double[n_of_conn_compr];
  pseudo_double *avPhiCCr=new pseudo_double[n_of_conn_compr];
  pseudo_double *GradEnCCr=new pseudo_double[n_of_conn_compr];
  pseudo_double *PotEnCCr=new pseudo_double[n_of_conn_compr];
  pseudo_double *KinEnCCr=new pseudo_double[n_of_conn_compr];
  pseudo_double *Peak_valueCCr=new pseudo_double[n_of_conn_compr];
  
  for (int i=0;i<n_of_conn_compl;i++)
  {
    VCCl[i]=0;avPhiCCl[i]=0;GradEnCCl[i]=0;PotEnCCl[i]=0;KinEnCCl[i]=0;Peak_valueCCl[i]=0;
  }
  for (int i=0;i<n_of_conn_compr;i++)
  {
    VCCr[i]=0;avPhiCCr[i]=0;GradEnCCr[i]=0;PotEnCCr[i]=0;KinEnCCr[i]=0;Peak_valueCCr[i]=0;
  }
  
  cl=0;
  cr=0;
  for (int i=0;i<N_t;i++)
    for (int j=0;j<N_t;j++)
    {
      if (b[pos(i,j,1)]!=-1 && b[pos(i,j,0)]!=-1)
      {
	if (adl[ b[pos(i,j,0)] ]==-2)
	{
	  adl[ b[pos(i,j,0)] ]=-3;
	  VCCl[ coll[N+b[pos(i,j,0)]] ]+=tempul[cl*6];
	  avPhiCCl[ coll[N+b[pos(i,j,0)]] ]+=tempul[cl*6]*tempul[cl*6+1];
	  GradEnCCl[ coll[N+b[pos(i,j,0)]] ]+=tempul[cl*6+2];
	  PotEnCCl[ coll[N+b[pos(i,j,0)]] ]+=tempul[cl*6+3];
	  KinEnCCl[ coll[N+b[pos(i,j,0)]] ]+=tempul[cl*6+4];
	  if (fabs(tempul[cl*6+5])>fabs( Peak_valueCCl[ coll[N+b[pos(i,j,0)] ]] ))
	    Peak_valueCCl[ coll[N+b[pos(i,j,0)] ] ]=tempul[cl*6+5];
	  cl++;
	}
      }
      if (b[pos(i,j,N_eta_loc)]!=-1 && b[pos(i,j,N_eta_loc+1)]!=-1)
      {
	if (adr[ b[pos(i,j,N_eta_loc+1)] ]==-2)
	{
	  adr[ b[pos(i,j,N_eta_loc+1)] ]=-3;
	  VCCr[ colr[N+b[pos(i,j,N_eta_loc+1)]] ]+=tempdr[cr*6];
	  avPhiCCr[ colr[N+b[pos(i,j,N_eta_loc+1)]] ]+=tempdr[cr*6]*tempdr[cr*6+1];
	  GradEnCCr[ colr[N+b[pos(i,j,N_eta_loc+1)]] ]+=tempdr[cr*6+2];
	  PotEnCCr[ colr[N+b[pos(i,j,N_eta_loc+1)]] ]+=tempdr[cr*6+3];
	  KinEnCCr[ colr[N+b[pos(i,j,N_eta_loc+1)]] ]+=tempdr[cr*6+4];
	  if (fabs(tempdr[cr*6+5])>fabs( Peak_valueCCr[ colr[N+b[pos(i,j,N_eta_loc+1)] ]] ))
	    Peak_valueCCr[ colr[N+b[pos(i,j,N_eta_loc+1)] ] ]=tempdr[cr*6+5];
	  cr++;
	}
      }
    }
  
    
  for (int i=0;i<found_islands.size();i++)
  {
    VCCl[ coll[i] ]+=found_islands[i].V.back();
    avPhiCCl[ coll[i] ]+=found_islands[i].avPhi.back()*found_islands[i].V.back();
    GradEnCCl[ coll[i] ]+=found_islands[i].GradEn.back();
    PotEnCCl[ coll[i] ]+=found_islands[i].PotEn.back();
    KinEnCCl[ coll[i] ]+=found_islands[i].KinEn.back();
   
    VCCr[ colr[i] ]+=found_islands[i].V.back();
    avPhiCCr[ colr[i] ]+=found_islands[i].avPhi.back()*found_islands[i].V.back();
    GradEnCCr[ colr[i] ]+=found_islands[i].GradEn.back();
    PotEnCCr[ colr[i] ]+=found_islands[i].PotEn.back();
    KinEnCCr[ colr[i] ]+=found_islands[i].KinEn.back();
    
    if (fabs( found_islands[i].Peak_value.back() )>fabs( Peak_valueCCl[ coll[i] ] ))
      Peak_valueCCl[ coll[i] ]=found_islands[i].Peak_value.back();
    if (fabs( found_islands[i].Peak_value.back() )>fabs( Peak_valueCCr[ colr[i] ] ))
      Peak_valueCCr[ colr[i] ]=found_islands[i].Peak_value.back();
  }
       
// assing info to the islands on the node  
  for (int i=0;i<found_islands.size();i++)
  {
    found_islands[i].avPhi.back() = ( avPhiCCl[coll[i]]+avPhiCCr[colr[i]] - found_islands[i].V.back()*found_islands[i].avPhi.back()) /( VCCl[coll[i]]+VCCr[colr[i]] -found_islands[i].V.back());
    found_islands[i].V.back() = VCCl[coll[i]]+VCCr[colr[i]] - found_islands[i].V.back();
    found_islands[i].GradEn.back() = GradEnCCl[coll[i]]+GradEnCCr[colr[i]] - found_islands[i].GradEn.back();
    found_islands[i].KinEn.back() = KinEnCCl[coll[i]]+KinEnCCr[colr[i]] - found_islands[i].KinEn.back();
    found_islands[i].PotEn.back() = PotEnCCl[coll[i]]+PotEnCCr[colr[i]] - found_islands[i].PotEn.back();
    found_islands[i].disp=0;
    if (( fabs(found_islands[i].Peak_value.back()-Peak_valueCCr[ colr[i] ]) <1e-6) && ( fabs(found_islands[i].Peak_value.back()-Peak_valueCCl[ coll[i] ])<1e-6))
      found_islands[i].disp=1;
    else
      found_islands[i].Peak_value.back() = ( fabs(Peak_valueCCl[coll[i]])>fabs(Peak_valueCCr[colr[i]]) )?(Peak_valueCCl[coll[i]]):(Peak_valueCCr[colr[i]]);
  }
  
  

  delete[] VCCl;
  delete[] avPhiCCl;
  delete[] GradEnCCl;
  delete[] PotEnCCl;
  delete[] KinEnCCl;
  delete[] Peak_valueCCl;
  
  delete[] VCCr;
  delete[] avPhiCCr;
  delete[] GradEnCCr;
  delete[] PotEnCCr;
  delete[] KinEnCCr;
  delete[] Peak_valueCCr;
  
  delete[] coll;
  delete[] colr;
  
  delete [] adl;
  delete[] adr;
    
  
//  if (addcountr!=0)
 // {
    delete[] tempdr;
    delete[] tempur;
//  }
 // if (addcountl!=0)
  //{
    delete[] tempdl;
    delete[] tempul;
 // }
  MPI_Barrier(MPI_COMM_WORLD);
}

//This function checks if there are island-like configurations which have shrunk, and if yes prints them out.
void Old_islands()
{
  using namespace LocalLattice;
  using namespace island_analysis;
  
  char * fname=new char[256];
  char * fname2=new char[256];
  snprintf(fname,256,"./output/%i/non_critical_islands.txt",LASTRUN);
  snprintf(fname2,256, "./output/%i/islands.txt", LASTRUN);
  
  for (int rID=0;rID<nodes;rID++)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ID==rID)
    {
      Safe_islands.open(fname,std::ifstream::app);
      Init_islands.open(fname2,std::ifstream::app);
      for (int i=0;i<found_islands.size();)
      {
        if (!is_island(found_islands[i].center.back().pos)) //if the island is no longer of interest
        {
            if (found_islands[i].disp)
            {
                for (int j=0;j<found_islands[i].V.size();j++)	    
                {
                    Safe_islands<<found_islands[i].at_time[j]<<' '<<found_islands[i].center[j].x()<<' '<<found_islands[i].center[j].y()<<' '<<found_islands[i].center[j].z()<<"        ";
                    Safe_islands<<found_islands[i].Peak_value[j]<<' '<<found_islands[i].avPhi[j]<<' '<<found_islands[i].V[j]<<' ';
                    Safe_islands<<found_islands[i].GradEn[j]<<' '<<found_islands[i].PotEn[j]<<' '<<found_islands[i].KinEn[j]<<' ' << found_islands[i].index << std::endl;
                }	  
                Safe_islands<<std::endl;
            }
            Init_islands<<pow(found_islands[i].V[0]*3./4./acos(-1.), 1./3.)<<' '<<found_islands[i].Peak_value[0]<<' '<<found_islands[i].at_time[0]<<' '<<found_islands[i].index << ' '<<found_islands[i].status<<' '<<pow(found_islands[i].V[0]*3./4./acos(-1.), 1./3.)/sqrt( log( found_islands[i].Peak_value[0]* found_islands[i].Peak_value[0] / thres / thres) ) <<std::endl;
//            int imax = std::max_element(found_islands[i].Peak_value.begin(), found_islands[i].Peak_value.end())-found_islands[i].Peak_value.begin();
 //           if (found_islands[i].Peak_value[imax]<0)
  //              imax = std::min_element(found_islands[i].Peak_value.begin(), found_islands[i].Peak_value.end())-found_islands[i].Peak_value.begin();
   //         Safe_islands2<<found_islands[i].at_time[imax]<<' '<<found_islands[i].Peak_value[imax]<<' '<<found_islands[i].V[imax]<<std::endl;
      
            found_islands[i].V.clear();
            found_islands[i].avPhi.clear();
            found_islands[i].GradEn.clear();
            found_islands[i].KinEn.clear();
            found_islands[i].PotEn.clear();
            found_islands[i].Peak_value.clear();
            found_islands[i].at_time.clear();
            found_islands[i].addr.clear();
            found_islands[i].addl.clear();
            found_islands[i].center.clear();
            found_islands.erase(found_islands.begin()+i);
        }
        else
            i++;
      }
      Safe_islands.close();
      Init_islands.close();
    }    
  }
  
  delete[] fname;
  delete[] fname2;
}

void print_bubble_history(int i)
{
	using namespace LocalLattice;
	using namespace island_analysis;

	char * fname = new char[256];
  char * fname2 = new char[256];
  snprintf(fname, 256, "./output/%i/expanding_bubbles.txt", LASTRUN);
  snprintf(fname2, 256, "./output/%i/islands.txt", LASTRUN);
	Expanding_bubbles.open(fname, std::ifstream::app);
  Init_islands.open(fname2, std::ifstream::app);

  if (found_islands[i].disp)
    Init_islands<<pow(found_islands[i].V[0]*3./4./acos(-1.), 1./3.)<<' '<<found_islands[i].Peak_value[0]<<' '<<found_islands[i].at_time[0]<<' '<<found_islands[i].index << ' '<<1<<' '<<pow(found_islands[i].V[0]*3./4./acos(-1.), 1./3.)/sqrt( log( found_islands[i].Peak_value[0]* found_islands[i].Peak_value[0] / thres / thres) ) <<std::endl;
   
	for (int j = 1;j<found_islands[i].V.size();j++)
	{
		Expanding_bubbles << found_islands[i].at_time[j] << ' ' << found_islands[i].center[j].x() << ' ' << found_islands[i].center[j].y() << ' ' << found_islands[i].center[j].z() << "        ";
		Expanding_bubbles << found_islands[i].Peak_value[j] << ' ' << found_islands[i].avPhi[j] << ' ' << found_islands[i].V[j] << ' ';
		Expanding_bubbles << found_islands[i].GradEn[j] << ' ' << found_islands[i].PotEn[j] << ' ' << found_islands[i].KinEn[j] << ' ' << found_islands[i].index << std::endl;
	}
	Expanding_bubbles << std::endl;

	Expanding_bubbles.close();
    Init_islands.close();
    delete[] fname;
    delete[] fname2;
    
	return;
}



void Oscillon_specific_observables()//must be generalized
{
	using namespace LocalLattice;
	using namespace island_analysis;

  pseudo_double r_osc=0, n_osc=0, avgM_osc=0.;
  for (int i=0;i<found_islands.size();i++)
  {
    if (found_islands[i].disp && is_oscillon(i))
    {
      r_osc+=( found_islands[i].PotEn.back() + found_islands[i].KinEn.back() + found_islands[i].GradEn.back() ); 
      n_osc = n_osc + 1.;
    }
  }
  avgM_osc = r_osc / n_osc;
  r_osc = r_osc/N_t/a_t/N_t/a_t/N_t/a_t/avgEnergy;
  n_osc = n_osc/N_t/a_t/N_t/a_t/N_t/a_t;
    
  ofstream osc_out;
  char * fname = new char[256];
  snprintf(fname, 256, "./output/%i/oscillons.txt", LASTRUN);
  osc_out.open(fname, std::ifstream::app);
  
  osc_out<<LocalLattice::tau<< ' '<<r_osc<<' '<<n_osc<<' '<<avgM_osc<<endl;

	osc_out.close();
  delete[] fname;
    
	return;
}



void Full_Analysis()
{
  using namespace LocalLattice;
  using namespace island_analysis;

#if (withEnergyPS==1 && EnergybasedID==1)
  CalculateVolumeAverages();
  thres = LocalLattice::avgEnergy*4.;
#endif
   
  for (int k = 0; k < Nc*N_t*N_t*(LocalLattice::N_eta_loc+2); k++)
      b[k]=-1;
  
  sort(found_islands.begin(), found_islands.end());
  
  Old_islands(); 
   
  for (int i=0;i<found_islands.size();i++)
  {
    if (is_island(found_islands[i].center.back().pos) && b[ found_islands[i].center.back().pos ]!=-1)   // this island has been merged to another one from the cycle, change the condition
    {
      found_islands[i].V.clear();
      found_islands[i].avPhi.clear();
      found_islands[i].GradEn.clear();
      found_islands[i].KinEn.clear();
      found_islands[i].PotEn.clear();
      found_islands[i].Peak_value.clear();
      found_islands[i].at_time.clear();
      found_islands[i].addr.clear();
      found_islands[i].addl.clear();
      found_islands[i].center.clear();
      found_islands.erase(found_islands.begin()+i);
      i--;
      continue;
    }

    if (is_island(found_islands[i].center.back().pos) && b[ found_islands[i].center.back().pos ]==-1)  //this bubble is still alive, update its information
		  BFS(found_islands[i].center.back().pos, found_islands[i].index, i, &found_islands[i]);
	
	  if (is_bubble(i) && found_islands[i].status == 0)//Expanding bubble
    {
      print_bubble_history(i);
      found_islands[i].status = 1;
    }

    if (is_oscillon(i) && found_islands[i].status == 0)
      found_islands[i].status = 2;
  }
     
  for (int K = pos(0, 0, 1); K < Nc*sqr(N_t)*(N_eta_loc + 1); K++)
	  if (is_island(K) && b[K] == -1)  // New bubble is born
	  {
		  island bbb;
		  BFS(K, number_of_islands + ID * 1000000, found_islands.size(), &bbb);
		  bbb.disp = 1;
		  number_of_islands++;
		  found_islands.push_back(bbb);
	  }
  
  if (nodes>1)
    Analysis_Neighbours();   
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif
