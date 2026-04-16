#include "parameters.cpp"

void fix_tadpole()
{
  if (renorm_tadpole == 0)
  {
//SET BARE MASSES
    LocalLattice::BareM_sqr_long = M_sqr;
    LocalLattice::BareM_sqr_trans = M_sqr;
#if (relicpockets==1)
    LocalLattice::BareM_sqr_trans = 0.;
#endif
#if (theory==3)
    LocalLattice::BareM_mono_sqr_long = M_mono_sqr;
    LocalLattice::BareM_mono_sqr_trans = M_mono_sqr;
#endif
    
#if (initialize==3 || initialize==4)

    M2_long = M_sqr;
    M2_trans = M_sqr;
#if (relicpockets==1)
    M2_trans = 0;
#endif
#if (theory==3)
    M2_long+=M_mono_sqr;
    M2_trans+=M_mono_sqr;  
#endif

#endif
    MPI_Barrier(MPI_COMM_WORLD);
    return;
    
  }
  
   MPI_Barrier(MPI_COMM_WORLD); 
}

 
