#ifndef _SEEDBUBBLES_CPP_
#define _SEEDBUBBLES_CPP_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "util.cpp"
#include "comm.cpp"
#include "lattice.cpp"
#include "clusterfft.cpp"

// ===============================
// Bounce profile
// ===============================
typedef struct {
    int N;
    double dr;
    double *r;
    double *phi;
    double* chi_factor;
} BounceProfile;

namespace randGen
{
	void init()
	{
		srand48(time(NULL)*(ID+1));
	}
	
	pseudo_double rng()
	{
		return drand48();
	}
	
	int on=0;
	pseudo_double gauss1,gauss2;
	
	pseudo_double Gauss()
	{
		if(on==1){on=0; return gauss2;}
		pseudo_double x,y,r,c,s;
		int accept=0;
		while(accept==0){
			x=2.0*rng()-1.0; y=2.0*rng()-1.0;
			r=sqrt(sqr(x)+sqr(y));
			if(r<1.0){accept=1;}
		}
		c=x/r; s=y/r;
		gauss1=c*sqrt(-2.0*log(sqr(r)));
		gauss2=(s/c)*gauss1;
		on=1; return gauss1;
	}
}

BounceProfile load_bounce_profile(const char *filename)
{
    BounceProfile bp;

    FILE *f = fopen(filename, "r");
    if (!f) 
    {
        printf("Error opening bounce file\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int count = 0;
    double r, phi, dphi;

    char line[256];
    fgets(line, sizeof(line), f);  // skip header
    // count lines
    while (fscanf(f, "%lf,%lf,%lf", &r, &phi, &dphi) == 3)
        count++;

    rewind(f);

    bp.N = count;

    // ✅ C++ allocation
    bp.r   = new double[count];
    bp.phi = new double[count];
    bp.chi_factor = new double[count];

    double dummy;
    fgets(line, sizeof(line), f); 
    for (int i = 0; i < count; i++)
        fscanf(f, "%lf,%lf,%lf", &bp.r[i], &bp.phi[i], &dummy);

    bp.dr = bp.r[1] - bp.r[0];
    bp.chi_factor[bp.N-1] = 1;
    for (int i=bp.N-2; i>=0; i--)
    {
        dummy = M_axion_sqr * 0.5*(1+tanh((bp.phi[i]-phimax)/Deltaphi));
        bp.chi_factor[i] = bp.chi_factor[i+1] * exp(-bp.dr * sqrt(dummy));
    }

    fclose(f);

    return bp;
}

// free memory
void free_bounce_profile(BounceProfile *bp)
{
    delete[] bp->r;
    delete[] bp->phi;
    delete[] bp->chi_factor;
}

// ===============================
// Interpolation
// ===============================
double interp_phi(double r, BounceProfile *bp)
{
    if (r <= bp->r[0]) return bp->phi[0];
    if (r >= bp->r[bp->N-1]) return 0.0;

    int i = (int)(r / bp->dr);

    // safety correction (in case of slight non-uniformity)
    while (i >= 0 && bp->r[i] > r)
        i--;

    // safety correction (in case of slight non-uniformity)
    while (i < bp->N - 1 && bp->r[i+1] < r)
        i++;

    if (i >= bp->N - 1) i = bp->N - 2;

    double t = (r - bp->r[i]) / bp->dr;

    return (1.0 - t) * bp->phi[i] + t * bp->phi[i+1];
}

double interp_chi_factor(double r, BounceProfile *bp)
{
    if (r <= bp->r[0]) return bp->chi_factor[0];
    if (r >= bp->r[bp->N-1]) return 1.;

    int i = (int)(r / bp->dr);

    // safety correction (in case of slight non-uniformity)
    while (i < bp->N - 1 && bp->r[i+1] < r)
        i++;

    if (i >= bp->N - 1) i = bp->N - 2;

    double t = (r - bp->r[i]) / bp->dr;

    return (1.0 - t) * bp->chi_factor[i] + t * bp->chi_factor[i+1];
}

// ===============================
// Points
// ===============================
typedef struct 
{
    double x, y, z;
} Point;

double periodic_distance(Point a, Point b, double L) 
{
    double dx = fabs(a.x - b.x);
    double dy = fabs(a.y - b.y);
    double dz = fabs(a.z - b.z);

    if (dx > L / 2.0) dx = L - dx;
    if (dy > L / 2.0) dy = L - dy;
    if (dz > L / 2.0) dz = L - dz;

    return sqrt(dx*dx + dy*dy + dz*dz);
}

double rand_uniform(double L)
{
    return randGen::rng() * L;
}

int generate_points(Point *points, int N, double L, double d_min) 
{
    int count = 0;
    int max_attempts = 10000000;
    int attempts = 0;

    while (count < N && attempts < max_attempts) 
    {
        attempts++;

        Point candidate = {
            randGen::rng() * L,
            randGen::rng() * L,
            randGen::rng() * L
        };

        int valid = 1;

        for (int i = 0; i < count; i++) 
        {
            if (periodic_distance(candidate, points[i], L) < d_min) 
            {
                valid = 0;
                break;
            }
        }

        if (valid)
            points[count++] = candidate;
    }

    if (count < N) {
        printf("Warning: Only placed %d/%d points\n", count, N);
        return 0;
    }

    return 1;
}

void generate_and_broadcast_points(Point *points, int N,
                                   double L, double d_min,
                                   int rank, MPI_Comm comm)
{
    if (rank == 0) 
    {
        if (!generate_points(points, N, L, d_min)) 
        {
            printf("Error generating points\n");
            MPI_Abort(comm, 1);
        }
    }

    MPI_Bcast(points, N * sizeof(Point), MPI_BYTE, 0, comm);
}

// ===============================
// Bubble initialization
// ===============================
double periodic_delta(double dx, double L)
{
    if (dx >  L/2.0) dx -= L;
    if (dx < -L/2.0) dx += L;
    return dx;
}

void initialize_bubbles(Point *points, int Np, BounceProfile *bp, double wall_r, double pocket_wall_r)
{
    double L = N_t * a_t;
    double phi_vac = bp->phi[bp->N - 1];   // ✅ important!

    for (int z = 1; z < N_eta_loc + 1; z++)
    for (int y = 0; y < N_t; y++)
    for (int x = 0; x < N_t; x++)
    {
        int Kooo = pos(x, y, z);

        double xp = x * a_t;
        double yp = y * a_t;
        double zp = (ID * LocalLattice::N_eta_loc + (z - 1)) * a_t;

        for (int p = 0; p < Np; p++)
        {
            double dx = periodic_delta(xp - points[p].x, L);
            double dy = periodic_delta(yp - points[p].y, L);
            double dz = periodic_delta(zp - points[p].z, L);

            double r = sqrt(dx*dx + dy*dy + dz*dz);
#if (seedbubbles == 1)
            double r_from_pocket_wall = pocket_wall_r - r;

            Phi[Kooo] += interp_phi(wall_r + r_from_pocket_wall, bp) - phi_vac;   // or overwrite if preferred
            Phi[Kooo+1] *= interp_chi_factor(wall_r + r_from_pocket_wall, bp);  // apply chi factor
            Pi[Kooo+1] *= interp_chi_factor(wall_r + r_from_pocket_wall, bp);  // apply chi factor
#else
            Phi[Kooo] += interp_phi(r, bp) - phi_vac;   // or overwrite if preferred
            Phi[Kooo+1] *= interp_chi_factor(r, bp);  // apply chi factor
            Pi[Kooo+1] *= interp_chi_factor(r, bp);  // apply chi factor
#endif
        }
    }
}

#endif