#include <cmath>

namespace ImportFilter
{
	pseudo_double dp_import;
	pseudo_double W_import[25000];
}

using namespace ImportFilter;

pseudo_double W_gaussian(pseudo_double r, pseudo_double p)
{
	return exp(-p*p*r*r/2.);
}

pseudo_double W_tophat(pseudo_double r, pseudo_double p)
{
	if (p*r == 0)
		return 1;
	return 3.*(sin(p*r) / (p*p*p*r*r*r) - cos(p*r) / (p*p*r*r));
}

pseudo_double W_ktophat(pseudo_double r, pseudo_double p)
{
    pseudo_double P = 2.*acos(-1.)/r;
	if (p<P)
		return 1;
	return 0;
}

void import_filter() //must be improved//
{
	ifstream fin("bounce.txt");
	for (int i = 0;i < 25000;i++)
	{
		double a, b, c, d, e;
		fin >> a >> b >> c >> d >> e;
        W_import[i] = e;
		if (i == 1)
			dp_import = c;
	}
	return;
}

pseudo_double W_import_bounce(pseudo_double p)
{
	if (p > 24999 * dp_import)
		return W_import[24999];
	return W_import[ int(p/dp_import) ];
}


void apply_filter(pseudo_double r_filter, int which_filter) //which_filter: 0 -> tophat, 1 -> gaussian
{
	int posLOC = 0, posARRAY = 0, numCounts = 0, pzGlobal = 0;
	pseudo_double ppx, ppy, nu, pabs;

	for (int a = 0; a < Nc; a++)
	{
		for (int k = 0;k < (N_t / 2 + 1)*N_t*LocalLattice::N_eta_loc; k++)
		{
			MomentumSpace::PhiP[k] = 0.;
			MomentumSpace::PiP[k] = 0.;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		FFT::fftXtoP_direct(a, LocalLattice::Phi, MomentumSpace::PhiP);
		MPI_Barrier(MPI_COMM_WORLD);
		FFT::fftXtoP_direct(a, LocalLattice::Pi, MomentumSpace::PiP);
		MPI_Barrier(MPI_COMM_WORLD);

		for (int pz = 0;pz < LocalLattice::N_eta_loc;pz++)
		{
			pzGlobal = pz + ID*LocalLattice::N_eta_loc;
			if (pzGlobal == 0)
				nu = 0;
			else
			{
				if (pzGlobal < N_eta / 2 + 1)
					nu = 2.0 / (a_eta)*sin((pzGlobal)*acos(-1.) / N_eta);
				else
					nu = -2.0 / (a_eta)*sin((pzGlobal)*acos(-1.) / N_eta);
			}

			for (int px = 0;px < N_t;px++)
			{
				if (px < N_t / 2 + 1)
					ppx = 2.0 / (a_t)*sin(px*acos(-1.) / N_t);
				else
					ppx = -2.0 / (a_t)*sin(px*acos(-1.) / N_t);

				for (int py = 0;py < N_t / 2 + 1; py++)
				{
					ppy = 2.0 / (a_t)*sin(py*acos(-1.) / N_t);

					pabs = sqrt(sqr(ppx) + sqr(ppy) + sqr(nu));

					//position of fields in momentum space
					posLOC = py + px*(N_t / 2 + 1) + pz*N_t*(N_t / 2 + 1);

					if (which_filter == 0)
					{
						PhiP[posLOC] *= W_tophat(r_filter, pabs);
						PiP[posLOC] *= W_tophat(r_filter, pabs);
					}
					if (which_filter == 1)
					{
						PhiP[posLOC] *= W_gaussian(r_filter, pabs);
						PiP[posLOC] *= W_gaussian(r_filter, pabs);
					}
					if (which_filter == 2)
                    {
						PhiP[posLOC] *= W_import_bounce( pabs);
						PiP[posLOC] *= W_import_bounce( pabs);
					}    
                    if (which_filter == 3)
                    {
						PhiP[posLOC] *= W_ktophat(r_filter, pabs);
						PiP[posLOC] *= W_ktophat(r_filter, pabs);
					}    
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		FFT::fftPtoX_direct(a, MomentumSpace::PhiP, LocalLattice::Phi);
		MPI_Barrier(MPI_COMM_WORLD);
		FFT::fftPtoX_direct(a, MomentumSpace::PiP, LocalLattice::Pi);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	return;
}
