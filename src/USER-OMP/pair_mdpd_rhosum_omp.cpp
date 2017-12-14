/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/*-----------------------------------------------------------------------
   Contributing author: Yidong Xia (Idaho National Laboratory)
------------------------------------------------------------------------- */

#include <math.h>
#include "pair_mdpd_rhosum_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "domain.h"
#include "random_mars.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMDPDRhoSumOMP::PairMDPDRhoSumOMP(LAMMPS *lmp) :
  PairMDPDRhoSum(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairMDPDRhoSumOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    eval(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
//    reduce_thr(this, eflag, vflag, thr); // no need to do this
  } // end of omp parallel region
}

void PairMDPDRhoSumOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double r, rsq, h, ih, ihsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double wf;

  const dbl3_t * _noalias const x    = (dbl3_t *) atom->x[0];
  const int    * _noalias const type = atom->type;
        double * _noalias const rho  = atom->rho;
  const double * _noalias const mass = atom->mass;

  // check consistency of pair coefficients
  if (first)
  {
    for (i = 1; i <= atom->ntypes; i++)
      for (j = 1; i <= atom->ntypes; i++)
        if (cutsq[i][j] > 0.0)
          if (!setflag[i][i] || !setflag[j][j])
            if (comm->me == 0)
              printf("mDPD particle types %d and %d interact, but not all of their single particle properties are set.\n", i, j);

    first = 0;
  }

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density with a full neighborlist here

  // initialize density with zero
  for (ii = iifrom; ii < iito; ++ii)
    rho[ilist[ii]] = 0;

  // add density at each atom via kernel function overlap
  for (ii = iifrom; ii < iito; ++ii)
  {
    i     = ilist[ii];
    xtmp  = x[i].x;
    ytmp  = x[i].y;
    ztmp  = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++)
    {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      delx  = xtmp - x[j].x;
      dely  = ytmp - x[j].y;
      delz  = ztmp - x[j].z;
      rsq   = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype])
      {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;

        // Lucy kernel, 3D
        r = sqrt(rsq);
        wf = (h - r) * ihsq;
        wf =  2.0889086280811262819e0 * (h + 3. * r) * wf * wf * wf * ih;
        rho[i] += mass[jtype] * wf;
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate densities
#if defined(_OPENMP)
#pragma omp master
#endif
  { comm->forward_comm_pair(this); }

  // wait until master thread is done with communication
  sync_threads();
}

/* ---------------------------------------------------------------------- */

double PairMDPDRhoSumOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMDPDRhoSum::memory_usage();

  return bytes;
}
