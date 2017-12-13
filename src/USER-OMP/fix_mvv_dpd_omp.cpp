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

/* ----------------------------------------------------------------------
   This is a time integrator for position and velocity (x and v) using the
   modified velocity-Verlet (MVV) algorithm.
   Setting verlet = 0.5 recovers the standard velocity-Verlet algorithm.

   Contributing author: Yidong Xia (Idaho National Laboratory)
   Email: yidong.xia@inl.gov
------------------------------------------------------------------------- */

#include "fix_mvv_dpd_omp.h"
#include "atom.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

FixMvvDPDOMP::FixMvvDPDOMP(LAMMPS *lmp, int narg, char **arg) :
  FixMvvDPD(lmp, narg, arg) { }

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixMvvDPDOMP::initial_integrate(int vflag)
{
  // update v, x and vest of atoms in group

        dbl3_t * _noalias const x    = (dbl3_t *) atom->x[0];
        dbl3_t * _noalias const v    = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f    = (dbl3_t *) atom->f[0];
        dbl3_t * _noalias const vest = (dbl3_t *) atom->vest[0];

  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  int i;

  if (atom->rmass)
  {
    const double * const rmass = atom->rmass;

#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
      {
        const double dtfm = dtf / rmass[i];

        vest[i].x = v[i].x + dtfm * f[i].x;
        vest[i].y = v[i].y + dtfm * f[i].y;
        vest[i].z = v[i].z + dtfm * f[i].z;

        x[i].x += dtv * vest[i].x;
        x[i].y += dtv * vest[i].y;
        x[i].z += dtv * vest[i].z;
        v[i].x += 2.0 * verlet * dtfm * f[i].x;
        v[i].y += 2.0 * verlet * dtfm * f[i].y;
        v[i].z += 2.0 * verlet * dtfm * f[i].z;
      }
  }
  else
  {
    const double * const mass = atom->mass;
    const int *    const type = atom->type;

#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
      {
        const double dtfm = dtf / mass[type[i]];

        vest[i].x = v[i].x + dtfm * f[i].x;
        vest[i].y = v[i].y + dtfm * f[i].y;
        vest[i].z = v[i].z + dtfm * f[i].z;

        x[i].x += dtv * vest[i].x;
        x[i].y += dtv * vest[i].y;
        x[i].z += dtv * vest[i].z;
        v[i].x += 2.0 * verlet * dtfm * f[i].x;
        v[i].y += 2.0 * verlet * dtfm * f[i].y;
        v[i].z += 2.0 * verlet * dtfm * f[i].z;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvDPDOMP::final_integrate()
{
        dbl3_t * _noalias const v    = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f    = (dbl3_t *) atom->f[0];
  const dbl3_t * _noalias const vest = (dbl3_t *) atom->vest[0];

  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  int i;

  if (atom->rmass)
  {
    const double * const rmass = atom->rmass;

#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
      {
        const double dtfm = dtf / rmass[i];

        v[i].x = vest[i].x + dtfm * f[i].x;
        v[i].y = vest[i].y + dtfm * f[i].y;
        v[i].z = vest[i].z + dtfm * f[i].z;
      }
  }
  else
  {
    const double * const mass = atom->mass;
    const int *    const type = atom->type;

#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
      {
        const double dtfm = dtf / mass[type[i]];

        v[i].x = vest[i].x + dtfm * f[i].x;
        v[i].y = vest[i].y + dtfm * f[i].y;
        v[i].z = vest[i].z + dtfm * f[i].z;
      }
  }
}
