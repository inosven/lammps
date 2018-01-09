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

   Contributing authors:
   Zhen Li (Brown University). Email: zhen_li@brown.edu
   Yidong Xia (Idaho National Laboratory). Email: yidong.xia@inl.gov
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "fix_mvv_mdpd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMvvMDPD::FixMvvMDPD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"mvv/mdpd") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix mvv/mdpd command");

  verlet = 0.5;
  if(narg > 3) verlet = force->numeric(FLERR,arg[3]);

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixMvvMDPD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMvvMDPD::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixMvvMDPD::initial_integrate(int vflag)
{
  double dtfm;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
  if (mask[i] & groupbit) {
    dtfm = dtf / mass[type[i]];
    for(int k = 0; k < 3; k++){
      vest[i][k] = v[i][k];
      vest[i][k+3] = f[i][k];
      x[i][k] += dtv * (v[i][k] + dtfm * f[i][k]);
      v[i][k] += 2.0 * verlet * dtfm * f[i][k];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvMDPD::final_integrate()
{
  double dtfm;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
  if (mask[i] & groupbit) {
    dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * ((1.0-2.0*verlet)*vest[i][3] + f[i][0]);
    v[i][1] += dtfm * ((1.0-2.0*verlet)*vest[i][4] + f[i][1]);
    v[i][2] += dtfm * ((1.0-2.0*verlet)*vest[i][5] + f[i][2]);
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvMDPD::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
