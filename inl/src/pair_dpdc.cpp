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
   Contributing author: Yidong Xia (Idaho National Laboratory)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_dpdc.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPDC::PairDPDC(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairDPDC::~PairDPDC()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(sigma);
    memory->destroy(A);
    memory->destroy(rc1);
    memory->destroy(B);
    memory->destroy(rc2);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPDC::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wc,wd,randnum,factor_dpdc;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpdc = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DPDC systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx*delvx + dely*delvy + delz*delvz;
        wd = 1.0 - r/cut[itype][jtype];
        randnum = random->gaussian();

        // cubic spline conservative force weight function
        wc = A[itype][jtype]*Wfunc(r,rc1[itype][jtype])
           - B[itype][jtype]*Wfunc(r,rc2[itype][jtype]);

        // conservative force = a0 * wc
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        fpair = a0[itype][jtype]*wc;
        fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
        fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
        fpair *= factor_dpdc*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
          evdwl *= factor_dpdc;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDPDC::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut  ,n+1,n+1,"pair:cut"  );
  memory->create(a0   ,n+1,n+1,"pair:a0"   );
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(A    ,n+1,n+1,"pair:A"    );
  memory->create(rc1  ,n+1,n+1,"pair:rc1"  );
  memory->create(B    ,n+1,n+1,"pair:B"    );
  memory->create(rc2  ,n+1,n+1,"pair:rc2"  );
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDPDC::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  temperature = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[1]);
  seed = force->inumeric(FLERR,arg[2]);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDC::coeff(int narg, char **arg)
{
  if (narg < 8 || narg > 9) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double    a0_one = force->numeric(FLERR,arg[2]);
  double gamma_one = force->numeric(FLERR,arg[3]);
  double     A_one = force->numeric(FLERR,arg[4]);
  double   rc1_one = force->numeric(FLERR,arg[5]);
  double     B_one = force->numeric(FLERR,arg[6]);
  double   rc2_one = force->numeric(FLERR,arg[7]);

  double cut_one = cut_global;
  if (narg == 9) cut_one = force->numeric(FLERR,arg[8]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
           a0[i][j] =    a0_one;
        gamma[i][j] = gamma_one;
            A[i][j] =     A_one;
          rc1[i][j] =   rc1_one;
            B[i][j] =     B_one;
          rc2[i][j] =   rc2_one;
          cut[i][j] =   cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDC::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpdc requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(FLERR,
      "Pair dpdc needs newton pair on for momentum conservation");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDPDC::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

    cut[j][i] =   cut[i][j];
     a0[j][i] =    a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];
      A[j][i] =     A[i][j];
    rc1[j][i] =   rc1[i][j];
      B[j][i] =     B[i][j];
    rc2[j][i] =   rc2[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDC::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(   &a0[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(    &A[i][j],sizeof(double),1,fp);
        fwrite(  &rc1[i][j],sizeof(double),1,fp);
        fwrite(    &B[i][j],sizeof(double),1,fp);
        fwrite(  &rc2[i][j],sizeof(double),1,fp);
        fwrite(  &cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDC::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(   &a0[i][j],sizeof(double),1,fp);
          fread(&gamma[i][j],sizeof(double),1,fp);
          fread(    &A[i][j],sizeof(double),1,fp);
          fread(  &rc1[i][j],sizeof(double),1,fp);
          fread(    &B[i][j],sizeof(double),1,fp);
          fread(  &rc2[i][j],sizeof(double),1,fp);
          fread(  &cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(   &a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(    &A[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(  &rc1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(    &B[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(  &rc2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(  &cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDC::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global ,sizeof(double),1,fp);
  fwrite(&seed       ,sizeof(int)   ,1,fp);
  fwrite(&mix_flag   ,sizeof(int)   ,1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDC::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&temperature,sizeof(double),1,fp);
    fread(&cut_global ,sizeof(double),1,fp);
    fread(&seed       ,sizeof(int)   ,1,fp);
    fread(&mix_flag   ,sizeof(int)   ,1,fp);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed       ,1,MPI_INT   ,0,world);
  MPI_Bcast(&mix_flag   ,1,MPI_INT   ,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDPDC::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g\n",          i,
                                           a0[i][i],
                                        gamma[i][i],
                                            A[i][i],
                                          rc1[i][i],
                                            B[i][i],
                                          rc2[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDPDC::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g\n",          i,
                                                          j,
                                                   a0[i][j],
                                                gamma[i][j],
                                                    A[i][j],
                                                  rc1[i][j],
                                                    B[i][j],
                                                  rc2[i][j],
                                                  cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairDPDC::single(int i, int j, int itype, int jtype, double rsq,
                       double factor_coul, double factor_dpdc, double &fforce)
{
  double r,rinv,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];
  fforce = a0[itype][jtype]*wd * factor_dpdc*rinv;

  phi = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
  return factor_dpdc*phi;
}

/* ---------------------------------------------------------------------- */

double PairDPDC::Wfunc(double r, double rc)
{
  double W;

  if (r < (rc/2.0))
    W = 12.0/pow(rc,2.0)*r - 18.0/pow(rc,3.0)*r*r;
  else if (r < rc)
    W = 6.0/rc*pow((1.0-r/rc),2.0);
  else
    W = 0.0;

  //return W/(0.75*PI*rc); // this can reproduce Poiseuille flow velocity profile
  //return W/(PI*rc);
  return W;
}
