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
   Contributing author: Yidong Xia

   Ref. P.B. Warren
        Vapour-liquid coexistence in many-body dissipative particle dynamics
        Phys. Rev. E 68, 066702 – Published 18 December 2003
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_mdpd.h"
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

PairMDPD::PairMDPD(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  random = NULL;

  nmax = 0;
  rho = NULL;

  // set comm size needed by this Pair

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

PairMDPD::~PairMDPD()
{
  memory->destroy(rho);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(gamma);
    memory->destroy(sigma);
    memory->destroy(A);
    memory->destroy(B);
    memory->destroy(rd);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairMDPD::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double wrho,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,randnum,factor_mdpd;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);
  double kbt = force->boltz*temperature;

  double PI = 4.0*atan(1.0);

  // grow local arrays if necessary, and need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  for (i = 0; i < nall; i++) rho[i] = 0.0;

  // first loop over neighbors of my atoms to calculate the density

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      r = sqrt(delx*delx + dely*dely + delz*delz);

      if (r < rd[itype][jtype]) {
        if (domain->dimension == 2) {
          error->all(FLERR,"Weight density function has not been implemented for 2D");
        }
        else {
          wrho = 7.5/PI/pow(rd[itype][jtype],3)*pow(1.0-r/rd[itype][jtype],2);
          rho[i] += wrho;
          if (newton_pair || j < nlocal) rho[j] += wrho;
        }
      }
    }
  }

  // communicate densities

  comm->forward_comm_pair(this);

  // second loop over neighbors of my atoms to calculate the force

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
      factor_mdpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx*delvx + dely*delvy + delz*delvz;
        wd = 1.0 - r/cut[itype][jtype];
        randnum = random->gaussian();

        // conservative force = Eq. (21) in Ref.
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        fpair = A[itype][jtype]*wd;
        if (r < rd[itype][jtype])
          fpair += B[itype][jtype]*(rho[i]+rho[j])*(1.0-r/rd[itype][jtype]);
        fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
        fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
        fpair *= factor_mdpd*rinv;

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

          evdwl = 1.0;
          evdwl *= factor_mdpd;
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

void PairMDPD::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut  ,n+1,n+1,"pair:cut"  );
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(A,n+1,n+1,"pair:A");
  memory->create(B,n+1,n+1,"pair:B");
  memory->create(rd,n+1,n+1,"pair:rd");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMDPD::settings(int narg, char **arg)
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

void PairMDPD::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double gamma_one = force->numeric(FLERR,arg[2]);
  double A_one = force->numeric(FLERR,arg[3]);
  double B_one = force->numeric(FLERR,arg[4]);
  double rd_one = force->numeric(FLERR,arg[5]);

  double cut_one = cut_global;
  if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
        gamma[i][j] = gamma_one;
            A[i][j] = A_one;
            B[i][j] = B_one;
           rd[i][j] = rd_one;
          cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMDPD::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair mdpd requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(FLERR,
      "Pair mdpd needs newton pair on for momentum conservation");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMDPD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

    cut[j][i] =   cut[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];
      A[j][i] =     A[i][j];
      B[j][i] =     B[i][j];
     rd[j][i] =    rd[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMDPD::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(    &A[i][j],sizeof(double),1,fp);
        fwrite(    &B[i][j],sizeof(double),1,fp);
        fwrite(   &rd[i][j],sizeof(double),1,fp);
        fwrite(  &cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMDPD::read_restart(FILE *fp)
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
          fread(&gamma[i][j],sizeof(double),1,fp);
          fread(    &A[i][j],sizeof(double),1,fp);
          fread(    &B[i][j],sizeof(double),1,fp);
          fread(   &rd[i][j],sizeof(double),1,fp);
          fread(  &cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(    &A[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(    &B[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(   &rd[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(  &cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMDPD::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global ,sizeof(double),1,fp);
  fwrite(&seed       ,sizeof(int)   ,1,fp);
  fwrite(&mix_flag   ,sizeof(int)   ,1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMDPD::read_restart_settings(FILE *fp)
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

void PairMDPD::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n", i, gamma[i][i],
                                          A[i][i],
                                          B[i][i],
                                         rd[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMDPD::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n", i, j, gamma[i][j],
                                                     A[i][j],
                                                     B[i][j],
                                                    rd[i][j],
                                                   cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairMDPD::single(int i, int j, int itype, int jtype, double rsq,
                       double factor_coul, double factor_mdpd, double &fforce)
{
  double r,rinv,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];
  //fforce = a0[itype][jtype]*wd * factor_mdpd*rinv;
  fforce = 0.0;

  //phi = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
  phi = 0.0;
  return factor_mdpd*phi;
}

/* ---------------------------------------------------------------------- */

int PairMDPD::pack_forward_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMDPD::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) rho[i] += buf[m++];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairMDPD::memory_usage()
{
  double bytes = Pair::memory_usage();
  bytes += nmax * sizeof(double);
  return bytes;
}
