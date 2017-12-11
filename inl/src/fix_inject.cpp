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
#include "stdlib.h"
#include "string.h"
#include "fix_inject.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_park.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{ATOM,MOLECULE};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

#define EPSILON 1.0e6

/* ---------------------------------------------------------------------- */

FixInject::FixInject(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix inject command");

  restart_global = 1;
  time_depend = 1;

  // required args

  ninject = force->inumeric(FLERR,arg[3]);
  ntype   = force->inumeric(FLERR,arg[4]);
  nfreq   = force->inumeric(FLERR,arg[5]);
  natom   = force->inumeric(FLERR,arg[6]);
  seed    = force->inumeric(FLERR,arg[7]);

  if (seed <= 0) error->all(FLERR,"Illegal fix inject command");

  // read options from end of input line

  options(narg-8,&arg[8]);

  // error check on type

  if (mode == ATOM && (ntype <= 0 || ntype > atom->ntypes))
    error->all(FLERR,"Invalid atom type in fix inject command");

  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all(FLERR,"Must specify a region in fix inject");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Fix inject region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all(FLERR,"Fix inject region cannot be dynamic");

  xlo = domain->regions[iregion]->extent_xlo;
  xhi = domain->regions[iregion]->extent_xhi;
  ylo = domain->regions[iregion]->extent_ylo;
  yhi = domain->regions[iregion]->extent_yhi;
  zlo = domain->regions[iregion]->extent_zlo;
  zhi = domain->regions[iregion]->extent_zhi;

  if (domain->triclinic == 0) {
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
        ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
        zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Injection region extends outside simulation box");
  } else {
    if (xlo < domain->boxlo_bound[0] || xhi > domain->boxhi_bound[0] ||
        ylo < domain->boxlo_bound[1] || yhi > domain->boxhi_bound[1] ||
        zlo < domain->boxlo_bound[2] || zhi > domain->boxhi_bound[2])
      error->all(FLERR,"Injection region extends outside simulation box");
  }

  // error check and further setup for mode = MOLECULE

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix_inject unless atoms have IDs");

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling to all input parameters with dist/vel units

  if (domain->dimension == 2) {
    lo *= yscale;
    hi *= yscale;
  } else {
    lo *= zscale;
    hi *= zscale;
  }
  deltasq *= xscale*xscale;
  nearsq *= xscale*xscale;
  vxlo *= xscale;
  vxhi *= xscale;
  vylo *= yscale;
  vyhi *= yscale;
  vzlo *= zscale;
  vzhi *= zscale;
  tx *= xscale;
  ty *= yscale;
  tz *= zscale;

  // find current max atom and molecule IDs if necessary

  if (idnext) find_maxid();

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninjected = 0;
}

/* ---------------------------------------------------------------------- */

FixInject::~FixInject()
{
  delete random;
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixInject::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInject::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix inject does not exist");

  // for finite size spherical particles:
  // warn if near < 2 * maxrad of existing and injected particles
  //   since may lead to overlaps
  // if injected molecule does not define diameters,
  //   use AtomVecSphere::create_atom() default radius = 0.5

  if (atom->radius_flag) {
    double *radius = atom->radius;
    int nlocal = atom->nlocal;

    double maxrad = 0.0;
    for (int i = 0; i < nlocal; i++)
      maxrad = MAX(maxrad,radius[i]);

    double maxradall;
    MPI_Allreduce(&maxrad,&maxradall,1,MPI_DOUBLE,MPI_MAX,world);

    double maxradinject = 0.0;
    maxradinject = 0.5;

    double separation = MAX(2.0*maxradinject,maxradall+maxradinject);
    if (sqrt(nearsq) < separation && comm->me == 0) {
      char str[128];
      sprintf(str,"Fix inject near setting < possible overlap separation %g",
              separation);
      error->warning(FLERR,str);
    }
  }
}

/* ----------------------------------------------------------------------
   perform particle injection
------------------------------------------------------------------------- */

void FixInject::pre_exchange()
{
  int i,m,n,nlocalprev,flag,flagall,success,attempt;
  imageint imageflag;
  double coord[3],lamda[3],delx,dely,delz,rsq;
  double r[3],vnew[3],rotmat[3][3],quat[4];
  double *newcoord;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  for (m = 0; m < natom; m++) {
    // find current max atom and molecule IDs if necessary

    if (!idnext) find_maxid();

    // attempt an injection until successful

    success = 0;
    attempt = 0;
    while (attempt < maxattempt) {
      attempt++;

      // choose random position for new particle within region

      coord[0] = xlo + random->uniform() * (xhi-xlo);
      coord[1] = ylo + random->uniform() * (yhi-ylo);
      coord[2] = zlo + random->uniform() * (zhi-zlo);
      while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
        coord[0] = xlo + random->uniform() * (xhi-xlo);
        coord[1] = ylo + random->uniform() * (yhi-ylo);
        coord[2] = zlo + random->uniform() * (zhi-zlo);
      }

      // apply PBC so final coords are inside box
      // also modify image flags due to PBC

      imageflag = ((imageint) IMGMAX << IMG2BITS) |
                  ((imageint) IMGMAX << IMGBITS)  | IMGMAX;

      // check distance between any existing atom and any injected atom
      // if less than near, try again
      // use minimum_image() to account for PBC

      flag = 0;
      for (i = 0; i < atom->nlocal; i++) {
        delx = coord[0] - atom->x[i][0];
        dely = coord[1] - atom->x[i][1];
        delz = coord[2] - atom->x[i][2];
        domain->minimum_image(delx,dely,delz);
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < nearsq) flag = 1;
      }
      MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
      if (flagall) continue;

      // proceed with injection

      // choose random velocity for new particle

      vnew[0] = vxlo + random->uniform() * (vxhi-vxlo);
      vnew[1] = vylo + random->uniform() * (vyhi-vylo);
      vnew[2] = vzlo + random->uniform() * (vzhi-vzlo);

      // check if new atoms are in my sub-box or above it if I am highest proc
      // if so, add atom to my list via create_atom()
      // initialize additional info about the atoms
      // set group mask to "all" plus fix group

      if (domain->triclinic) {
        domain->x2lamda(coord,lamda);
        newcoord = lamda;
      } else newcoord = coord;

      flag = 0;
      if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
          newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
          newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
      else if (domain->dimension == 3 && newcoord[2] >= domain->boxhi[2]) {
        if (comm->layout != LAYOUT_TILED) {
          if (comm->myloc[2] == comm->procgrid[2]-1 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
              newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
        } else {
          if (comm->mysplit[2][1] == 1.0 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
              newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
        }
      } else if (domain->dimension == 2 && newcoord[1] >= domain->boxhi[1]) {
        if (comm->layout != LAYOUT_TILED) {
          if (comm->myloc[1] == comm->procgrid[1]-1 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
        } else {
          if (comm->mysplit[1][1] == 1.0 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
        }
      }

      if (flag) {
        atom->avec->create_atom(ntype,coord);
        n = atom->nlocal - 1;
        atom->tag[n] = maxtag_all + m+1;
        atom->mask[n] = 1 | groupbit;
        atom->image[n] = imageflag;
        atom->v[n][0] = vnew[0];
        atom->v[n][1] = vnew[1];
        atom->v[n][2] = vnew[2];
        modify->create_attribute(n);
      }

      success = 1;
      break;
    }

    // warn if not successful b/c too many attempts

    if (!success && comm->me == 0)
      error->warning(FLERR,"Particle injection was unsuccessful",0);

    // reset global natoms,nbonds,etc
    // increment maxtag_all and maxmol_all if necessary
    // if global map exists, reset it now instead of waiting for comm
    // since adding atoms messes up ghosts

    if (success) {
      atom->natoms ++;
      if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
        error->all(FLERR,"Too many total atoms");
      maxtag_all ++;
      if (maxtag_all >= MAXTAGINT)
        error->all(FLERR,"New atom IDs exceed maximum allowed ID");
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
      ninjected++;
    }

    if (ninjected >= ninject) break;
  }

  // next timestep to inject
  // next_reneighbor = 0 if done

  if (ninjected < ninject) next_reneighbor += nfreq;
  else next_reneighbor = 0;
}

/* ----------------------------------------------------------------------
   maxtag_all = current max atom ID for all atoms
   maxmol_all = current max molecule ID for all atoms
------------------------------------------------------------------------- */

void FixInject::find_maxid()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixInject::options(int narg, char **arg)
{
  // defaults

  iregion = -1;
  idregion = NULL;
  mode = ATOM;
  idnext = 0;
  lo = hi = deltasq = 0.0;
  nearsq = 0.0;
  maxattempt = 10;
  vxlo = vxhi = vylo = vyhi = vzlo = vzhi = 0.0;
  scaleflag = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inject command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix inject does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"id") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inject command");
      if (strcmp(arg[iarg+1],"max") == 0) idnext = 0;
      else if (strcmp(arg[iarg+1],"next") == 0) idnext = 1;
      else error->all(FLERR,"Illegal fix inject command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"near") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inject command");
      nearsq = force->numeric(FLERR,arg[iarg+1]) *
        force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"attempt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inject command");
      maxattempt = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix inject command");
      vxlo = force->numeric(FLERR,arg[iarg+1]);
      vxhi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix inject command");
      vylo = force->numeric(FLERR,arg[iarg+1]);
      vyhi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix inject command");
      vzlo = force->numeric(FLERR,arg[iarg+1]);
      vzhi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inject command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix inject command");
      iarg += 2;

    } else error->all(FLERR,"Illegal fix inject command");
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixInject::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = random->state();
  list[n++] = ninjected;
  list[n++] = nfirst;
  list[n++] = next_reneighbor;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixInject::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  ninjected = static_cast<int> (list[n++]);
  nfirst = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);

  random->reset(seed);
}

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

void *FixInject::extract(const char *str, int &itype)
{
  if (strcmp(str,"radius") == 0) {
    if (mode == ATOM) {
      if (itype == ntype) oneradius = 0.5;
      else oneradius = 0.0;
    }
    itype = 0;
    return &oneradius;
  }

  return NULL;
}
