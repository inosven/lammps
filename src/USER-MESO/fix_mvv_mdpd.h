/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mvv/mdpd,FixMvvMDPD)

#else

#ifndef LMP_FIX_MVV_MDPD_H
#define LMP_FIX_MVV_MDPD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMvvMDPD : public Fix {
 public:
  FixMvvMDPD(class LAMMPS *, int, char **);
  virtual ~FixMvvMDPD() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void reset_dt();

 protected:
  double dtv, dtf;
  double verlet;
};

}

#endif
#endif
