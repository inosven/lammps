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

FixStyle(mvv/dpd/omp,FixMvvDPDOMP)

#else

#ifndef LMP_FIX_MVV_DPD_OMP_H
#define LMP_FIX_MVV_DPD_OMP_H

#include "fix_mvv_dpd.h"

namespace LAMMPS_NS {

class FixMvvDPDOMP : public FixMvvDPD {
 public:
  FixMvvDPDOMP(class LAMMPS *, int, char **);
  virtual ~FixMvvDPDOMP() {}
  virtual void initial_integrate(int);
  virtual void final_integrate();
};

}

#endif
#endif
