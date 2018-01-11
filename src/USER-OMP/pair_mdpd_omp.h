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

/* ----------------------------------------------------------------------
   Contributing author: Yidong Xia (Idaho National Laboratory)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(mdpd/omp,PairMDPDOMP)

#else

#ifndef LMP_PAIR_MDPD_OMP_H
#define LMP_PAIR_MDPD_OMP_H

#include "pair_mdpd.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairMDPDOMP : public PairMDPD, public ThrOMP {

 public:
  PairMDPDOMP(class LAMMPS *);
  virtual ~PairMDPDOMP();

  virtual void compute(int, int);
  virtual double memory_usage();

 protected:
  class RanMars **random_thr;
  int nthreads;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair dpd requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

W: Pair dpd needs newton pair on for momentum conservation

Self-explanatory.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
