#!/bin/bash

/home/xiay/gitprojects/lammps/src/lmp_icc_omp                   -i in.mdpd_liquid_drop
mv log.lammps log.serial

/home/xiay/gitprojects/lammps/src/lmp_icc_omp -sf omp -pk omp 1 -i in.mdpd_liquid_drop
mv log.lammps log.omp01

/home/xiay/gitprojects/lammps/src/lmp_icc_omp -sf omp -pk omp 2 -i in.mdpd_liquid_drop
mv log.lammps log.omp02
