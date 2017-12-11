#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <omp.h>

using namespace std;

struct atoms_t {
  int    *  type;
  int    ** flag;
  double ** x;
};

struct tet4_t {
  double volMin,volMax;
  double *  volume;
  int    ** vertex;
};

struct coord_t {
  double xmin,xmax,ymin,ymax,zmin,zmax;
  double ** x;
};

double computeVolTet4(double xp1, double yp1, double zp1,
                      double xp2, double yp2, double zp2,
                      double xp3, double yp3, double zp3,
                      double xp4, double yp4, double zp4)
{
  double x21 = xp2 - xp1;
  double y21 = yp2 - yp1;
  double z21 = zp2 - zp1;

  double x31 = xp3 - xp1;
  double y31 = yp3 - yp1;
  double z31 = zp3 - zp1;

  double x41 = xp4 - xp1;
  double y41 = yp4 - yp1;
  double z41 = zp4 - zp1;

  return (x21*(y31*z41-z31*y41) +
          x31*(z21*y41-y21*z41) +
          x41*(y21*z31-z21*y31) )/6.0;
};



int main(int argc, char **argv)
{
  // ===========================================================================
  // A header message
  // ===========================================================================

  cout << endl;
  cout << "################################################################\n";
  cout << endl;
  cout << endl;
  cout << "  An OpenMP Parallel Code to Color Solid Grains in Shale Rocks\n";
  cout << endl;
  cout << "                          Yidong Xia \n";
  cout << endl;
  cout << "################################################################\n";
  cout << endl;

  // ===========================================================================
  // define data types
  // ===========================================================================

  // integers
  unsigned const int ndims  = 3;
  unsigned const int ndeli  = 18;
  unsigned const int nvert  = 4;
  unsigned int nthreads, tid;
  unsigned int i,j;
  unsigned int natoms = 0;
  unsigned int ntypes = 0;
  unsigned int nsolid = 0;
  unsigned int nfluid = 0;
  unsigned int matoms = 0;
  unsigned int mtypes = 0;
  unsigned int imass  = 0;
  unsigned int ipair  = 0;
  unsigned int iatom  = 0;
  unsigned int itype  = 0;
  unsigned int ntet4  = 0;
  unsigned int npoin  = 0;
  unsigned int ip1,ip2,ip3,ip4;
  int xflag,yflag,zflag;

  // floats
  const double eps = 1e-14;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xat,yat,zat;
  double xvelo,yvelo,zvelo;
  double vmass;
  double coef1,coef2;
  double vol;
  double xp1,yp1,zp1;
  double xp2,yp2,zp2;
  double xp3,yp3,zp3;
  double xp4,yp4,zp4;
  double volMin,volMax;

  // characters
  char cinit[80];
  char ctet4[80];

  // strings
  string skipl;

  // structs
  atoms_t atom_in;
  atoms_t atom_ou;

  tet4_t tet4_in;

  coord_t coord_in;

  // files
  ifstream finit;
  ifstream ftet4;
  ofstream fpost;

  // ===========================================================================
  // Open LAMMPS written data
  // ===========================================================================

  // type in file name
  do {
    cout << "Enter the name of LAMMPS written data file:\n";
    cin >> cinit;
    finit.open(cinit, ios::in);
  } while (!finit.is_open());
  cout << "The LAMMPS dump file is "  << cinit << endl << endl;

  // skip the first line
  getline(finit,skipl);

  // header information
  getline(finit,skipl);
  finit >> natoms; getline(finit,skipl);
  finit >> ntypes; getline(finit,skipl);
  cout << left << setw(ndeli) << "natoms" << setw(ndeli) << natoms << endl;
  cout << left << setw(ndeli) << "ntypes" << setw(ndeli) << ntypes << endl;
  cout << endl;

  // create dynamic arrays
  //atom = new atoms_t [natoms];

  atom_in.type = new int     [natoms];
  atom_in.flag = new int*    [natoms];
  atom_in.x    = new double* [natoms];
  for (unsigned int n=0; n<natoms; n++) {
    atom_in.flag[n] = new int    [ndims];
    atom_in.x[n]    = new double [ndims];
  }

  // box range
  getline(finit,skipl);
  finit >> xlo >> xhi; getline(finit,skipl);
  finit >> ylo >> yhi; getline(finit,skipl);
  finit >> zlo >> zhi; getline(finit,skipl);
  cout << left << setw(ndeli) << "xlo" << setw(ndeli) << xlo << endl;
  cout << left << setw(ndeli) << "xhi" << setw(ndeli) << xhi << endl;
  cout << left << setw(ndeli) << "ylo" << setw(ndeli) << ylo << endl;
  cout << left << setw(ndeli) << "yhi" << setw(ndeli) << yhi << endl;
  cout << left << setw(ndeli) << "zlo" << setw(ndeli) << zlo << endl;
  cout << left << setw(ndeli) << "zhi" << setw(ndeli) << zhi << endl;
  cout << endl;

  // masses
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);
  finit >> imass >> vmass; getline(finit,skipl);
  cout << left << setw(ndeli) << "mass ID"    << setw(ndeli) << imass << endl;
  cout << left << setw(ndeli) << "mass value" << setw(ndeli) << vmass << endl;
  cout << endl;

  // pair coefficients
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);
  finit >> ipair >> coef1 >> coef2; getline(finit,skipl);
  cout << left << setw(ndeli) << "pair ID"    << setw(ndeli) << ipair << endl;
  cout << left << setw(ndeli) << "pair coef1" << setw(ndeli) << coef1 << endl;
  cout << left << setw(ndeli) << "pair coef2" << setw(ndeli) << coef2 << endl;
  cout << endl;

  // list of atom coordinates
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);
  for (unsigned int n=0; n<natoms; n++) {
    finit >> iatom >> itype >> xat >> yat >> zat >> xflag >> yflag >> zflag;
    getline(finit,skipl);
    atom_in.type[iatom-1] = itype;
    atom_in.flag[iatom-1][0] = xflag;
    atom_in.flag[iatom-1][1] = yflag;
    atom_in.flag[iatom-1][2] = zflag;
    atom_in.x[iatom-1][0] = xat;
    atom_in.x[iatom-1][1] = yat;
    atom_in.x[iatom-1][2] = zat;
  }

  getline(finit,skipl);
  finit.close();

  cout << "Finished reading the particle data file\n" << endl;

  // ===========================================================================
  // Open Tecplot ASCII block type data file
  // ===========================================================================

  /*
      #-of-nodes #-of-tetrahedron
      x-coordinates
      y-coordinates
      z-coordinates
      list-of-connectivity
  */

  // type in file name
  do {
    cout << "Enter the name of Tecplot ASCII block-type data file:\n";
    cin >> ctet4;
    ftet4.open(ctet4, ios::in);
  } while (!ftet4.is_open());
  cout << "The Tecplot ASCII block-type data file is "  << ctet4 << endl;

  ftet4 >> npoin >> ntet4; getline(ftet4,skipl);

  // print header information
  cout << left << setw(ndeli) << "npoin" << setw(ndeli) << npoin << endl;
  cout << left << setw(ndeli) << "ntet4" << setw(ndeli) << ntet4 << endl;

  // create dynamic arrays
  //tet4_in = new tet4_t [ntet4];
  //coord_in = new coord_t [npoin];

  tet4_in.volume = new double [ntet4];
  tet4_in.vertex = new int*   [ntet4];
  for (unsigned int n=0; n<ntet4; n++)
    tet4_in.vertex[n] = new int [nvert];

  coord_in.x = new double* [npoin];
  for (unsigned int n=0; n<npoin; n++)
    coord_in.x[n] = new double [ndims];

  // read in point coordinates
  for (unsigned int n=0; n<npoin; n++)
    ftet4 >> coord_in.x[n][0];
  getline(ftet4,skipl);
  for (unsigned int n=0; n<npoin; n++)
    ftet4 >> coord_in.x[n][1];
  getline(ftet4,skipl);
  for (unsigned int n=0; n<npoin; n++)
    ftet4 >> coord_in.x[n][2];
  getline(ftet4,skipl);

  // read in element-vertex connectivity
  cout << "Read in element-vertex connectivity per tet4\n";
  for (unsigned int n=0; n<ntet4; n++) {
    ftet4 >> tet4_in.vertex[n][0]
          >> tet4_in.vertex[n][1]
          >> tet4_in.vertex[n][2]
          >> tet4_in.vertex[n][3];
    getline(ftet4,skipl);
    cout << "\r" << n+1 << " of " << ntet4 << flush;
  }
  cout << endl;

  // close the file if it is the end
  getline(ftet4,skipl);
  if(!ftet4.eof()) cout << "Not the end of file??????"  << endl;
  ftet4.close();

  // calculate volume per tet4
  for (unsigned int n=0; n<ntet4; n++) {
    ip1 = tet4_in.vertex[n][0];
    ip2 = tet4_in.vertex[n][1];
    ip3 = tet4_in.vertex[n][2];
    ip4 = tet4_in.vertex[n][3];

    xp1 = coord_in.x[ip1-1][0];
    yp1 = coord_in.x[ip1-1][1];
    zp1 = coord_in.x[ip1-1][2];

    xp2 = coord_in.x[ip2-1][0];
    yp2 = coord_in.x[ip2-1][1];
    zp2 = coord_in.x[ip2-1][2];

    xp3 = coord_in.x[ip3-1][0];
    yp3 = coord_in.x[ip3-1][1];
    zp3 = coord_in.x[ip3-1][2];

    xp4 = coord_in.x[ip4-1][0];
    yp4 = coord_in.x[ip4-1][1];
    zp4 = coord_in.x[ip4-1][2];

    vol = computeVolTet4(xp1, yp1, zp1,
                         xp2, yp2, zp2,
                         xp3, yp3, zp3,
                         xp4, yp4, zp4);

    if (vol < eps)
      cout << "Negative volume! No. of tet4 = " << n+1 << ", Vol = " << vol << endl;

    tet4_in.volume[n] = fabs(vol);
  }

  // Find the min and max volume of the tet4 mesh
  tet4_in.volMin = tet4_in.volume[0];
  tet4_in.volMax = tet4_in.volume[0];
  for (unsigned int n=1; n<ntet4; n++) {
    if (tet4_in.volume[n] < tet4_in.volMin) tet4_in.volMin = tet4_in.volume[n];
    if (tet4_in.volume[n] > tet4_in.volMax) tet4_in.volMax = tet4_in.volume[n];
  }
  cout << "Original tet4 min volume = " << tet4_in.volMin << endl;
  cout << "Original tet4 max volume = " << tet4_in.volMax << endl;

  // Find the min and max coordinates of the tet4 mesh
  coord_in.xmin = coord_in.x[0][0];
  coord_in.xmax = coord_in.x[0][0];
  coord_in.ymin = coord_in.x[0][1];
  coord_in.ymax = coord_in.x[0][1];
  coord_in.zmin = coord_in.x[0][2];
  coord_in.zmax = coord_in.x[0][2];

  for (unsigned int n=1; n<npoin; n++) {
    if (coord_in.x[n][0] < coord_in.xmin) coord_in.xmin = coord_in.x[n][0];
    if (coord_in.x[n][0] > coord_in.xmax) coord_in.xmax = coord_in.x[n][0];
    if (coord_in.x[n][1] < coord_in.ymin) coord_in.ymin = coord_in.x[n][1];
    if (coord_in.x[n][1] > coord_in.ymax) coord_in.ymax = coord_in.x[n][1];
    if (coord_in.x[n][2] < coord_in.zmin) coord_in.zmin = coord_in.x[n][2];
    if (coord_in.x[n][2] > coord_in.zmax) coord_in.zmax = coord_in.x[n][2];
  }
  cout << "Original mesh x bound = " << coord_in.xmin << ", " << coord_in.xmax << endl;
  cout << "Original mesh y bound = " << coord_in.ymin << ", " << coord_in.ymax << endl;
  cout << "Original mesh z bound = " << coord_in.zmin << ", " << coord_in.zmax << endl;

  // Scale the mesh to match LAMMPS data
  double xhilo = xhi - xlo;
  double yhilo = yhi - ylo;
  double zhilo = zhi - zlo;
  double xmaxmin = coord_in.xmax - coord_in.xmin;
  double ymaxmin = coord_in.ymax - coord_in.ymin;
  double zmaxmin = coord_in.zmax - coord_in.zmin;
  for (unsigned int n=0; n<npoin; n++) {
    coord_in.x[n][0] = (coord_in.x[n][0] - coord_in.xmin)*xhilo/xmaxmin + xlo;
    coord_in.x[n][1] = (coord_in.x[n][1] - coord_in.ymin)*yhilo/ymaxmin + ylo;
    coord_in.x[n][2] = (coord_in.x[n][2] - coord_in.zmin)*zhilo/zmaxmin + zlo;
  }

  // Find the min and max coordinates of the tet4 mesh after rescaling
  coord_in.xmin = coord_in.x[0][0];
  coord_in.xmax = coord_in.x[0][0];
  coord_in.ymin = coord_in.x[0][1];
  coord_in.ymax = coord_in.x[0][1];
  coord_in.zmin = coord_in.x[0][2];
  coord_in.zmax = coord_in.x[0][2];

  for (unsigned int n=1; n<npoin; n++) {
    if (coord_in.x[n][0] < coord_in.xmin) coord_in.xmin = coord_in.x[n][0];
    if (coord_in.x[n][0] > coord_in.xmax) coord_in.xmax = coord_in.x[n][0];
    if (coord_in.x[n][1] < coord_in.ymin) coord_in.ymin = coord_in.x[n][1];
    if (coord_in.x[n][1] > coord_in.ymax) coord_in.ymax = coord_in.x[n][1];
    if (coord_in.x[n][2] < coord_in.zmin) coord_in.zmin = coord_in.x[n][2];
    if (coord_in.x[n][2] > coord_in.zmax) coord_in.zmax = coord_in.x[n][2];
  }
  cout << "Rescaled mesh x bound = " << coord_in.xmin << ", " << coord_in.xmax << endl;
  cout << "Rescaled mesh y bound = " << coord_in.ymin << ", " << coord_in.ymax << endl;
  cout << "Rescaled mesh z bound = " << coord_in.zmin << ", " << coord_in.zmax << endl;

  // ===========================================================================
  // Color the atoms that are in pore network according to the mesh info
  // ===========================================================================

#if defined _OPENMP
  #pragma omp parallel private(tid)
  {
   tid = omp_get_thread_num();
   if (tid == 0)
     cout << "OpenMP mode, number of threads = " << omp_get_num_threads() << endl;
  }
#endif
  cout << "Loop over all atoms and color those in pore networks\n";

  // setup wall-clock timer
#if defined _OPENMP
  double wtime = omp_get_wtime( );
#else
  clock_t wtime = clock();
#endif

#pragma omp parallel shared(nfluid)
{
  #pragma omp for private(j,xat,yat,zat,ip1,ip2,ip3,ip4,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3,xp4,yp4,zp4)
  for (i=0; i<natoms; i++) {
    xat = atom_in.x[i][0];
    yat = atom_in.x[i][1];
    zat = atom_in.x[i][2];

    for (j=0; j<ntet4; j++) {
      ip1 = tet4_in.vertex[j][0];
      ip2 = tet4_in.vertex[j][1];
      ip3 = tet4_in.vertex[j][2];
      ip4 = tet4_in.vertex[j][3];

      xp1 = coord_in.x[ip1-1][0];
      yp1 = coord_in.x[ip1-1][1];
      zp1 = coord_in.x[ip1-1][2];

      xp2 = coord_in.x[ip2-1][0];
      yp2 = coord_in.x[ip2-1][1];
      zp2 = coord_in.x[ip2-1][2];

      xp3 = coord_in.x[ip3-1][0];
      yp3 = coord_in.x[ip3-1][1];
      zp3 = coord_in.x[ip3-1][2];

      xp4 = coord_in.x[ip4-1][0];
      yp4 = coord_in.x[ip4-1][1];
      zp4 = coord_in.x[ip4-1][2];

      if (computeVolTet4(xp2,yp2,zp2,xp4,yp4,zp4,xp3,yp3,zp3,xat,yat,zat) < 0.0) continue;
      if (computeVolTet4(xp1,yp1,zp1,xp3,yp3,zp3,xp4,yp4,zp4,xat,yat,zat) < 0.0) continue;
      if (computeVolTet4(xp1,yp1,zp1,xp4,yp4,zp4,xp2,yp2,zp2,xat,yat,zat) < 0.0) continue;
      if (computeVolTet4(xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3,xat,yat,zat) < 0.0) continue;

      atom_in.type[i] ++;
      #pragma omp atomic
      nfluid ++;
      break;
    }
#if defined _OPENMP
    if (omp_get_thread_num() == 0)
      cout << "\r" << "On master thread: " << i+1 << " of " << natoms << " particles, found "
           << nfluid << " pore particles" << flush;
#else
    cout << "\r" << i+1 << " of " << natoms << ", found " << nfluid << " pore particles" << flush;
#endif
  }
}
  cout << endl;

  // print elapsed wall-clock time
#if defined _OPENMP
  wtime = omp_get_wtime() - wtime;
  cout << "Execution time for coloring = " << wtime << " (seconds)" << endl;
#else
  wtime = clock() - wtime;
  cout << "Execution time for coloring = " << ((float)wtime)/CLOCKS_PER_SEC << " (seconds)" << endl;
#endif

  nsolid = natoms - nfluid;
  cout << "Result: found nfluid = " << nfluid << ", nsolid = " << nsolid << endl;



  matoms = natoms;
  mtypes = 2;

  // ===========================================================================
  // Write to the new LAMMPS data file
  // ===========================================================================

  fpost.open("lammps_restart.dat", ios::out);

  fpost.precision(16);
  fpost << scientific;

  // first line
  fpost << "LAMMPS data file for read_data" << endl;

  // head information block
  fpost << endl;
  fpost << matoms << " atoms" << endl;
  fpost << mtypes << " atom types" << endl;

  // range box block
  fpost << endl;
  fpost << xlo <<  " " << xhi << " xlo xhi" << endl;
  fpost << ylo <<  " " << yhi << " ylo yhi" << endl;
  fpost << zlo <<  " " << zhi << " zlo zhi" << endl;

  // coordinates
  fpost << endl;
  fpost << "Atoms" << endl;
  fpost << endl;
  for (unsigned int n=0; n<natoms; n++) {

    fpost << n+1 << " "
          << atom_in.type[n] << " "
          << atom_in.x[n][0] << " "
          << atom_in.x[n][1] << " "
          << atom_in.x[n][2] << " "
          << atom_in.flag[n][0] << " "
          << atom_in.flag[n][1] << " "
          << atom_in.flag[n][2] << endl;
  }
  fpost.close();
  cout << endl;
  cout << "Done processing the LAMMPS restart file..." << endl;


  // release dynamic arrays
  delete [] atom_in.type;
  delete [] atom_in.flag;
  delete [] atom_in.x;
  delete [] tet4_in.volume;
  delete [] tet4_in.vertex;
  delete [] coord_in.x;

  cout << "Exit program." << endl;
  cout << endl;
}
