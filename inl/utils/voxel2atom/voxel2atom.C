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
  bool   *  wall;
  double ** x;
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
  cout << "  An OpenMP parallel code to convert voxel map to LAMMPS data\n";
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
  unsigned int i,j,k,n,ii,jj,kk;
  unsigned int nx,ny,nz;
  unsigned int natoms = 0;
  unsigned int ntypes = 0;
  unsigned int nphase = 0;
  unsigned int nsolid = 0;
  unsigned int nfluid = 0;
  unsigned int nwalls = 0;
  unsigned int nsplit = 0;
  unsigned int npower = 0;
  unsigned int nneigh = 0;
  unsigned int matoms = 0;
  unsigned int mtypes = 3;
  unsigned int icount = 0;
  unsigned int imass  = 0;
  unsigned int ipair  = 0;
  unsigned int itype  = 0;
  unsigned int npoin  = 0;
  unsigned int iphase_solid = 0;
  unsigned int iphase_fluid = 1;
  unsigned int ip1,ip2,ip3,ip4;
  unsigned int nip1,nim1,njp1,njm1,nkp1,nkm1;
  int m;
  int xflag,yflag,zflag;

  // floats
  const double eps = 1e-14;
  double scale;
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

  bool ifound;

  // characters
  char cinit[80];

  // strings
  string skipl;

  // structs
  atoms_t atom_in;
  atoms_t atom_ou;

  // files
  ifstream finpu;
  ifstream finit;
  ofstream fpost;

  // ===========================================================================
  // Open configuration input file
  // ===========================================================================

  finpu.open("voxel2atom.inp", ios::in);
  if (!finpu.is_open()) {
    cout << "Error ... the file 'voxel2atom.inp' does not exist! " << endl;
    exit(0);
  }

  // input the type of particles to dump
  getline(finpu,skipl);
  finpu >> nphase; getline(finpu,skipl);
  if (!((nphase==iphase_solid)||(nphase==iphase_fluid))) {
    cout << "Error ... in file 'voxel2atom.inp', nphase = " << nphase << " is not a valid value" << endl;
    exit(0);
  }

  // input the number of neighboring voxcels in each direction for searching wall particles
  getline(finpu,skipl);
  finpu >> nneigh; getline(finpu,skipl);
  if (!((nneigh==1)||(nneigh==2)||(nneigh==3)||(nneigh==4))) {
    cout << "Error ... in file 'voxel2atom.inp', nneigh = " << nneigh << " is not a valid value" << endl;
    exit(0);
  }

  // input the number of split refinement for each voxcel
  getline(finpu,skipl);
  finpu >> nsplit; getline(finpu,skipl);
  if (!((nsplit==0)||(nsplit==1)||(nsplit==2)||(nsplit==3))) {
    cout << "Error ... in file 'voxel2atom.inp', nsplit = " << nsplit << " is not a valid value" << endl;
    exit(0);
  }
  npower = (int)(pow(2.0,(float)(nsplit)));

  // input the scaling factor
  getline(finpu,skipl);
  finpu >> scale; getline(finpu,skipl);
  if (scale <= 0.0) {
    cout << "Error ... in file 'voxel2atom.inp', scale = " << scale << " is not a valid value" << endl;
    exit(0);
  }

  finpu.close();

  // ===========================================================================
  // Open ASCII 3D colormap data file
  // ===========================================================================

  // input the file name
  do {
    cout << "Enter the name of ASCII 3D colormap data file:\n";
    cin >> cinit;
    finit.open(cinit, ios::in);
  } while (!finit.is_open());
  cout << "The ASCII 3D colormap data file is "  << cinit << endl << endl;

  // define Lattice
  getline(finit,skipl);
  finit >> nx >> ny >> nz; getline(finit,skipl);

  // bounding box
  getline(finit,skipl);
  finit >> xlo >> xhi >> ylo >> yhi >> zlo >> zhi; getline(finit,skipl);

  // types of atoms
  getline(finit,skipl);
  finit >> ntypes; getline(finit,skipl);

  // calculate the number of initial atoms
  natoms = nx * ny * nz;

  // create dynamic arrays
  atom_in.type = new int [natoms+1];
  atom_in.wall = new bool [natoms+1];

  // print header information
  cout << left << setw(ndeli) << "natoms" << setw(ndeli) << natoms << endl;
  cout << left << setw(ndeli) << "ntypes" << setw(ndeli) << ntypes << endl;
  cout << left << setw(ndeli) << "xlo"    << setw(ndeli) << xlo    << endl;
  cout << left << setw(ndeli) << "xhi"    << setw(ndeli) << xhi    << endl;
  cout << left << setw(ndeli) << "ylo"    << setw(ndeli) << ylo    << endl;
  cout << left << setw(ndeli) << "yhi"    << setw(ndeli) << yhi    << endl;
  cout << left << setw(ndeli) << "zlo"    << setw(ndeli) << zlo    << endl;
  cout << left << setw(ndeli) << "zhi"    << setw(ndeli) << zhi    << endl;
  cout << endl;

  // I*J*K lattice data
  getline(finit,skipl);
  for (n=1; n<=natoms; n++) {
    finit >> atom_in.type[n]; getline(finit,skipl);
    atom_in.wall[n] = false;
    if (atom_in.type[n] == iphase_fluid) nfluid ++;
    cout << "\r" << n << " of " << natoms << ", found " << nfluid << " pore voxcels" << flush;
  }
  cout << endl;
  cout << "Sample porosity is " << (float)(nfluid)/(float)(natoms) << endl;
  cout << endl;

  getline(finit,skipl);
  finit.close();
  cout << "Finished reading the ASCII 3D colormap data file\n" << endl;

  // calculate the number of solid particles
  nsolid = natoms - nfluid;

  // loop over all voxcels
  for (k=1; k<=nz; k++)
     for (j=1; j<=ny; j++)
       for (i=1; i<=nx; i++) {
         // current index
         n = (k-1)*ny*nx + (j-1)*nx + i;
         if (atom_in.type[n] == iphase_solid && atom_in.wall[n] == false) {
           ifound = false;
           // (2*nneigh+1)^3 -vexcel block to emcapsulate the fluid particles
           for (kk=k-nneigh; kk<=k+nneigh; kk++) {
             for (jj=j-nneigh; jj<=j+nneigh; jj++) {
               for (ii=i-nneigh; ii<=i+nneigh; ii++) {
                 // a coloring scheme to find the required 'wall' particles
                 m = (kk-1)*ny*nx + (jj-1)*nx + ii;
                 if (m < 1 || m > natoms) continue;
                 if (atom_in.type[m] == iphase_fluid) {
                   atom_in.wall[n] = true;
                   nwalls ++;
                   cout << "\r" << nwalls << " wall voxcels found" << flush;
                   ifound = true;
                   break;
                 }
               }
               if (ifound == true) break;
             }
             if (ifound == true) break;
           }
         }
       }
  cout << endl;
  // re-assign the total number of particles to be written in output
  nwalls = nwalls*npower*npower*npower;
  matoms = nwalls;

  if (nphase == iphase_fluid) {
    nfluid = nfluid*(npower-1)*(npower-1)*(npower-1);
    matoms+= nfluid;
  }

  // ===========================================================================
  // Write to the new LAMMPS data file
  // ===========================================================================

  xlo = 0.0; xhi = (float)(nx)*scale;
  ylo = 0.0; yhi = (float)(ny)*scale;
  zlo = 0.0; zhi = (float)(nz)*scale;

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

  icount = 0;

  for (k=1; k<=nz; k++)
    for (j=1; j<=ny; j++)
      for (i=1; i<=nx; i++) {
        n = (k-1)*ny*nx + (j-1)*nx + i;
        if (atom_in.wall[n] == true) {
          for (kk=1; kk<=npower; kk++)
            for (jj=1; jj<=npower; jj++)
              for (ii=1; ii<=npower; ii++) {
                icount ++;
                fpost << icount                     << " "
                      << atom_in.type[n]+1          << " "
                      << ((float)(i) - 1.0 + ((float)(ii)-0.5)/(float)(npower))*scale << " "
                      << ((float)(j) - 1.0 + ((float)(jj)-0.5)/(float)(npower))*scale << " "
                      << ((float)(k) - 1.0 + ((float)(kk)-0.5)/(float)(npower))*scale << " "
                      << endl;
                cout << "\r No. " << icount << " of " << matoms << " (solid) particles written in lammps_restart.dat" << flush;
              }
        }
      }
  cout << endl;

  if (nphase == iphase_fluid) {
    for (k=1; k<=nz; k++)
      for (j=1; j<=ny; j++)
        for (i=1; i<=nx; i++) {
          n = (k-1)*ny*nx + (j-1)*nx + i;
          if (atom_in.type[n] == iphase_fluid) {
            for (kk=1; kk<=(npower-1); kk++)
              for (jj=1; jj<=(npower-1); jj++)
                for (ii=1; ii<=(npower-1); ii++) {
                  icount ++;
                  fpost << icount                     << " "
                        << atom_in.type[n]+1          << " "
                        << ((float)(i) - 1.0 + ((float)(ii)-0.5)/(float)(npower-1))*scale << " "
                        << ((float)(j) - 1.0 + ((float)(jj)-0.5)/(float)(npower-1))*scale << " "
                        << ((float)(k) - 1.0 + ((float)(kk)-0.5)/(float)(npower-1))*scale << " "
                        << endl;
                  cout << "\r No. " << icount << " of " << matoms << " (fluid) particles written in lammps_restart.dat" << flush;
                }
          }
        }
  }
  cout << endl;

  fpost.close();

  // release dynamic arrays
  delete [] atom_in.type;
  delete [] atom_in.wall;

  cout << "Done processing the LAMMPS restart file...\n" << endl;
  exit(0);
}
