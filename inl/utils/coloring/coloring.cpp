#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

struct atoms_t {
  int    *  type;
  int    ** flag;
  double ** x;
  double ** v;
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
  cout << "    A Pre-Process Code to Color Solid Grains in Porous Media\n";
  cout << endl;
  cout << "            Yidong Xia (Idaho National Laboratory)\n";
  cout << endl;
  cout << "################################################################\n";
  cout << endl;

  // ===========================================================================
  // define data types
  // ===========================================================================

  // integers
  unsigned const int ndims  =  3;
  unsigned const int ndeli  = 18;
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
  int xflag,yflag,zflag;

  // floats
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xcoor,ycoor,zcoor;
  double xvelo,yvelo,zvelo;
  double vmass;
  double coef1,coef2;

  // characters
  char cinit[80];

  // strings
  string skipl;

  // structs
  atoms_t atom;

  // files
  ifstream finit;
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

  atom.type = new int     [natoms];
  atom.flag = new int*    [natoms];
  atom.x    = new double* [natoms];
  atom.v    = new double* [natoms];
  for (unsigned int n=0; n<natoms; n++) {
    atom.flag[n] = new int    [ndims];
    atom.x[n]    = new double [ndims];
    atom.v[n]    = new double [ndims];
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
    finit >> iatom >> itype >> xcoor >> ycoor >> zcoor >> xflag >> yflag >> zflag;
    //cout << setw(18) << n << setw(18) << iatom << endl;
    getline(finit,skipl);
    atom.type[iatom-1] = itype;
    atom.flag[iatom-1][0] = xflag;
    atom.flag[iatom-1][1] = yflag;
    atom.flag[iatom-1][2] = zflag;
    atom.x[iatom-1][0] = xcoor;
    atom.x[iatom-1][1] = ycoor;
    atom.x[iatom-1][2] = zcoor;
  }

  // list of atom velocities
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);
  for (unsigned int n=0; n<natoms; n++) {
    finit >> iatom >> xvelo >> yvelo >> zvelo;
    //cout << setw(18) << n << setw(18) << iatom << endl;
    getline(finit,skipl);
    atom.v[iatom-1][0] = xvelo;
    atom.v[iatom-1][1] = yvelo;
    atom.v[iatom-1][2] = zvelo;
  }

  // close the file if it is the end
  getline(finit,skipl);
  if(!finit.eof()) cout << "Not the end of file??????"  << endl;
  finit.close();

  // ===========================================================================
  // Color the solid grain atoms according to some certain rules
  // ===========================================================================

  // a simple example as to color a circular solid grain

  for (unsigned int n=0; n<natoms; n++) {

    double radii = sqrt(atom.x[n][0]*atom.x[n][0]+
                        atom.x[n][1]*atom.x[n][1]+
                        atom.x[n][2]*atom.x[n][2]);

    if(radii > 10.0) {
      atom.type[n] ++;
      nfluid ++;
    }
  }
  nsolid = natoms - nfluid;

  matoms = natoms;
  mtypes = 2;

  // ===========================================================================
  // Write to the new LAMMPS data file
  // ===========================================================================

  fpost.open("restart.dat", ios::out);

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
          << atom.type[n] << " "
          << atom.x[n][0] << " "
          << atom.x[n][1] << " "
          << atom.x[n][2] << " "
          << atom.flag[n][0] << " "
          << atom.flag[n][1] << " "
          << atom.flag[n][2] << endl;
  }

  fpost.close();


  // release dynamic arrays
  delete [] atom.type;
  delete [] atom.flag;
  delete [] atom.x;
  delete [] atom.v;

  cout << endl;
  cout << "Done processing the LAMMPS restart file..." << endl;
  cout << "Exit program." << endl;
  cout << endl;
}
