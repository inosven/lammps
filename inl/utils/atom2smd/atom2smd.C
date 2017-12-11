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
  int    *  tag;
  int    *  type;
  int    ** flag;
  double ** x;
};

int main(int argc, char **argv)
{
  // A header message
  cout << endl;
  cout << "################################################################\n";
  cout << endl;
  cout << endl;
  cout << "   A code to convert the standard LAMMPS data to the SMD data\n";
  cout << endl;
  cout << "                          Yidong Xia \n";
  cout << endl;
  cout << "################################################################\n";
  cout << endl;

  // integers
  unsigned const int ndeli  = 18;
  unsigned int natoms = 0;
  unsigned int ntypes = 0;
  unsigned int tag = 0;
  unsigned int type = 0;
  unsigned int mol = 0;

  // floats
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double vfrac,rmass,radius,contact_radius,x,y,z;

  // characters
  char cinit[80];

  // strings
  string skipl;

  // files
  ifstream finpu;
  ifstream finit;
  ofstream fpost;

  // open configuration input file
  finpu.open("atom2smd.inp", ios::in);
  if (!finpu.is_open()) {
    cout << "Error: no 'atom2smd.inp' is found in the current directory!" << endl;
    exit(0);
  }
  getline(finpu,skipl);
  finpu >> mol >> vfrac >> rmass >> radius >> contact_radius; getline(finpu,skipl);
  finpu.close();

  // open an LAMMPS data file written by OVITO
  do {
    cout << "Enter the name of an LAMMPS data file written by OVITO:\n";
    cin >> cinit;
    finit.open(cinit, ios::in);
  } while (!finit.is_open());
  cout << "The LAMMPS data file written by OVITO is "  << cinit << endl << endl;

  getline(finit,skipl);
  finit >> natoms; getline(finit,skipl);
  finit >> ntypes; getline(finit,skipl);
  finit >> xlo >> xhi; getline(finit,skipl);
  finit >> ylo >> yhi; getline(finit,skipl);
  finit >> zlo >> zhi; getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);


  // open the new LAMMPS-SMD data file
  fpost.open("lammps_smd_restart.dat", ios::out);

  fpost.precision(6);
  fpost << scientific;

  // first line
  fpost << "LAMMPS-SMD data file" << endl;

  // head information section
  fpost << natoms << " atoms" << endl;
  fpost << ntypes << " atom types" << endl;

  // range box section
  fpost << xlo <<  " " << xhi << " xlo xhi" << endl;
  fpost << ylo <<  " " << yhi << " ylo yhi" << endl;
  fpost << zlo <<  " " << zhi << " zlo zhi" << endl;

  // Atoms section
  fpost << endl;
  fpost << "Atoms" << endl;
  fpost << endl;

  fpost.precision(8);

  for (unsigned int n=1; n<=natoms; n++) {
    finit >> tag >> type >> x >> y >> z; getline(finit,skipl);

    //fpost << right << setw(9) << tag
    fpost << right << setw(9) << n
          << right << setw(3) << type
          << right << setw(3) << mol
          << right << setw(16) << vfrac
          << right << setw(16) << rmass
          << right << setw(16) << radius
          << right << setw(16) << contact_radius
          << right << setw(16) << x
          << right << setw(16) << y
          << right << setw(16) << z
          << right << setw(16) << x
          << right << setw(16) << y
          << right << setw(16) << z
          << endl;
  }
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);

  // close the files
  finit.close();
  cout << "Closed the standard LAMMPS data file\n" << endl;
  fpost.close();
  cout << "Closed the LAMMPS-SMD data file\n" << endl;

  exit(0);
}
