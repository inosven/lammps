#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
  // ===========================================================================
  // A header message
  // ===========================================================================

  cout << endl;
  cout << "################################################################\n";
  cout << endl;
  cout << endl;
  cout << "  Calculate Solution Distribution along a Sourced line in 3-D\n";
  cout << endl;
  cout << "                         Yidong Xia\n";
  cout << endl;
  cout << "################################################################\n";
  cout << endl;

  // ===========================================================================
  // define data types
  // ===========================================================================

  // integers
  unsigned const int ndeli = 16;
  unsigned int timestep,nchunks,tcount,chunk,nvars;

  // floats
  double ncount,var1;
  double x,y,z;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xmin,xmax,ymin,ymax,zmin,zmax;

  // characters
  char cinit[80];

  // strings
  string skipl;

  // files
  ifstream finit;
  ofstream fpost;

  fpost.open("sprofile.dat", ios::out);
  //fpost.precision(16);
  //fpost << scientific;

  // ===========================================================================
  // Open LAMMPS profile data
  // ===========================================================================

  // type in file name
  do {
    cout << "Enter the name of LAMMPS profile data:\n";
    cin >> cinit;
    finit.open(cinit, ios::in);
  } while (!finit.is_open());
  cout << "The LAMMPS profile data is " << cinit << endl << endl;

  cout << "Enter how many variables to operate:\n";
  cin >> nvars;
  cout << "nvars = " << nvars << endl << endl;

  cout << "Define the screen box range: xlo xhi ylo yhi zlo zhi\n";
  cin >> xlo >> xhi >> ylo >> yhi >> zlo >> zhi;
  if (xlo >= xhi) {
    cout << "Error: xlo >= xhi\n" << endl;
    exit(0);
  }
  if (ylo >= yhi) {
    cout << "Error: ylo >= yhi\n" << endl;
    exit(0);
  }
  if (zlo >= zhi) {
    cout << "Error: zlo >= zhi\n" << endl;
    exit(0);
  }
  cout << "xlo  = " << xlo << ", xhi = " << xhi << endl;
  cout << "ylo  = " << ylo << ", yhi = " << yhi << endl;
  cout << "zlo  = " << zlo << ", zhi = " << zhi << endl << endl;

  // Set header of statistical information
  cout << left << setw(ndeli) << "Timestep"
       << left << setw(ndeli) << "Nchunks"
       << left << setw(ndeli) << "Total count"
       << endl;
  cout << left << setw(ndeli) << "=============="
       << left << setw(ndeli) << "=============="
       << left << setw(ndeli) << "=============="
       << endl << endl;

  // skip the first three lines
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);

  while (!finit.eof()) {
    // header information
    finit >> timestep >> nchunks >> tcount; getline(finit,skipl);

    xmin = 0.0; xmax = 0.0;
    ymin = 0.0; ymax = 0.0;
    zmin = 0.0; zmax = 0.0;

    // loop over chunks
    for (unsigned int n=0; n<nchunks; n++) {
      // read the line
      finit >> chunk >> x >> y >> z >> ncount >> var1; getline(finit,skipl);

      if (x < xlo || x > xhi) continue;
      if (y < ylo || y > yhi) continue;
      if (z < zlo || z > zhi) continue;

      // write to file
      fpost << setw(ndeli) << x
            << setw(ndeli) << y
            << setw(ndeli) << z
            << setw(ndeli) << var1
            << endl;
    }
  }

  fpost.close();
  finit.close();

  cout << endl;
  cout << "Done processing the LAMMPS profile data..." << endl;
  cout << "Exit program." << endl;
  cout << endl;
}
