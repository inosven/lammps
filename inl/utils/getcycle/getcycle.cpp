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
  cout << "   A Post-Process Code to Calculate a Droplet Dumping Histgram\n";
  cout << endl;
  cout << "            Yidong Xia (Idaho National Laboratory)\n";
  cout << endl;
  cout << "################################################################\n";
  cout << endl;

  // ===========================================================================
  // define data types
  // ===========================================================================

  // integers
  unsigned const int ndeli = 16;
  unsigned int timestep,nchunks,tcount,chunk,ncount;

  // floats
  double coord,coordmin,coordmax,density,width;

  // characters
  char cinit[80];

  // strings
  string skipl;

  // files
  ifstream finit;
  ofstream fpost;

  fpost.open("cycle.dat", ios::out);
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

  // Set header of statistical information
  cout << left << setw(ndeli) << "Timestep"
       << left << setw(ndeli) << "Nchunks"
       << left << setw(ndeli) << "Total count"
       << left << setw(ndeli) << "Width"
       << endl;
  cout << left << setw(ndeli) << "=============="
       << left << setw(ndeli) << "=============="
       << left << setw(ndeli) << "=============="
       << left << setw(ndeli) << "=============="
       << endl << endl;

  // skip the first three lines
  getline(finit,skipl);
  getline(finit,skipl);
  getline(finit,skipl);
  //getline(finit,skipl);

  while (!finit.eof()) {
    // header information
    finit >> timestep >> nchunks >> tcount;getline(finit,skipl);

    // note that there should not be single particles escaping the droplet
    coordmin = 0.0;
    coordmax = 0.0;

    // loop over chunks
    for (unsigned int n=0; n<nchunks; n++) {
      // read the line
      finit >> chunk >> coord >> ncount >> density; getline(finit,skipl);

      if (density > 1e-3) {
        if (coord < coordmin) coordmin = coord;
        if (coord > coordmax) coordmax = coord;
      }

      width = 0.5*(coordmax - coordmin);
    }

    // print the results at this timestep
    //cout << left << setw(ndeli) << timestep
    //     << left << setw(ndeli) << nchunks
    //     << left << setw(ndeli) << tcount
    //     << left << setw(ndeli) << width
    //     << endl;

    // write to file
    fpost << setw(ndeli) << timestep << setw(ndeli) << width << endl;
  }

  fpost.close();
  finit.close();

  cout << endl;
  cout << "Done processing the LAMMPS restart file..." << endl;
  cout << "Exit program." << endl;
  cout << endl;
}
