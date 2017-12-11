#include <fstream>
#include <iomanip>
#include <iostream>
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
  cout << "            A Simple Code to Do Quantitative Analysis\n";
  cout << "                for LAMMPS Default 'dump' Data\n";
  cout << endl;
  cout << "              Author: Yidong Xia (yidong.xia@inl.gov)\n";
  cout << "                     Idaho National Laboratory\n";
  cout << endl;
  cout << endl;
  cout << " Limit: the particle search path can only be parallel to axes\n";
  cout << endl;
  cout << " GNUPLOT is suggested to draw the x-y plot\n";
  cout << endl;
  cout << "################################################################\n";
  cout << endl;

  // ===========================================================================
  // I/O operation of input script (the file name is fixed for convenience)
  // ===========================================================================

  // 1). Open
  // ========

  char cinput[] = "xyinp";
  ifstream finput;
  finput.open(cinput, ios::in);
  if (!finput.is_open())
  {
    cout << "No 'xyinp' script file is found in the current directory!\n";
    exit(0);
  }

  // 2). Read
  // ========

  // Source box range
  double xmin = 0, xmax = 0;
  double ymin = 0, ymax = 0;
  double zmin = 0, zmax = 0;

  string skipl;
  getline(finput,skipl);
  finput >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
  cout << "Source box range (normalized):\n";
  cout << setw(12) << "xmin = " << setw(16) << xmin
       << setw(12) << "xmax = " << setw(16) << xmax << endl;
  cout << setw(12) << "ymin = " << setw(16) << ymin
       << setw(12) << "ymax = " << setw(16) << ymax << endl;
  cout << setw(12) << "zmin = " << setw(16) << zmin
       << setw(12) << "zmax = " << setw(16) << zmax << endl;
  cout << endl;
  getline(finput,skipl);
  if (xmin >= xmax)
  {
    cout << "Error in 'xyinp': xmin >= xmax\n";
    cout << endl;
    exit(0);
  }
  if (ymin >= ymax)
  {
    cout << "Error in 'xyinp': ymin >= ymax\n";
    cout << endl;
    exit(0);
  }
  if (zmin >= zmax)
  {
    cout << "Error in 'xyinp': zmin >= zmax\n";
    cout << endl;
    exit(0);
  }

  int idire = -1;
  getline(finput,skipl);
  finput >> idire;
  getline(finput,skipl);
  if (idire == 0)
    cout << "Source line direction: x\n";
  else if (idire == 1)
    cout << "Source line direction: y\n";
  else if (idire == 2)
    cout << "Source line direction: z\n";
  else
  {
    cout << "Error in 'xyinp': idire must be 0 or 1 or 2\n";
    cout << "idire = " << idire << endl;
    cout << endl;
    exit(0);
  }

  cout << endl;

  // 3). Close
  // =========

  finput.close();

  // ===========================================================================
  // Open LAMMPS dump data
  // ===========================================================================

  cout << "Enter the name of LAMMPS dump file:\n";
  char cdump[80];
  cin >> cdump;
  ifstream fdump;
  fdump.open (cdump, ios::in);
  while (!fdump.is_open())
  {
    cout << "No file named " << cdump << " is found!\n";
    cin >> cdump;
    fdump.open(cdump, ios::in);
  }
  cout << "The LAMMPS dump file is "  << cdump << endl << endl;

  // ===========================================================================
  // Read and process LAMMPS dump data
  // ===========================================================================

  // No. of timestep
  int istep = -1;
  int jstep = -1;
  // No. of particles
  int natom = -1;

  // Simulation box range
  double xlo = 0, xhi = 0;
  double ylo = 0, yhi = 0;
  double zlo = 0, zhi = 0;

  int ndeli0 =  5;
  int ndeli1 =  7;
  int ndeli2 = 10;
  int ndeli3 = 10;
  int ndeli4 = 18;
  int ndeli5 = 14;

  // Set header of statistical information
  cout << setw(ndeli1) << "Frame"
       << setw(ndeli2) << "Timestep"
       << setw(ndeli3) << "Atoms"
       << setw(ndeli3) << "Inbox"
       << setw(ndeli4) << "Xlo"
       << setw(ndeli4) << "Xhi"
       << setw(ndeli4) << "Ylo"
       << setw(ndeli4) << "Yhi"
       << setw(ndeli4) << "Zlo"
       << setw(ndeli4) << "Zhi"
       << endl;
  cout << setw(ndeli1) << "======"
       << setw(ndeli2) << "========"
       << setw(ndeli3) << "========"
       << setw(ndeli3) << "========"
       << setw(ndeli4) << "================"
       << setw(ndeli4) << "================"
       << setw(ndeli4) << "================"
       << setw(ndeli4) << "================"
       << setw(ndeli4) << "================"
       << setw(ndeli4) << "================"
       << endl;

  // Loop over timesteps
  int ifram = 0;
  while (!fdump.eof())
  {
    ifram += 1;
    if (ifram > 1) getline(fdump,skipl);
    getline(fdump,skipl);
    fdump >> istep;

    if (jstep == istep) break;

    getline(fdump,skipl);
    getline(fdump,skipl);
    fdump >> natom;

    getline(fdump,skipl);
    getline(fdump,skipl);
    fdump >> xlo >> xhi;
    fdump >> ylo >> yhi;
    fdump >> zlo >> zhi;

    getline(fdump,skipl);
    getline(fdump,skipl);

    // Initialize post data file
    string fname;
    if (ifram < 10)
      fname = "xypos.0000" + to_string(ifram) + ".txt";
    else if (ifram < 100)
      fname = "xypos.000"  + to_string(ifram) + ".txt";
    else if (ifram < 1000)
      fname = "xypos.00"   + to_string(ifram) + ".txt";
    else if (ifram < 10000)
      fname = "xypos.0"    + to_string(ifram) + ".txt";
    else if (ifram < 100000)
      fname = "xypos."     + to_string(ifram) + ".txt";
    else
    {
      cout << endl;
      cout << "Error: running out of post file sequence limit: 99999\n";
      cout << endl;
      exit(0);
    }
    ofstream fpost;
    fpost.open(fname, ios::out);
    fpost << setw(ndeli5) << "coord"
          << setw(ndeli5) << "velox" << setw(ndeli5) << "veloy" << setw(ndeli5) << "veloz"
          << endl;
    fpost.precision(5);
    fpost << scientific;


    int iatom = 0;
    int jatom = 0;
    for (int ii = 0; ii < natom; ii += 1)
    {
      iatom += 1;

      //cout << "iatom = " << iatom << endl;

      int ipoin = -1;
      int itype = -1;
      double coord = 0;
      double xcoor = 0, ycoor = 0, zcoor = 0;
      double velox = 0, veloy = 0, veloz = 0;

      fdump >> ipoin
            >> itype
            >> xcoor >> ycoor >> zcoor
            >> velox >> veloy >> veloz;

      if (idire == 0)
        coord = xcoor;
      else if (idire == 1)
        coord = ycoor;
      else if (idire == 2)
        coord = zcoor;

      // Check if this atom is within the search box
      if ( xcoor >= xmin && xcoor <= xmax &&
           ycoor >= ymin && ycoor <= ymax &&
           zcoor >= zmin && zcoor <= zmax )
      {
        jatom += 1;

        fpost << setw(ndeli5) << coord
              << setw(ndeli5) << velox << setw(ndeli5) << veloy << setw(ndeli5) << veloz
              << endl;
      }
    }

    fpost.close();

    if (iatom != natom)
    {
      cout << "Error: iatom NOT equal to natom at timestep = " << istep << endl;
      cout << "iatom = " << iatom << endl;
      cout << "natom = " << natom << endl;
      cout << endl;
      exit(0);
    }

    cout
         << setw(ndeli1) << ifram
         << setw(ndeli2) << istep
         << setw(ndeli3) << iatom
         << setw(ndeli3) << jatom
         << setw(ndeli4) << xlo
         << setw(ndeli4) << xhi
         << setw(ndeli4) << ylo
         << setw(ndeli4) << yhi
         << setw(ndeli4) << zlo
         << setw(ndeli4) << zhi
         << endl;

    if (jatom == 0)
    {
      cout << endl;
      cout << "Error: NO atom is found! Please adjust the search box range...\n";
      cout << endl;
      exit(0);
    }

    jstep = istep;
  }
  cout << endl;
  cout << "Done reading the LAMMPS dump file..." << endl;
  cout << "Exit program." << endl;
  cout << endl;
}
