#include <fstream>
#include <iomanip>
#include <random>
#include <iostream>
#include <time.h>
#include <vector>

//struct _t {
//  int    *  type;
//  int    ** flag;
//  bool   *  wall;
//  double ** x;
//};


void print_message()
{
  std::cout << "\n";
  std::cout << "#############################################################\n";
  std::cout << "\n";
  std::cout << "  Code Name: LATTICE2DPD v2\n";
  std::cout << "  Author(s): Yidong Xia (Idaho National Lab)\n";
  std::cout << "  Functions: Generate initial DPD bead locations from voxels\n";
  std::cout << "\n";
  std::cout << "#############################################################\n";
  std::cout << "\n";
}

int main(int argc, char **argv)
{
  print_message();

  clock_t timeStart = clock();

  // ========================
  // Parameter initialization
  // ========================

  bool outSolid = false;
  bool outFluid = false;

  unsigned int flowDirection = 0;

  unsigned int typeVoid = 0;
  unsigned int typeWall = 1;
  unsigned int typePore = 2;

  unsigned long ix,iy,iz;

  unsigned int cropLo = 0, cropHi = 0;

  const unsigned long nMaxTypes = 12;
  unsigned long nTypes = 0;
  unsigned long nSolidsPerLattice = 0;
  unsigned long nFluidsPerLattice = 0;

  unsigned long numSolidBeads = 0, numFluidBeads = 0;

  unsigned long numSolidVox = 0, numFluidVox = 0;

  unsigned long nxSolidVoxLo = 0, nxSolidVoxHi = 0;
  unsigned long nySolidVoxLo = 0, nySolidVoxHi = 0;
  unsigned long nzSolidVoxLo = 0, nzSolidVoxHi = 0;
  unsigned long nxFluidVoxLo = 0, nxFluidVoxHi = 0;
  unsigned long nyFluidVoxLo = 0, nyFluidVoxHi = 0;
  unsigned long nzFluidVoxLo = 0, nzFluidVoxHi = 0;

  unsigned long hx  = 0, hy  = 0, hz  = 0;

  double xrnd, yrnd, zrnd;

  // strings

  std::string skipLine;
  std::string infoFileName = "info";
  std::string ctrlFileName = "control";
  std::string inpSolidFileName = "inp_solid.dat";
  std::string inpFluidFileName = "inp_fluid.dat";
  std::string outSolidFileName = "out_solid.dat";
  std::string outFluidFileName = "out_fluid.dat";

  // I/O files

  std::ifstream infoFile;
  std::ifstream ctrlFile;
  std::ifstream inpSolidFile;
  std::ifstream inpFluidFile;
  std::ofstream outSolidFile;
  std::ofstream outFluidFile;

  // ==========================
  // Open and read control file
  // ==========================

  ctrlFile.open(ctrlFileName, std::ios::in);

  if (!ctrlFile.is_open())
  {
    std::cout << "\n" << "Fatal: file " << ctrlFileName << " does not exist.\n";
    std::exit(0);
  }

  // Read the flow direction indicator

  std::getline(ctrlFile,skipLine);
  ctrlFile >> flowDirection; std::getline(ctrlFile,skipLine);

  if (flowDirection > 3)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " flowDirection = " << flowDirection << " is not a valid value.\n";
    std::exit(0);
  }

  // Read the lower and higher bounds for region cropping

  std::getline(ctrlFile,skipLine);
  ctrlFile >> cropLo >> cropHi; std::getline(ctrlFile,skipLine);

  if (cropLo >= cropHi || cropHi > 100)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName
      << " cropLo = " << cropLo << ", "
      << " cropHi = " << cropHi << "\n";
    std::exit(0);
  }


  // Read the scaling coefficient in each direction

  std::getline(ctrlFile,skipLine);
  ctrlFile >> hx >> hy >> hz; std::getline(ctrlFile,skipLine);

  // Read the solid and fluid number densities

  std::getline(ctrlFile,skipLine);
  ctrlFile >> nSolidsPerLattice >> nFluidsPerLattice; std::getline(ctrlFile,skipLine);

  // Read the number of types of beads

  std::getline(ctrlFile,skipLine);
  ctrlFile >> nTypes; std::getline(ctrlFile,skipLine);

  if (nTypes > nMaxTypes)
  {
    std::cout
      << "\n" << "Fatal: in file " << ctrlFileName << "\n"
      << "  nTypes = " << nTypes << "\n"
      << "  is larger than\n"
      << "  nMaxTypes = " << nMaxTypes << "\n";
    std::exit(0);
  }

  // Read whether to output solid and fluid

  std::getline(ctrlFile,skipLine);
  ctrlFile >> outSolid >> outFluid; std::getline(ctrlFile,skipLine);

  ctrlFile.close();

  // =======================
  // Open and read info file
  // =======================

  // See if the file exists

  infoFile.open(infoFileName, std::ios::in);

  if (!infoFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << infoFileName << " does not exist.\n";
    std::exit(0);
  }

  // Get the number of solid voxels, and low & high ends of indices

  std::getline(infoFile,skipLine);
  infoFile
    >> numSolidVox;
  std::getline(infoFile,skipLine);

  std::getline(infoFile,skipLine);
  infoFile
    >> nxSolidVoxLo >> nxSolidVoxHi
    >> nySolidVoxLo >> nySolidVoxHi
    >> nzSolidVoxLo >> nzSolidVoxHi;
  std::getline(infoFile,skipLine);

  // Get the number of fluid voxels, and low & high ends of indices

  std::getline(infoFile,skipLine);
  infoFile
    >> numFluidVox;
  std::getline(infoFile,skipLine);

  std::getline(infoFile,skipLine);
  infoFile
    >> nxFluidVoxLo >> nxFluidVoxHi
    >> nyFluidVoxLo >> nyFluidVoxHi
    >> nzFluidVoxLo >> nzFluidVoxHi;
  std::getline(infoFile,skipLine);

  infoFile.close();

  // ========================================================================
  // A full 3D matrix for storing the region info is created to help identify
  // inner sides of wall voxels adjacent to pore space
  // ... void = 0 (default)
  // ... wall = 1
  // ... pore = 2
  // ========================================================================

  unsigned long nxVox = nxSolidVoxHi - nxSolidVoxLo + 1;
  unsigned long nyVox = nySolidVoxHi - nySolidVoxLo + 1;
  unsigned long nzVox = nzSolidVoxHi - nzSolidVoxLo + 1;

  unsigned long nxLattice = nxVox * hx;
  unsigned long nyLattice = nyVox * hy;
  unsigned long nzLattice = nzVox * hz;

  std::vector<std::vector<std::vector<unsigned int> > > LATTICE;
  LATTICE.resize(nxLattice+2);
  for (unsigned long i = 0; i < (nxLattice+2); i++)
  {
    LATTICE[i].resize(nyLattice+2);
    for (unsigned long j = 0; j < (nyLattice+2); j++)
      LATTICE[i][j].resize(nzLattice+2);
  }

  // ===========================
  // Display on screen and check
  // ===========================

  unsigned long xlo = 0;
  unsigned long ylo = 0;
  unsigned long zlo = 0;
  unsigned long xhi = nxLattice;
  unsigned long yhi = nyLattice;
  unsigned long zhi = nzLattice;

  switch ( flowDirection )
  {
    case 0:
      break;

    case 1:

      xlo = nxVox * cropLo / 100; xlo = xlo * hx;
      xhi = nxVox * cropHi / 100; xhi = xhi * hx;
      break;

    case 2:

      ylo = nyVox * cropLo / 100; ylo = ylo * hy;
      yhi = nyVox * cropHi / 100; yhi = yhi * hy;
      break;

    case 3:

      zlo = nzVox * cropLo / 100; zlo = zlo * hz;
      zhi = nzVox * cropHi / 100; zhi = zhi * hz;
      break;
  }

  unsigned long numSolidLattice = numSolidVox * hx * hy * hz;
  unsigned long numFluidLattice = numFluidVox * hx * hy * hz;

  std::cout
    << std::left << std::setfill('.')
    << std::setw(40) << "numSolidVox  "       << "  " << numSolidVox  << "\n"
    << std::setw(40) << "numFluidVox  "       << "  " << numFluidVox  << "\n"
    << std::setw(40) << "nxSolidVoxLo  "      << "  " << nxSolidVoxLo << "\n"
    << std::setw(40) << "nxSolidVoxHi  "      << "  " << nxSolidVoxHi << "\n"
    << std::setw(40) << "nySolidVoxLo  "      << "  " << nySolidVoxLo << "\n"
    << std::setw(40) << "nySolidVoxHi  "      << "  " << nySolidVoxHi << "\n"
    << std::setw(40) << "nzSolidVoxLo  "      << "  " << nzSolidVoxLo << "\n"
    << std::setw(40) << "nzSolidVoxHi  "      << "  " << nzSolidVoxHi << "\n"
    << std::setw(40) << "nxFluidVoxLo  "      << "  " << nxFluidVoxLo << "\n"
    << std::setw(40) << "nxFluidVoxHi  "      << "  " << nxFluidVoxHi << "\n"
    << std::setw(40) << "nyFluidVoxLo  "      << "  " << nyFluidVoxLo << "\n"
    << std::setw(40) << "nyFluidVoxHi  "      << "  " << nyFluidVoxHi << "\n"
    << std::setw(40) << "nzFluidVoxLo  "      << "  " << nzFluidVoxLo << "\n"
    << std::setw(40) << "nzFluidVoxHi  "      << "  " << nzFluidVoxHi << "\n"
    << std::setw(40) << "nxVox  "             << "  " << nxVox        << "\n"
    << std::setw(40) << "nyVox  "             << "  " << nyVox        << "\n"
    << std::setw(40) << "nzVox  "             << "  " << nzVox        << "\n"
    << std::setw(40) << "flowDirection  "     << "  " << flowDirection<< "\n"
    << std::setw(40) << "hx  "                << "  " << hx           << "\n"
    << std::setw(40) << "hy  "                << "  " << hy           << "\n"
    << std::setw(40) << "hz  "                << "  " << hz           << "\n"
    << std::setw(40) << "xlo  "               << "  " << xlo          << "\n"
    << std::setw(40) << "xhi  "               << "  " << xhi          << "\n"
    << std::setw(40) << "ylo  "               << "  " << ylo          << "\n"
    << std::setw(40) << "yhi  "               << "  " << yhi          << "\n"
    << std::setw(40) << "zlo  "               << "  " << zlo          << "\n"
    << std::setw(40) << "zhi  "               << "  " << zhi          << "\n"
    << std::setw(40) << "nTypes  "            << "  " << nTypes       << "\n"
    << std::setw(40) << "nSolidsPerLattice  " << "  " << nSolidsPerLattice<< "\n"
    << std::setw(40) << "nFluidsPerLattice  " << "  " << nFluidsPerLattice<< "\n"
    << std::setw(40) << "outSolid  "          << "  " << (outSolid ? "Yes" : "No") << "\n"
    << std::setw(40) << "outFluid  "          << "  " << (outFluid ? "Yes" : "No") << "\n"
    << "\n";

  // =====================================
  // Open and read input solid voxel files
  // =====================================

  inpSolidFile.open(inpSolidFileName, std::ios::in);

  if (!inpSolidFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << inpSolidFileName << " does not exist.\n";
    std::exit(0);
  }

  for (unsigned long i = 0; i < numSolidVox; i++)
  {
    inpSolidFile >> ix >> iy >> iz;
    std::getline(inpSolidFile,skipLine);

    bool isOutOfRange = false;

    if (ix < nxSolidVoxLo || ix > nxSolidVoxHi) isOutOfRange = true;
    if (iy < nySolidVoxLo || iy > nySolidVoxHi) isOutOfRange = true;
    if (iz < nzSolidVoxLo || iz > nzSolidVoxHi) isOutOfRange = true;
    if (isOutOfRange)
    {
      std::cout << "\n\n"
        << "Error: in file " << inpSolidFileName << ", Line " << (i+1) << ", "
        << "(ix,iy,iz) = (" << ix << "," << iy << "," << iz << ") out of range!\n";
      std::exit(0);
    }

    // Remember each index starts from ONE
    unsigned long I = ix - nxSolidVoxLo + 1;
    unsigned long J = iy - nySolidVoxLo + 1;
    unsigned long K = iz - nzSolidVoxLo + 1;

    unsigned long I0 = (I-1) * hx + 1;
    unsigned long I1 =     I * hx;
    unsigned long J0 = (J-1) * hy + 1;
    unsigned long J1 =     J * hy;
    unsigned long K0 = (K-1) * hz + 1;
    unsigned long K1 =     K * hz;

    for (unsigned long ii = I0; ii <= I1; ii++)
      for (unsigned long jj = J0; jj <= J1; jj++)
        for (unsigned long kk = K0; kk <= K1; kk++)
          LATTICE[ii][jj][kk] = typeWall;
  }

  inpSolidFile.close();

  // =====================================
  // Open and read input fluid voxel files
  // =====================================

  inpFluidFile.open(inpFluidFileName, std::ios::in);

  if (!inpFluidFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << inpFluidFileName << " does not exist.\n";
    std::exit(0);
  }

  for (unsigned long i = 0; i < numFluidVox; i++)
  {
    inpFluidFile >> ix >> iy >> iz;
    std::getline(inpFluidFile,skipLine);

    bool isOutOfRange = false;

    if (ix < nxFluidVoxLo || ix > nxFluidVoxHi) isOutOfRange = true;
    if (iy < nyFluidVoxLo || iy > nyFluidVoxHi) isOutOfRange = true;
    if (iz < nzFluidVoxLo || iz > nzFluidVoxHi) isOutOfRange = true;
    if (isOutOfRange)
    {
      std::cout << "\n\n"
        << "Error: in file " << inpFluidFileName << ", Line " << (i+1) << ", "
        << "(ix,iy,iz) = (" << ix << "," << iy << "," << iz << ") out of range!\n";
      std::exit(0);
    }

    // Remember each index starts from ONE
    unsigned long I = ix - nxSolidVoxLo + 1;
    unsigned long J = iy - nySolidVoxLo + 1;
    unsigned long K = iz - nzSolidVoxLo + 1;

    unsigned long I0 = (I-1) * hx + 1;
    unsigned long I1 =     I * hx;
    unsigned long J0 = (J-1) * hy + 1;
    unsigned long J1 =     J * hy;
    unsigned long K0 = (K-1) * hz + 1;
    unsigned long K1 =     K * hz;

    for (unsigned long ii = I0; ii <= I1; ii++)
      for (unsigned long jj = J0; jj <= J1; jj++)
        for (unsigned long kk = K0; kk <= K1; kk++)
          LATTICE[ii][jj][kk] = typePore;
  }

  inpFluidFile.close();

  // ==========================================================
  // At the two region boundaries normal to the flow direction,
  // replace previously added wall voxels with void voxels,
  // and replace original void voxels with wall voxels.
  // So fluid flow can be contained only in pore network.
  // ==========================================================

  switch ( flowDirection )
  {
    // In all 6 boundary planes: convert wall voxels to void voxels
    case 0:

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long j = 1; j <= nyLattice; j++)
          for (unsigned long k = 1; k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      numSolidBeads += nSolidsPerLattice;
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if (LATTICE[i][j][k] == typePore)
              numFluidBeads += nFluidsPerLattice;
          }

      break;

    // In xlo and xhi planes:
    //   * switch over between wall voxels <--> void voxels;
    //   * count solid particle numbers in the region boundary voxels.
    case 1:

      for (unsigned long j = 1; j <= nyLattice; j++)
        for (unsigned long k = 1; k <= nzLattice; k++)
        {
          for (unsigned long i = 1; i <= hx; i++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              LATTICE[i][j][k] = typeVoid;
              numSolidLattice -= 1;
            }
            else
            {
              LATTICE[i][j][k] = typeWall;
              numSolidLattice += 1;
              if (i == hx) numSolidBeads += nSolidsPerLattice;
            }
          }

          for (unsigned long i = (nxLattice-hx+1); i <= nxLattice; i++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              LATTICE[i][j][k] = typeVoid;
              numSolidLattice -= 1;
            }
            else
            {
              LATTICE[i][j][k] = typeWall;
              numSolidLattice += 1;
              if (i == (nxLattice-hx+1)) numSolidBeads += nSolidsPerLattice;
            }
          }
        }

      for (unsigned long i = (hx+1); i <= (nxLattice-hx); i++)
        for (unsigned long j = 1; j <= nyLattice; j++)
          for (unsigned long k = 1; k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      numSolidBeads += nSolidsPerLattice;
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if (LATTICE[i][j][k] == typePore)
              numFluidBeads += nFluidsPerLattice;
          }

      break;

    // In ylo and yhi planes:
    //   * switch over between wall voxels <--> void voxels;
    //   * count solid particle numbers in the region boundary voxels.
    case 2:

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long k = 1; k <= nzLattice; k++)
        {
          for (unsigned long j = 1; j <= hy; j++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              LATTICE[i][j][k] = typeVoid;
              numSolidLattice -= 1;
            }
            else
            {
              LATTICE[i][j][k] = typeWall;
              numSolidLattice += 1;
              if (j == hy) numSolidBeads += nSolidsPerLattice;
            }
          }

          for (unsigned long j = (nyLattice-hy+1); j <= nyLattice; j++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              LATTICE[i][j][k] = typeVoid;
              numSolidLattice -= 1;
            }
            else
            {
              LATTICE[i][j][k] = typeWall;
              numSolidLattice += 1;
              if (j == (nyLattice-hy+1)) numSolidBeads += nSolidsPerLattice;
            }
          }
        }

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long j = (hy+1); j <= (nyLattice-hy); j++)
          for (unsigned long k = 1; k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      numSolidBeads += nSolidsPerLattice;
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if (LATTICE[i][j][k] == typePore)
              numFluidBeads += nFluidsPerLattice;
          }

      break;

    // In zlo and zhi planes:
    //   * switch over between wall voxels <--> void voxels;
    //   * count solid particle numbers in the region boundary voxels.
    case 3:

/*
      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long j = 1; j <= nyLattice; j++)
        {
          for (unsigned long k = 1; k <= hz; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              LATTICE[i][j][k] = typeVoid;
              numSolidLattice -= 1;
            }
            else
            {
              LATTICE[i][j][k] = typeWall;
              numSolidLattice += 1;
              if (k == hz) numSolidBeads += nSolidsPerLattice;
            }
          }

          for (unsigned long k = (nzLattice-hz+1); k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              LATTICE[i][j][k] = typeVoid;
              numSolidLattice -= 1;
            }
            else
            {
              LATTICE[i][j][k] = typeWall;
              numSolidLattice += 1;
              if (k == (nzLattice-hz+1)) numSolidBeads += nSolidsPerLattice;
            }
          }
        }

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long j = 1; j <= nyLattice; j++)
          for (unsigned long k = (hz+1); k <= (nzLattice-hz); k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      numSolidBeads += nSolidsPerLattice;
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if (LATTICE[i][j][k] == typePore)
              numFluidBeads += nFluidsPerLattice;
          }
*/

      for (unsigned long i = 1; i <= xhi; i++)
        for (unsigned long j = 1; j <= yhi; j++)
        {
          unsigned long k;

          k = zlo+1;

          if (LATTICE[i][j][k] == typeWall)
            numSolidBeads += nSolidsPerLattice;
          else if (LATTICE[i][j][k] == typePore)
            numFluidBeads += nFluidsPerLattice;
          else
          {
            LATTICE[i][j][k] = typeWall;
            numSolidBeads += nSolidsPerLattice;
          }

          k = zhi;

          if (LATTICE[i][j][k] == typeWall)
            numSolidBeads += nSolidsPerLattice;
          else if (LATTICE[i][j][k] == typePore)
            numFluidBeads += nFluidsPerLattice;
          else
          {
            LATTICE[i][j][k] = typeWall;
            numSolidBeads += nSolidsPerLattice;
          }
        }

      for (unsigned long i = 1; i <= xhi; i++)
        for (unsigned long j = 1; j <= yhi; j++)
          for (unsigned long k = (zlo+2); k <= (zhi-1); k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      numSolidBeads += nSolidsPerLattice;
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if (LATTICE[i][j][k] == typePore)
              numFluidBeads += nFluidsPerLattice;
          }

      break;

    default:
      std::cout << "Error: flowDirection = " << flowDirection << "\n";
      std::exit(0);
      break;
  }

  numSolidVox = numSolidLattice / hx / hy / hz;
  numFluidVox = numFluidLattice / hx / hy / hz;

  // ===========================
  // Output the LAMMPS data file
  // ===========================

  outSolidFile.open(outSolidFileName, std::ios::out);
  outSolidFile.precision(16);
  outSolidFile << std::scientific;

  outSolidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numSolidBeads << " atoms\n"
    << nTypes << " atom types\n"
    << "\n"
    << xlo << " " << xhi << " xlo xhi\n"
    << ylo << " " << yhi << " ylo yhi\n"
    << zlo << " " << zhi << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  outFluidFile.open(outFluidFileName, std::ios::out);
  outFluidFile.precision(16);
  outFluidFile << std::scientific;

  outFluidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numFluidBeads << " atoms\n"
    << nTypes << " atom types\n"
    << "\n"
    << xlo << " " << xhi << " xlo xhi\n"
    << ylo << " " << yhi << " ylo yhi\n"
    << zlo << " " << zhi << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  // Obtain a seed for the random number engine
  // Standard mersenne_twister_engine seeded with rd()

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.01, 0.99);

  unsigned long indexSolidBead = 0;
  unsigned long indexFluidBead = 0;

  switch ( flowDirection )
  {
    case 0:

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long j = 1; j <= nyLattice; j++)
          for (unsigned long k = 1; k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
                      {
                        // Each call to dis(gen) generates a new random double

                        //xrnd = double(i-1) + dis(gen);
                        //yrnd = double(j-1) + dis(gen);
                        //zrnd = double(k-1) + dis(gen);
                        xrnd = double(i-1) + 0.5;
                        yrnd = double(j-1) + 0.5;
                        zrnd = double(k-1) + 0.5;

                        indexSolidBead ++;

                        outSolidFile
                          << indexSolidBead << " " << typeWall << " "
                          << xrnd << " " << yrnd << " " << zrnd << "\n";
                      }
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if ((LATTICE[i][j][k] == typePore) && outFluid)
              for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                //xrnd = double(i-1) + dis(gen);
                //yrnd = double(j-1) + dis(gen);
                //zrnd = double(k-1) + dis(gen);
                xrnd = double(i-1) + 0.5;
                yrnd = double(j-1) + 0.5;
                zrnd = double(k-1) + 0.5;

                indexFluidBead ++;

                outFluidFile
                  << indexFluidBead << " " << typePore << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
          }

      break;

    // In xlo and xhi planes:
    case 1:

      for (unsigned long j = 1; j <= nyLattice; j++)
        for (unsigned long k = 1; k <= nzLattice; k++)
        {
          unsigned long i;

          i = hx;
          if (LATTICE[i][j][k] == typeWall)
          {
            for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
            {
              // Each call to dis(gen) generates a new random double

              xrnd = double(i-1) + dis(gen);
              yrnd = double(j-1) + dis(gen);
              zrnd = double(k-1) + dis(gen);

              indexSolidBead ++;

              outSolidFile
                << indexSolidBead << " " << typeWall << " "
                << xrnd << " " << yrnd << " " << zrnd << "\n";
            }
          }

          i = nxLattice-hx+1;
          if (LATTICE[i][j][k] == typeWall)
          {
            for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
            {
              // Each call to dis(gen) generates a new random double

              xrnd = double(i-1) + dis(gen);
              yrnd = double(j-1) + dis(gen);
              zrnd = double(k-1) + dis(gen);

              indexSolidBead ++;

              outSolidFile
                << indexSolidBead << " " << typeWall << " "
                << xrnd << " " << yrnd << " " << zrnd << "\n";
            }
          }
        }

      for (unsigned long i = hx+1; i <= (nxLattice-hx); i++)
        for (unsigned long j = 1; j <= nyLattice; j++)
          for (unsigned long k = 1; k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
                      {
                        // Each call to dis(gen) generates a new random double

                        //xrnd = double(i-1) + dis(gen);
                        //yrnd = double(j-1) + dis(gen);
                        //zrnd = double(k-1) + dis(gen);
                        xrnd = double(i-1) + 0.5;
                        yrnd = double(j-1) + 0.5;
                        zrnd = double(k-1) + 0.5;

                        indexSolidBead ++;

                        outSolidFile
                          << indexSolidBead << " " << typeWall << " "
                          << xrnd << " " << yrnd << " " << zrnd << "\n";
                      }
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if ((LATTICE[i][j][k] == typePore) && outFluid)
              for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                //xrnd = double(i-1) + dis(gen);
                //yrnd = double(j-1) + dis(gen);
                //zrnd = double(k-1) + dis(gen);
                xrnd = double(i-1) + 0.5;
                yrnd = double(j-1) + 0.5;
                zrnd = double(k-1) + 0.5;

                indexFluidBead ++;

                outFluidFile
                  << indexFluidBead << " " << typePore << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
          }

      break;

    // In ylo and yhi planes:
    case 2:

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long k = 1; k <= nzLattice; k++)
        {
          unsigned long j;

          j = hy;
          if (LATTICE[i][j][k] == typeWall)
          {
            for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
            {
              // Each call to dis(gen) generates a new random double

              xrnd = double(i-1) + dis(gen);
              yrnd = double(j-1) + dis(gen);
              zrnd = double(k-1) + dis(gen);

              indexSolidBead ++;

              outSolidFile
                << indexSolidBead << " " << typeWall << " "
                << xrnd << " " << yrnd << " " << zrnd << "\n";
            }
          }

          j = nyLattice-hy+1;
          if (LATTICE[i][j][k] == typeWall)
          {
            for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
            {
              // Each call to dis(gen) generates a new random double

              xrnd = double(i-1) + dis(gen);
              yrnd = double(j-1) + dis(gen);
              zrnd = double(k-1) + dis(gen);

              indexSolidBead ++;

              outSolidFile
                << indexSolidBead << " " << typeWall << " "
                << xrnd << " " << yrnd << " " << zrnd << "\n";
            }
          }
        }

      for (unsigned long i = 1; i <= nxLattice; i++)
        for (unsigned long j = (hy+1); j <= (nyLattice-hy); j++)
          for (unsigned long k = 1; k <= nzLattice; k++)
          {
            if (LATTICE[i][j][k] == typeWall)
            {
              bool fillSolids = false;

              for (unsigned long ii = i-1; ii <= (i+1); ii++)
              {
                for (unsigned long jj = j-1; jj <= (j+1); jj++)
                {
                  for (unsigned long kk = k-1; kk <= (k+1); kk++)
                  {
                    if(LATTICE[ii][jj][kk] == typePore)
                    {
                      fillSolids = true;
                      for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
                      {
                        // Each call to dis(gen) generates a new random double

                        //xrnd = double(i-1) + dis(gen);
                        //yrnd = double(j-1) + dis(gen);
                        //zrnd = double(k-1) + dis(gen);
                        xrnd = double(i-1) + 0.5;
                        yrnd = double(j-1) + 0.5;
                        zrnd = double(k-1) + 0.5;

                        indexSolidBead ++;

                        outSolidFile
                          << indexSolidBead << " " << typeWall << " "
                          << xrnd << " " << yrnd << " " << zrnd << "\n";
                      }
                      break;
                    }
                  }
                  if (fillSolids) break;
                }
                if (fillSolids) break;
              }
            }
            else if ((LATTICE[i][j][k] == typePore) && outFluid)
              for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                //xrnd = double(i-1) + dis(gen);
                //yrnd = double(j-1) + dis(gen);
                //zrnd = double(k-1) + dis(gen);
                xrnd = double(i-1) + 0.5;
                yrnd = double(j-1) + 0.5;
                zrnd = double(k-1) + 0.5;

                indexFluidBead ++;

                outFluidFile
                  << indexFluidBead << " " << typePore << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
          }

      break;

    // In zlo and zhi planes:
    case 3:

      if (zlo == 0 && zhi == nzLattice)
      {
        for (unsigned long i = 1; i <= nxLattice; i++)
          for (unsigned long j = 1; j <= nyLattice; j++)
          {
            unsigned long k;

            k = hz;
            if (LATTICE[i][j][k] == typeWall)
            {
              for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                //xrnd = double(i-1) + dis(gen);
                //yrnd = double(j-1) + dis(gen);
                //zrnd = double(k-1) + dis(gen);
                xrnd = double(i-1) + 0.5;
                yrnd = double(j-1) + 0.5;
                zrnd = double(k-1) + 0.5;

                indexSolidBead ++;

                outSolidFile
                  << indexSolidBead << " " << typeWall << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
            }

            k = nzLattice-hz+1;
            if (LATTICE[i][j][k] == typeWall)
            {
              for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                //xrnd = double(i-1) + dis(gen);
                //yrnd = double(j-1) + dis(gen);
                //zrnd = double(k-1) + dis(gen);
                xrnd = double(i-1) + 0.5;
                yrnd = double(j-1) + 0.5;
                zrnd = double(k-1) + 0.5;

                indexSolidBead ++;

                outSolidFile
                  << indexSolidBead << " " << typeWall << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
            }
          }

        for (unsigned long i = 1; i <= nxLattice; i++)
          for (unsigned long j = 1; j <= nyLattice; j++)
            for (unsigned long k = (hz+1); k <= (nzLattice-hz); k++)
            {
              if (LATTICE[i][j][k] == typeWall)
              {
                bool fillSolids = false;

                for (unsigned long ii = i-1; ii <= (i+1); ii++)
                {
                  for (unsigned long jj = j-1; jj <= (j+1); jj++)
                  {
                    for (unsigned long kk = k-1; kk <= (k+1); kk++)
                    {
                      if(LATTICE[ii][jj][kk] == typePore)
                      {
                        fillSolids = true;
                        for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
                        {
                          // Each call to dis(gen) generates a new random double

                          //xrnd = double(i-1) + dis(gen);
                          //yrnd = double(j-1) + dis(gen);
                          //zrnd = double(k-1) + dis(gen);
                          xrnd = double(i-1) + 0.5;
                          yrnd = double(j-1) + 0.5;
                          zrnd = double(k-1) + 0.5;

                          indexSolidBead ++;

                          outSolidFile
                            << indexSolidBead << " " << typeWall << " "
                            << xrnd << " " << yrnd << " " << zrnd << "\n";
                        }
                        break;
                      }
                    }
                    if (fillSolids) break;
                  }
                  if (fillSolids) break;
                }
              }
              else if ((LATTICE[i][j][k] == typePore) && outFluid)
                for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
                {
                  // Each call to dis(gen) generates a new random double

                  //xrnd = double(i-1) + dis(gen);
                  //yrnd = double(j-1) + dis(gen);
                  //zrnd = double(k-1) + dis(gen);
                  xrnd = double(i-1) + 0.5;
                  yrnd = double(j-1) + 0.5;
                  zrnd = double(k-1) + 0.5;

                  indexFluidBead ++;

                  outFluidFile
                    << indexFluidBead << " " << typePore << " "
                    << xrnd << " " << yrnd << " " << zrnd << "\n";
                }
            }

      }
      else
      {
        for (unsigned long i = 1; i <= xhi; i++)
          for (unsigned long j = 1; j <= yhi; j++)
          {
            unsigned long k;

            k = zlo+1;

            if (LATTICE[i][j][k] == typeWall)
            {
              for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                xrnd = double(i-1) + dis(gen);
                yrnd = double(j-1) + dis(gen);
                zrnd = double(k-1) + dis(gen);

                indexSolidBead ++;

                outSolidFile
                  << indexSolidBead << " " << typeWall << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
            }
            else if ((LATTICE[i][j][k] == typePore) && outFluid)
              for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                xrnd = double(i-1) + dis(gen);
                yrnd = double(j-1) + dis(gen);
                zrnd = double(k-1) + dis(gen);

                indexFluidBead ++;

                outFluidFile
                  << indexFluidBead << " " << typePore << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }

            k = zhi;

            if (LATTICE[i][j][k] == typeWall)
            {
              for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                xrnd = double(i-1) + dis(gen);
                yrnd = double(j-1) + dis(gen);
                zrnd = double(k-1) + dis(gen);

                indexSolidBead ++;

                outSolidFile
                  << indexSolidBead << " " << typeWall << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
            }
            else if ((LATTICE[i][j][k] == typePore) && outFluid)
              for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
              {
                // Each call to dis(gen) generates a new random double

                xrnd = double(i-1) + dis(gen);
                yrnd = double(j-1) + dis(gen);
                zrnd = double(k-1) + dis(gen);

                indexFluidBead ++;

                outFluidFile
                  << indexFluidBead << " " << typePore << " "
                  << xrnd << " " << yrnd << " " << zrnd << "\n";
              }
          }

        for (unsigned long i = 1; i <= xhi; i++)
          for (unsigned long j = 1; j <= yhi; j++)
            for (unsigned long k = (zlo+2); k <= (zhi-1); k++)
            {
              if (LATTICE[i][j][k] == typeWall)
              {
                bool fillSolids = false;

                for (unsigned long ii = i-1; ii <= (i+1); ii++)
                {
                  for (unsigned long jj = j-1; jj <= (j+1); jj++)
                  {
                    for (unsigned long kk = k-1; kk <= (k+1); kk++)
                    {
                      if(LATTICE[ii][jj][kk] == typePore)
                      {
                        fillSolids = true;
                        for (unsigned int irnd = 0; irnd < nSolidsPerLattice; irnd++)
                        {
                          // Each call to dis(gen) generates a new random double

                          xrnd = double(i-1) + dis(gen);
                          yrnd = double(j-1) + dis(gen);
                          zrnd = double(k-1) + dis(gen);

                          indexSolidBead ++;

                          outSolidFile
                            << indexSolidBead << " " << typeWall << " "
                            << xrnd << " " << yrnd << " " << zrnd << "\n";
                        }
                        break;
                      }
                    }
                    if (fillSolids) break;
                  }
                  if (fillSolids) break;
                }
              }
              else if ((LATTICE[i][j][k] == typePore) && outFluid)
                for (unsigned int irnd = 0; irnd < nFluidsPerLattice; irnd++)
                {
                  // Each call to dis(gen) generates a new random double

                  xrnd = double(i-1) + dis(gen);
                  yrnd = double(j-1) + dis(gen);
                  zrnd = double(k-1) + dis(gen);

                  indexFluidBead ++;

                  outFluidFile
                    << indexFluidBead << " " << typePore << " "
                    << xrnd << " " << yrnd << " " << zrnd << "\n";
                }
            }
      }
      break;
  }

  std::cout << "After voxel swap at the region boundaries:\n\n";
  std::cout
    << std::left << std::setfill('.')
    << std::setw(40) << "numSolidVox  "    << "  " << numSolidVox    << "\n"
    << std::setw(40) << "numFluidVox  "    << "  " << numFluidVox    << "\n"
    << std::setw(40) << "numSolidBeads  "  << "  " << numSolidBeads  << "\n"
    << std::setw(40) << "indexSolidBead  " << "  " << indexSolidBead << "\n"
    << std::setw(40) << "numFluidBeads  "  << "  " << numFluidBeads  << "\n"
    << std::setw(40) << "indexFluidBead  " << "  " << indexFluidBead << "\n";

/*
  // =========================================================================
  // Read voxel data for the 1st time to obtain the number of beads for output
  // =========================================================================

  inpSolidFile.open(inpSolidFileName, std::ios::in);

  if (!inpSolidFile.is_open())
  {
    std::cout << "\n Fatal: file " << inpSolidFileName << " does not exist.\n";
    std::exit(0);
  }

  myVoxelIndex = 0;

  std::cout << "1st time looping through input file " << inpSolidFileName << " to check data sanity\n";

  while (!inpSolidFile.eof())
  {
    inpSolidFile >> ix >> iy >> iz >> ivoxel; std::getline(inpSolidFile,skipLine);

    if (inpSolidFile.eof()) break;

    myVoxelIndex ++;
    std::cout << std::setw(41) << "\r  IMG: voxels looped  " << "  " << myVoxelIndex << std::flush;

    // check bound

    if (!checkRange) continue;
    if (ix > imgNx || iy > imgNy || iz > imgNz)
    {
      std::cout
        << "\n\n" << "Error: in file " << inpSolidFileName << ", Line " << myVoxelIndex << ", "
        << "(ix,iy,iz) = (" << ix << "," << iy << "," << iz << ") out of range!\n";
      std::exit(0);
    }

    // check ROI

    bool isROI = true;

    if (ix < roiNxLo || ix > roiNxHi) isROI = false;
    if (iy < roiNyLo || iy > roiNyHi) isROI = false;
    if (iz < roiNzLo || iz > roiNzHi) isROI = false;

    // count numbers

    if (ivoxel == 0)
    {
      numSolidVoxels ++;
      if (isROI)
      {
        numVoxelsROI ++;
        numSolidVoxelsROI ++;
        numSolidBeads += nSolidsPerLattice;
        numBeads      += nSolidsPerLattice;
      }
    }
    else if (ivoxel < nTypes)
    {
      numFluidVoxels ++;
      if (isROI)
      {
        numVoxelsROI ++;
        numFluidVoxelsROI ++;
        numFluidBeads += nFluidsPerLattice;
        numBeads      += nFluidsPerLattice;
      }
    }
    else if (ivoxel >= nMaxTypes)
    {
      std::cout
        << "\n Error: in file " << inpSolidFileName << "\n"
        << "  ix = " << ix << ", iy = " << iy << ", iz = " << iz << ", ivoxel = " << ivoxel << "\n"
        << "contains an invalid voxel value ...\n";
      std::exit(0);
    }
  }
  std::cout << "\n";
  inpSolidFile.close();

  // Check data consistancy

  if (myVoxelIndex != numVoxels)
  {
    std::cout
      << "\n" << "Error: myVoxelIndex not equal to numVoxels!\n"
      << "myVoxelIndex = " << myVoxelIndex << "\n"
      << "numVoxels    = " << numVoxels    << "\n";
    std::exit(0);
  }

  numBeads = numSolidBeads + numFluidBeads;

  std::cout << std::setw(40) << "  IMG: solid voxels  " << "  " << numSolidVoxels << "\n";
  std::cout << std::setw(40) << "  IMG: fluid voxels  " << "  " << numFluidVoxels << "\n";
  std::cout << std::setw(40) << "  ROI: total voxels  " << "  " << numVoxelsROI << "\n";
  std::cout << std::setw(40) << "  ROI: solid voxels  " << "  " << numSolidVoxelsROI << "\n";
  std::cout << std::setw(40) << "  ROI: fluid voxels  " << "  " << numFluidVoxelsROI << "\n";
  std::cout << std::setw(40) << "  ROI: total beads  "  << "  " << numBeads << "\n";
  std::cout << std::setw(40) << "  ROI: solid beads  "  << "  " << numSolidBeads << "\n";
  std::cout << std::setw(40) << "  ROI: fluid beads  "  << "  " << numFluidBeads << "\n";

  // ===============================================================
  // Read voxel data for the 2nd time to output the LAMMPS data file
  // ===============================================================

  roixlo *= myScale;
  roixhi *= myScale;
  roiylo *= myScale;
  roiyhi *= myScale;
  roizlo *= myScale;
  roizhi *= myScale;

  outSolidFile.open(outSolidFileName, std::ios::out);
  outSolidFile.precision(16);
  outSolidFile << std::scientific;

  outSolidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numSolidBeads << " atoms\n"
    << nTypes << " atom types\n"
    << "\n"
    << roixlo << " " << roixhi << " xlo xhi\n"
    << roiylo << " " << roiyhi << " ylo yhi\n"
    << roizlo << " " << roizhi << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  outFluidFile.open(outFluidFileName, std::ios::out);
  outFluidFile.precision(16);
  outFluidFile << std::scientific;

  outFluidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numFluidBeads << " atoms\n"
    << nTypes << " atom types\n"
    << "\n"
    << roixlo << " " << roixhi << " xlo xhi\n"
    << roiylo << " " << roiyhi << " ylo yhi\n"
    << roizlo << " " << roizhi << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  // Obtain a seed for the random number engine
  // Standard mersenne_twister_engine seeded with rd()

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.01, 0.99);

  // Re-initialize counters

  myVoxelIndex = 0;
  mySolidBeadIndex = 0;
  myFluidBeadIndex = 0;

  inpSolidFile.open(inpSolidFileName, std::ios::in);

  std::cout << "2nd time looping through input file " << inpSolidFileName << " to write output\n";

  while (!inpSolidFile.eof())
  {
    inpSolidFile >> ix >> iy >> iz >> ivoxel; std::getline(inpSolidFile,skipLine);

    if (inpSolidFile.eof()) break;

    myVoxelIndex ++;
    std::cout << std::setw(41) << "\r  IMG: voxels looped  " << "  " << myVoxelIndex << std::flush;

    bool isROI = true;

    if (ix < roiNxLo || ix > roiNxHi) isROI = false;
    if (iy < roiNyLo || iy > roiNyHi) isROI = false;
    if (iz < roiNzLo || iz > roiNzHi) isROI = false;
    if (!isROI) continue;

    if (ivoxel == 0)
    {
      if (outSolid)
      {
        for (unsigned int i = 1; i <= nSolidsPerLattice; i++)
        {
          mySolidBeadIndex ++;

          // Each call to dis(gen) generates a new random double

          xrnd = myScale * (double(ix - 1) + dis(gen)) * hx;
          yrnd = myScale * (double(iy - 1) + dis(gen)) * hy;
          zrnd = myScale * (double(iz - 1) + dis(gen)) * hz;

          outSolidFile << mySolidBeadIndex << " " << (ivoxel+1) << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
        }
      }
    }
    else if (ivoxel < nTypes)
    {
      if (outFluid)
      {
        for (unsigned int i = 1; i <= nFluidsPerLattice; i++)
        {
          myFluidBeadIndex ++;

          // Each call to dis(gen) generates a new random double

          xrnd = myScale * (double(ix - 1) + dis(gen)) * hx;
          yrnd = myScale * (double(iy - 1) + dis(gen)) * hy;
          zrnd = myScale * (double(iz - 1) + dis(gen)) * hz;

          outFluidFile << myFluidBeadIndex << " " << (ivoxel+1) << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
        }
      }
    }
    else if (ivoxel >= nMaxTypes)
    {
      std::cout
        << "\n Error: in file " << inpSolidFileName << "\n"
        << "  ix = " << ix << ", iy = " << iy << ", iz = " << iz << ", ivoxel = " << ivoxel << "\n"
        << "contains an invalid voxel value ...\n";
      std::exit(0);
    }
  }
  std::cout << "\n";

*/

  inpSolidFile.close();
  outSolidFile.close();
  outFluidFile.close();

  double timeElapsed = (double)(clock() - timeStart) / CLOCKS_PER_SEC;
  std::cout << "Time elapsed in seconds: " << timeElapsed << "\n";
}
