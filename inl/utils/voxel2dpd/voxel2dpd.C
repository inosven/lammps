#include <fstream>
#include <iomanip>
#include <random>
#include <iostream>
#include <time.h>

void print_message()
{
  std::cout << "\n";
  std::cout << "#########################################################################\n";
  std::cout << "\n";
  std::cout << " __     __   ___   __  __  _____   _       ____    ____    ____    ____   \n";
  std::cout << " \\ \\   / /  / _ \\  \\ \\/ / | ____| | |     |___ \\  |  _ \\  |  _ \\  |  _ \\  \n";
  std::cout << "  \\ \\ / /  | | | |  \\  /  |  _|   | |       __) | | | | | | |_) | | | | | \n";
  std::cout << "   \\ V /   | |_| |  /  \\  | |___  | |___   / __/  | |_| | |  __/  | |_| | \n";
  std::cout << "    \\_/     \\___/  /_/\\_\\ |_____| |_____| |_____| |____/  |_|     |____/  \n";
  std::cout << "\n";
  std::cout << "  A mini code to generate initial particle distribution from voxel data\n";
  std::cout << "\n";
  std::cout << "                        Author: Yidong Xia \n";
  std::cout << "\n";
  std::cout << "#########################################################################\n";
  std::cout << "\n";
}

int main(int argc, char **argv)
{
  print_message();

  clock_t timeStart = clock();

  // ========================
  // Parameter initialization
  // ========================

  bool checkRange  = false;
  bool outputSolid = false;
  bool outputFluid = false;

  const unsigned int typeSolidVoxel = 0;
  const unsigned int typeFluidVoxel = 1;

  const unsigned int typeSolid = 1;
  const unsigned int typeFluid = 2;

  unsigned int ix,iy,iz,ivoxel;
  unsigned int imgNx, imgNy, imgNz;
  unsigned int roiNxLo, roiNxHi;
  unsigned int roiNyLo, roiNyHi;
  unsigned int roiNzLo, roiNzHi;

  unsigned int numTypes = 0;
  unsigned int numSolidBeadsPerVoxel = 0;
  unsigned int numFluidBeadsPerVoxel = 0;

  unsigned int numBeads = 0;
  unsigned int numSolidBeads = 0;
  unsigned int numFluidBeads = 0;

  unsigned int numVoxels = 0;
  unsigned int numSolidVoxels = 0;
  unsigned int numFluidVoxels = 0;

  unsigned int numVoxelsROI = 0;
  unsigned int numSolidVoxelsROI = 0;
  unsigned int numFluidVoxelsROI = 0;

  unsigned int myVoxelIndex;
  unsigned int mySolidBeadIndex;
  unsigned int myFluidBeadIndex;

  double myScale = 0.0;
  double h = 0.0;
  double xlo = 0.0;
  double xhi = 0.0;
  double ylo = 0.0;
  double yhi = 0.0;
  double zlo = 0.0;
  double zhi = 0.0;
  double roixlo = 0.0;
  double roixhi = 0.0;
  double roiylo = 0.0;
  double roiyhi = 0.0;
  double roizlo = 0.0;
  double roizhi = 0.0;
  double xrnd, yrnd, zrnd;

  // strings

  std::string controlFileName = "control";
  std::string inputDataFileName = "inp.dat";
  std::string outputSolidFileName = "out_solid.dat";
  std::string outputFluidFileName = "out_fluid.dat";
  std::string skipLine;

  // I/O files

  std::ifstream controlFile;
  std::ifstream inputDataFile;
  std::ofstream outputSolidFile;
  std::ofstream outputFluidFile;

  // ==================================
  // Open and read input control script
  // ==================================

  controlFile.open(controlFileName, std::ios::in);

  if (!controlFile.is_open())
  {
    std::cout
      << "\n" << "Fatal: file " << controlFileName
      << " does not exist.\n";
    std::exit(0);
  }

  // Read the voxel range check

  std::getline(controlFile,skipLine);
  controlFile >> checkRange; std::getline(controlFile,skipLine);

  // Read the scaling factor

  std::getline(controlFile,skipLine);
  controlFile >> myScale; std::getline(controlFile,skipLine);

  if (myScale <= 0.0)
  {
    std::cout
      << "\n" << "Fatal: in file " << controlFileName
      << " myScale = " << myScale << " is not a valid value.\n";
    std::exit(0);
  }

  // Read the physical length of unit resolution in nm

  std::getline(controlFile,skipLine);
  controlFile >> h; std::getline(controlFile,skipLine);

  if (h <= 0.0)
  {
    std::cout
      << "\n" << "Fatal: in file " << controlFileName
      << " h = " << h << " is not a valid value.\n";
    std::exit(0);
  }

  // Read resolution in three directions

  std::getline(controlFile,skipLine);
  controlFile >> imgNx >> imgNy >> imgNz; std::getline(controlFile,skipLine);

  // Read user-defined ROI in three directions

  std::getline(controlFile,skipLine);
  controlFile >> roiNxLo >> roiNxHi; std::getline(controlFile,skipLine);
  std::getline(controlFile,skipLine);
  controlFile >> roiNyLo >> roiNyHi; std::getline(controlFile,skipLine);
  std::getline(controlFile,skipLine);
  controlFile >> roiNzLo >> roiNzHi; std::getline(controlFile,skipLine);

  if (roiNxLo == 0) roiNxLo = 1;
  if (roiNxHi == 0) roiNxHi = imgNx;
  if (roiNxHi > imgNx || roiNxHi <= roiNxLo)
  {
    std::cout
      << "\n" << "Fatal: in file " << controlFileName << "\n"
      << "  roiNxLo = " << roiNxLo << "\n"
      << "  roiNxHi = " << roiNxHi << "\n"
      << "  are not valid values.\n";
    std::exit(0);
  }

  if (roiNyLo == 0) roiNyLo = 1;
  if (roiNyHi == 0) roiNyHi = imgNy;
  if (roiNyHi > imgNy || roiNyHi <= roiNyLo)
  {
    std::cout
      << "\n" << "Fatal: in file " << controlFileName << "\n"
      << "  roiNyLo = " << roiNyLo << "\n"
      << "  roiNyHi = " << roiNyHi << "\n"
      << "  are not valid values.\n";
    std::exit(0);
  }

  if (roiNzLo == 0) roiNzLo = 1;
  if (roiNzHi == 0) roiNzHi = imgNz;
  if (roiNzHi > imgNz || roiNzHi <= roiNzLo)
  {
    std::cout
      << "\n" << "Fatal: in file " << controlFileName << "\n"
      << "  roiNzLo = " << roiNzLo << "\n"
      << "  roiNzHi = " << roiNzHi << "\n"
      << "  are not valid values.\n";
    std::exit(0);
  }

  // Read the number of wall voxels

  std::getline(controlFile,skipLine);
  controlFile >> numVoxels; std::getline(controlFile,skipLine);

  // Read the solid number density

  std::getline(controlFile,skipLine);
  controlFile >> numSolidBeadsPerVoxel; std::getline(controlFile,skipLine);

  // Read the fluid number density

  std::getline(controlFile,skipLine);
  controlFile >> numFluidBeadsPerVoxel; std::getline(controlFile,skipLine);

  // Read the number of types of beads

  std::getline(controlFile,skipLine);
  controlFile >> numTypes; std::getline(controlFile,skipLine);

  // Read whether to output solid

  std::getline(controlFile,skipLine);
  controlFile >> outputSolid; std::getline(controlFile,skipLine);

  // Read whether to output solid

  std::getline(controlFile,skipLine);
  controlFile >> outputFluid; std::getline(controlFile,skipLine);

  controlFile.close();

  // ==================================
  // Screen print of header information
  // ==================================

  xhi = double(imgNx) * h;
  yhi = double(imgNy) * h;
  zhi = double(imgNz) * h;

  roixlo = double(roiNxLo - 1) * h;
  roixhi = double(roiNxHi)     * h;
  roiylo = double(roiNyLo - 1) * h;
  roiyhi = double(roiNyHi)     * h;
  roizlo = double(roiNzLo - 1) * h;
  roizhi = double(roiNzHi)     * h;

  std::cout << std::left << std::setfill('.')
    << std::setw(40) << "IMG: imgNx,imgNy,imgNz  "       << "  " << imgNx << "," << imgNy << "," << imgNz << "\n"
    << std::setw(40) << "IMG: check voxel range  "       << "  " << (checkRange ? "Yes" : "No") << "\n"
    << std::setw(40) << "IMG: h [nm]  "                  << "  " << h << "\n"
    << std::setw(40) << "IMG: xlo,xhi [nm]  "            << "  " << xlo << "," << xhi << "\n"
    << std::setw(40) << "IMG: ylo,yhi [nm]  "            << "  " << ylo << "," << yhi << "\n"
    << std::setw(40) << "IMG: zlo,zhi [nm]  "            << "  " << zlo << "," << zhi << "\n"
    << std::setw(40) << "IMG: voxels  "                  << "  " << numVoxels << "\n"
    << std::setw(40) << "ROI: NxLo,NxHi  "               << "  " << roiNxLo << "," << roiNxHi << "\n"
    << std::setw(40) << "ROI: NyLo,NyHi  "               << "  " << roiNyLo << "," << roiNyHi << "\n"
    << std::setw(40) << "ROI: NzLo,NzHi  "               << "  " << roiNzLo << "," << roiNzHi << "\n"
    << std::setw(40) << "ROI: xlo,xhi [nm]  "            << "  " << roixlo << "," << roixhi << "\n"
    << std::setw(40) << "ROI: ylo,yhi [nm]  "            << "  " << roiylo << "," << roiyhi << "\n"
    << std::setw(40) << "ROI: zlo,zhi [nm]  "            << "  " << roizlo << "," << roizhi << "\n"
    << std::setw(40) << "ROI: number of bead types  "    << "  " << numTypes << "\n"
    << std::setw(40) << "ROI: beads per solid voxel  "   << "  " << numSolidBeadsPerVoxel << "\n"
    << std::setw(40) << "ROI: beads per fluid voxel  "   << "  " << numFluidBeadsPerVoxel << "\n"
    << std::setw(40) << "ROI: length scaling  "          << "  " << myScale << "\n"
    << std::setw(40) << "ROI: output solid beads  "      << "  " << (outputSolid ? "Yes" : "No") << "\n"
    << std::setw(40) << "ROI: output fluid beads  "      << "  " << (outputFluid ? "Yes" : "No") << "\n"
    << "\n\n";

  // =========================================================================
  // Read voxel data for the 1st time to obtain the number of beads for output
  // =========================================================================

  inputDataFile.open(inputDataFileName, std::ios::in);

  if (!inputDataFile.is_open())
  {
    std::cout << "\n Fatal: file " << inputDataFileName << " does not exist.\n";
    std::exit(0);
  }

  myVoxelIndex = 0;

  std::cout << "1st time looping through input file " << inputDataFileName << " to check data sanity\n";

  while (!inputDataFile.eof())
  {
    inputDataFile >> ix >> iy >> iz >> ivoxel; std::getline(inputDataFile,skipLine);

    if (inputDataFile.eof()) break;

    myVoxelIndex ++;
    std::cout << std::setw(41) << "\r  IMG: voxels looped  " << "  " << myVoxelIndex << std::flush;

    // check bound

    if (!checkRange) continue;
    if (ix > imgNx || iy > imgNy || iz > imgNz)
    {
      std::cout
        << "\n\n" << "Error: in file " << inputDataFileName << ", Line " << myVoxelIndex << ", "
        << "(ix,iy,iz) = (" << ix << "," << iy << "," << iz << ") out of range!\n";
      std::exit(0);
    }

    // check ROI

    bool isROI = true;

    if (ix < roiNxLo || ix > roiNxHi) isROI = false;
    if (iy < roiNxLo || iy > roiNxHi) isROI = false;
    if (iz < roiNxLo || iz > roiNxHi) isROI = false;

    // count numbers

    switch(ivoxel)
    {
      case (typeSolidVoxel) :
        numSolidVoxels ++;
        if (isROI)
        {
          numVoxelsROI ++;
          numSolidVoxelsROI ++;
          numSolidBeads += numSolidBeadsPerVoxel;
          numBeads      += numSolidBeadsPerVoxel;
        }
        break;
      case (typeFluidVoxel) :
        numFluidVoxels ++;
        if (isROI)
        {
          numVoxelsROI ++;
          numFluidVoxelsROI ++;
          numFluidBeads += numFluidBeadsPerVoxel;
          numBeads      += numFluidBeadsPerVoxel;
        }
        break;
      default :
        std::cout
          << "\n Error: in file " << inputDataFileName << "\n"
          << "  ix = " << ix << ", iy = " << iy << ", iz = " << iz << ", ivoxel = " << ivoxel << "\n"
          << "contains an invalid voxel value ...\n";
        std::exit(0);
        break;
    }
  }
  std::cout << "\n";
  inputDataFile.close();

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

  outputSolidFile.open(outputSolidFileName, std::ios::out);
  outputSolidFile.precision(16);
  outputSolidFile << std::scientific;

  outputSolidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numSolidBeads << " atoms\n"
    << numTypes << " atom types\n"
    << "\n"
    << roixlo << " " << roixhi << " xlo xhi\n"
    << roiylo << " " << roiyhi << " ylo yhi\n"
    << roizlo << " " << roizhi << " zlo zhi\n"
    << "\n"
    << "Atoms\n"
    << "\n";

  outputFluidFile.open(outputFluidFileName, std::ios::out);
  outputFluidFile.precision(16);
  outputFluidFile << std::scientific;

  outputFluidFile
    << "LAMMPS data file for read_data\n"
    << "\n"
    << numFluidBeads << " atoms\n"
    << numTypes << " atom types\n"
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

  inputDataFile.open(inputDataFileName, std::ios::in);

  std::cout << "2nd time looping through input file " << inputDataFileName << " to write output\n";

  while (!inputDataFile.eof())
  {
    inputDataFile >> ix >> iy >> iz >> ivoxel; std::getline(inputDataFile,skipLine);

    if (inputDataFile.eof()) break;

    myVoxelIndex ++;
    std::cout << std::setw(41) << "\r  IMG: voxels looped  " << "  " << myVoxelIndex << std::flush;

    bool isROI = true;

    if (ix < roiNxLo || ix > roiNxHi) isROI = false;
    if (iy < roiNxLo || iy > roiNxHi) isROI = false;
    if (iz < roiNxLo || iz > roiNxHi) isROI = false;
    if (!isROI) continue;

    switch(ivoxel)
    {
      case (typeSolidVoxel) :
        if (outputSolid)
        {
          for (unsigned int i = 1; i <= numSolidBeadsPerVoxel; i++)
          {
            mySolidBeadIndex ++;

            // Each call to dis(gen) generates a new random double

            xrnd = myScale * (double(ix - 1) * h + dis(gen) * h);
            yrnd = myScale * (double(iy - 1) * h + dis(gen) * h);
            zrnd = myScale * (double(iz - 1) * h + dis(gen) * h);

            outputSolidFile << mySolidBeadIndex << " " << typeSolid << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
          }
        }
        break;
      case (typeFluidVoxel) :
        if (outputFluid)
        {
          for (unsigned int i = 1; i <= numFluidBeadsPerVoxel; i++)
          {
            myFluidBeadIndex ++;

            // Each call to dis(gen) generates a new random double

            xrnd = myScale * (double(ix - 1) * h + dis(gen) * h);
            yrnd = myScale * (double(iy - 1) * h + dis(gen) * h);
            zrnd = myScale * (double(iz - 1) * h + dis(gen) * h);

            outputFluidFile << myFluidBeadIndex << " " << typeFluid << " " << xrnd << " " << yrnd << " " << zrnd << "\n";
          }
        }
        break;
      default :
        std::cout
          << "\n Error: in file " << inputDataFileName << "\n"
          << "  ix = " << ix << ", iy = " << iy << ", iz = " << iz << ", ivoxel = " << ivoxel << "\n"
          << "contains an invalid voxel value ...\n";
        std::exit(0);
        break;
    }
  }
  std::cout << "\n";

  inputDataFile.close();
  outputSolidFile.close();
  outputFluidFile.close();

  double timeElapsed = (double)(clock() - timeStart) / CLOCKS_PER_SEC;
  std::cout << "Time elapsed in seconds: " << timeElapsed << "\n";
}
