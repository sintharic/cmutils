#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <cstring>

using namespace std;
string version = "V1.0 (2020-09-15)";

typedef ptrdiff_t Lint;

// Default values
string inputPath = "equilPos0.dat", outputPath = "-";
double xStart = 0., xEnd = 1., yStart = 0., yEnd = 1.;


int main(int argc, char **argv) {
  cout << "----- Surface Range Extraction " << version << " -----\n\n";
  // Process input arguments
  vector<string> args(argc);
  for (int iArg=0; iArg<argc; ++iArg) args[iArg] = argv[iArg];
  if (argc > 1) inputPath = args[1];
  if (inputPath == "--help" || inputPath == "-h") {
    cout << "Usage: " << args[0] << " infile xStart xEnd yStart yEnd [outfile]\n\n";
    cout << "Extracts from <infile> the surface points in the \n";
    cout << "interval given by <xStart>, <xEnd>, <yStart>, <yEnd> \n";
    cout << "and dumps it to [outfile] (default: '<infile>-Red').\n";
    return(0);
  }
  if (argc > 2) xStart = stod(args[2]);
  if (argc > 3) xEnd = stod(args[3]);
  if (argc > 4) yStart = stod(args[4]);
  if (argc > 5) yEnd = stod(args[5]);
  else {cout << "ERROR: Not enough input arguments.\n"; return(1);}
  if (argc > 6) outputPath = args[6];
  if (outputPath == "-") outputPath = inputPath + "-Red";

  // Open files
  ifstream configIn(inputPath);
  if ( (!configIn.is_open()) ) return(1);
  ofstream configOut(outputPath);

  // Surface dimensions
  Lint nxOld, nyOld;
  char dummyChar; //used to hold the "#" character in the input file
  configIn >> dummyChar >> nxOld >> nyOld; // First line contains dimensions
  Lint nyOldM2 = nyOld - 2;

  // Usable for files exported with dumpReal() [elaDim=0] or with dumpReal2() [elaDim=1]
  string ROL;

  // Read file
  Lint nxnyNew = 0, nyNew = 0;
  double xOld = -1.e38, xNow = -1.e37, yNow = -1.e36;
  while (configIn.peek() != EOF) {
    configIn >> xNow >> yNow;

    if (xNow < xStart) {
      // Skip nyOldM2 lines
      for (Lint iyOld=0; iyOld<nyOldM2; ++iyOld) getline(configIn,ROL);
      xOld = xNow;
      continue;
    }
    if (xNow > xEnd) break;

    if (yNow < yStart) {
      getline(configIn,ROL);
      continue;
    }
    if (yNow > yEnd) {
      getline(configIn,ROL);
      continue;
    }

    // xNow and yNow are in the interval
    if (xNow != xOld) {
      configOut << endl;
      nyNew = 0;
    }
    nyNew += 1;
    nxnyNew += 1;

    // Print to file
    getline(configIn,ROL);
    configOut << xNow << "\t" << yNow << ROL << "\n";
    xOld = xNow;
  }
  cout << "New file dimensions (copy to first line if needed):\n";
  cout << "#" << 1.*nxnyNew/nyNew << "\t" << nyNew << "\n\n";
  cout << "Result exported to " << outputPath << "\n";

  configIn.close();
  configOut.close();
};
