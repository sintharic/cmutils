#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <algorithm>
#include <bitset>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <ctime>

using namespace std;

string version = "V1.0 (2020-09-15)";

typedef ptrdiff_t Lint;
typedef std::complex<double> Complex;

#ifndef _PI_
  #define _PI_
  const double PI = 4.*atan(1.), TWOPI = 2*PI;
#endif

const Lint N1024 = 1024;

template <class T> int sign(T a){
  return(a<0 ? -1 : 1);
}


// Dimensions
double lengthX = 1, lengthY = 1; // Only necessary for lateral information
int nx, ny, nyHP1;
double dx, dy;
Lint sizeReal, sizeCompl;
bool convert_to_mm = 0; // Scales z-values by 1000

// Fourier space coordinates
Lint iqx, iqy;
double q, dqx2, dqy2;
inline void get_iqxiqy(Lint k){ iqx = k/nyHP1; iqy = k%nyHP1; };
double getQ(Lint);

// Input/Output
string inputPath;
void readINM(double *, const string&, int, int, double);
string outputPath = "-";
void dumpImage(double *, const string&, int, int);
double normFac = 0.001;

// Timer
double tStart, tFinish;


int main(int argc, char **argv) {
  cout << "----- INM Image Conversion " << version << " -----\n\n";

  // Process arguments
  vector<string> args(argc);
  for (int iArg=0; iArg<argc; ++iArg) args[iArg] = argv[iArg];
  if (argc > 1) inputPath = args[1];
  if (inputPath == "-h" || inputPath == "--help") {
    cout << "Usage: " << args[0] << " inputPath nx ny [outputPath] [normFac]\n\n";
    cout << "Converts an INM microscope file from flat array notation \n";
    cout << "to matrix/image format without redundant lateral information \n";
    cout << "and saves it to <outputPath> (default: '<inputPath>-Img').\n";
    cout << "The surface can be multiplied by <normFac> (default: 0.001, \n";
    cout << "i.e. µm -> mm). If you want to specify <normFac>, set <outputPath> \n";
    cout << "to '-' to get the default.\n";
    return(0);
  }
  if (argc > 2) nx = stoi(args[2]);
  if (argc > 3) ny = stoi(args[3]);
  if (argc > 4) outputPath = args[4];
  if (outputPath == "-") outputPath = inputPath + "-Img";
  if (argc > 5) normFac = stod(args[5]);


  // Calculate array sizes
  dx = lengthX/nx; dy=dx;
  if (nx != ny) lengthY = dx*ny;
  nyHP1 = (ny/2) + 1;
  sizeCompl = nx*nyHP1;
  sizeReal = 2*sizeCompl;
  dqx2 = TWOPI / lengthX; dqx2 *= dqx2;
  dqy2 = TWOPI / lengthY; dqy2 *= dqy2;

  // Version without FFTW
  vector<double> surfR(sizeReal);
  double *ptrSurfR = surfR.data();

  // Read file
  readINM(ptrSurfR, inputPath, nx, ny, normFac);
  dumpImage(ptrSurfR, outputPath, nx, ny);//
}



double getQ(Lint k) {
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx) ? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
}


void readINM(double* array, const string& fileName, int nxOld, int nyOld, double norm) {

  ifstream configIn(fileName);
  if (!configIn.is_open())  return;

  double dummyX, dummyY, dispZ;
  string dispZstr;

  // Stats of data file
  double minX = 1.e+38, maxX = -1.e+38;
  double minY = 1.e+38, maxY = -1.e+38;
  double minHeight = 1.e+38, maxHeight = -1.e+38;

  // Read Z values from file
  int prog = 0;
  cout << "Importing file:\n";
  tStart = (double) clock();
  bool someNaN = 0;
  for (Lint iyOld = 0; iyOld < nyOld; iyOld++){

    // NEW1 DEBUG
    std::string msg;
    if ( (int)(100*iyOld)/nyOld > prog ) {
      if (prog < 10) msg = "% :";
      else msg = "%:";
      cout << "\r" << prog << msg << flush;
      prog+=5;
    } // ENDNEW DEBUG

    for (Lint ixOld = 0; ixOld < nxOld; ixOld++){

      // Read line
      configIn >> dummyX >> dummyY >> dispZstr;

      //if(ixOld < nx && iyOld < ny) {// HERE
      // X and Y value stats
      if (dummyX < minX) minX = dummyX;
      else if (dummyX > maxX) maxX = dummyX;
      if (dummyY < minY) minY = dummyY;
      else if (dummyY > maxY) maxY = dummyY;

      // convert Z value
      if (dispZstr != "***" && !dispZstr.empty()) {//Fix2020-12-01
        dispZ = stof(dispZstr);

        // Z value stats
        if (dispZ < minHeight) minHeight = dispZ;
        else if (dispZ > maxHeight) maxHeight = dispZ;
      }
      else {
        dispZ = -1.e+38;
        someNaN = 1;
      }
      Lint ii = ixOld*nyOld+iyOld;// HERE
      array[ii] = dispZ;
      //}// HERE
    }
  }
  configIn.close();
  tFinish = (double) clock() - tStart;
  cout << "\r100% (Import took " << tFinish/CLOCKS_PER_SEC << " s)\n";

  lengthX = nx*(maxX - minX)/(nx-1);
  lengthY = ny*(maxY - minY)/(ny-1);
  dx = lengthX/nx;
  dy = lengthY/ny;
  cout << " lengthX: " << lengthX << ", lengthY: " << lengthY << "\n";
  //cout << " minX: " << minX << ", maxX: " << maxX << "\n";
  //cout << " minY: " << minY << ", maxY: " << maxY << "\n";
  cout << " minZ: " << minHeight << ", maxZ: " << maxHeight << " (unnormalized)\n";
  cout << " minZ: " << minHeight*norm << ", maxZ: " << maxHeight*norm << " (normalized)\n";

  // Set NaN values
  if (someNaN) {
    tStart = (double) clock();
    double borderVal = 2*minHeight - maxHeight;
    for (Lint ixOld = 2; ixOld < nxOld-2; ++ixOld) {
      for (Lint iyOld = 2; iyOld < nyOld-2; ++iyOld) {
        Lint ii = ixOld*nyOld+iyOld;// HERE
        if (array[ii] == -1.e+38) {
          // Check neighbors
          int goodNeighbors = 0;
          double meanNeighbors = 0.;
          vector<Lint> neighIdcs = {(ixOld-1)*nyOld+iyOld, (ixOld+1)*nyOld+iyOld, 
                                    ixOld*nyOld+iyOld-1, ixOld*nyOld+iyOld+1,
                                    (ixOld-1)*nyOld+iyOld-1, (ixOld-1)*nyOld+iyOld+1,
                                    (ixOld+1)*nyOld+iyOld-1, (ixOld+1)*nyOld+iyOld+1,
                                    (ixOld-2)*nyOld+iyOld-2, (ixOld-2)*nyOld+iyOld-1,
                                    (ixOld-2)*nyOld+iyOld, (ixOld-2)*nyOld+iyOld+1,
                                    (ixOld+2)*nyOld+iyOld-2, (ixOld+2)*nyOld+iyOld-1,
                                    (ixOld+2)*nyOld+iyOld, (ixOld+2)*nyOld+iyOld+1,
                                    (ixOld-2)*nyOld+iyOld+2, (ixOld+2)*nyOld+iyOld+2, 
                                    (ixOld+1)*nyOld+iyOld-2, (ixOld-1)*nyOld+iyOld-2,
                                    (ixOld+1)*nyOld+iyOld+2, (ixOld-1)*nyOld+iyOld+2,
                                    ixOld*nyOld+iyOld-2, ixOld*nyOld+iyOld+2 };
          for (Lint iNeigh : neighIdcs) {
            if (array[iNeigh]!=-1.e+38 && array[iNeigh]>borderVal) {
              goodNeighbors ++;
              meanNeighbors += array[iNeigh];
            }
          }
          meanNeighbors /= goodNeighbors;
          if (goodNeighbors >=12) array[ii] = meanNeighbors;
          else array[ii] = borderVal;//*/
        }
      }
    }

    // Borders
    for (Lint ixOld = 0; ixOld < nxOld; ++ixOld) {
      for (Lint kyOld : {0,1,nyOld-2,nyOld-1}) {
        Lint ii = ixOld*nyOld + kyOld;
        if (array[ii] == -1.e+38) array[ii] = borderVal;
      }
    }
    for (Lint iyOld = 0; iyOld < nyOld; ++iyOld) {
      for (Lint kxOld : {0,1,nxOld-2,nxOld-1}) {
        Lint ii = kxOld*nyOld + iyOld;
        if (array[ii] == -1.e+38) array[ii] = borderVal;
      }
    }

    tFinish = (double) clock() - tStart;
    cout << "Setting NaN values took " << tFinish/CLOCKS_PER_SEC << " s\n";
  }

  // convert to mm
  for (int k=0; k<nx*ny; ++k) array[k] *= norm;


  // Calculate msHeight, Grad, Curv in Real space
  double msGradRx = 0.0, msCurvRx = 0.0;
  double msGradRy = 0.0, msCurvRy = 0.0;
  double msHeightR = 0.0;
  for (Lint ixOld=0; ixOld < nxOld; ixOld++){
    double msGradRxLoc = 0.0, msCurvRxLoc = 0.0;
    double msGradRyLoc = 0.0, msCurvRyLoc = 0.0;
    double msHeightLoc = 0.0;
    for (Lint iyOld=0; iyOld<nyOld; iyOld++){
      if (ixOld > 0){
        msGradRxLoc += pow(array[ixOld*nyOld+iyOld] - array[(ixOld-1)*nyOld+iyOld],2); // Norm dx missing!
        if (ixOld < nxOld-1) {
          msCurvRxLoc += pow(array[(ixOld+1)*nyOld+iyOld] + array[(ixOld-1)*nyOld+iyOld] - 2*array[ixOld*nyOld+iyOld], 2); // Norm dx2 missing!
        }
      }
      msHeightLoc += array[ixOld*nyOld+iyOld]*array[ixOld*nyOld+iyOld];
      if (iyOld > 0){
        msGradRyLoc += pow(array[ixOld*nyOld+iyOld] - array[ixOld*nyOld+iyOld-1],2); // Norm dy missing!
        if (iyOld < nyOld-1) {
          msCurvRyLoc += pow(array[ixOld*nyOld+iyOld+1] + array[ixOld*nyOld+iyOld-1] - 2*array[ixOld*nyOld+iyOld], 2); // Norm dx2 missing!
        }
      }
    }
    msHeightR += msHeightLoc/nyOld;
    msGradRx += msGradRxLoc/nyOld;
    msGradRy += msGradRyLoc/(nyOld-1);
    msCurvRx += msCurvRxLoc/nyOld;
    msCurvRy += msCurvRyLoc/(nyOld-2);
  }
  msHeightR /= nxOld;
  msGradRx /= (nxOld-1)*dx*dx; // Norm dx
  msGradRy /= nxOld*dy*dy; // Norm dy
  msCurvRx /= (nxOld-2)*dx*dx*dx*dx; // Norm dx2
  msCurvRy /= nxOld*dy*dy*dy*dy; // Norm dy2
  cout << " rmsHeight (imported): " << sqrt(msHeightR) << "\n";
  cout << " rmsGrad   (imported): " << sqrt(msGradRx+msGradRy) << "\n";
  cout << " rmsCurv   (imported): " << sqrt(msCurvRx+msCurvRy) << "\n";
  cout << "(Caution: accuracy of Grad/Curv is limited by discretization!)\n\n";
}

void dumpImage(double* array, const string& fileName, int nxExp, int nyExp) {
  ofstream output;
  output.open(fileName);

  int prog = 0;
  tStart = (double) clock();
  cout << "Exporting image file:\n";
  for (int ix=0; ix<nxExp; ++ix) {
    // NEW1 DEBUG
    std::string msg;
    if ( (int)(100*ix)/nxExp > prog ) {
      if (prog < 10) msg = "% :";
      else msg = "%:";
      cout << "\r" << prog << msg << flush;
      prog+=5;
    } // ENDNEW DEBUG

    for (int iy=0; iy<nyExp; ++iy) {
      output << array[ix*ny+iy] << "\t";
    }
    output << "\n";
  }
  tFinish = (double) clock() - tStart;
  cout << "\r100% (Export took " << tFinish/CLOCKS_PER_SEC << " s)\n";
  output.close();
}