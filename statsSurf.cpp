#include "header.h"
string version = "V1.0 (2020-09-15)";

// Dimensions
double lengthX = 1, lengthY = 1;
int nx=0, ny=0, nyHP1;
double dx, dy;
Lint sizeReal, sizeCompl;

// Fourier space parameters
Lint iqx, iqy;
double q, dqx2, dqy2;
inline void get_iqxiqy(Lint k){ iqx = k/nyHP1; iqy = k%nyHP1; };
double getQ(Lint k) {
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx) ? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
};

// Input / Output
string inputPath, outputPath = "-";
void readReal(double *, const string&);
void readINM(double *, const string&, int, int, double);
unsigned short mode = 1; // 1=CM, 2=INM
double normFac = 0.001;

// Timer
double tStart, tFinish;


int main(int argc, char **argv) {
  cout << "----- Surface Analysis (PSD) " << version << " -----\n\n";

  // Process input arguments
  vector<string> args(argc);
  for (int iArg=0; iArg<argc; ++iArg) args[iArg] = argv[iArg];
  if (argc > 1) inputPath = args[1];
  if (inputPath == "--help" || inputPath == "-h") {
    cout << "Usage 1: " << args[0] << " inputPath [outputPath]\n";
    cout << "  Use this mode to import equilPos files from the contMech code.\n\n";
    cout << "Usage 2: " << args[0] << " inputPath nx ny [outputPath] [normFac]\n";
    cout << "  Use this mode to import INM microscope image files.\n\n";
    cout << "The 2D PSD of the surface is calculated and saved to <outputPath> \n";
    cout << "(Default: '<inputPath>-PSD'). The INM surface can be multiplied \n";
    cout << "by <normFac> (default: 0.001, i.e. µm -> mm). If you want to \n";
    cout << "specify <normFac>, set <outputPath> to '-' to get the default.\n";
    return(0);
  }
  if (argc > 3) mode = 2;
  else mode = 1;

  // Mode 1: CM
  if (mode == 1) {
    if (argc > 2) outputPath = args[2];
    else outputPath = inputPath + "-PSD";
    // Read dimensions from file
    ifstream configIn(inputPath);
    if (!configIn.is_open()) {
      cout << "ERROR: File '" << inputPath << "' does not exist.\n";
      return(1);
    }
    char dummyChar; //used to hold the "#" character
    configIn >> dummyChar >> nx >> ny; // First line contains dimensions
    configIn.close();
    dx = lengthX/nx;
    dy = lengthY/ny;
  }

  // Mode 2: INM
  else if (mode == 2) {
    if (argc > 2) nx = stoi(args[2]);
    if (argc > 3) ny = stoi(args[3]);
    if (argc > 4) outputPath = args[4];
    if (argc > 5) normFac = stod(args[5]);
    if (outputPath=="-") outputPath = inputPath + "-PSD";
  }

  if (nx==0 || ny==0) {
    cout << "ERROR: Surface dimensions nx,ny unknown. Run '" << args[0] << " -h'.\n";
    return(1);
  }

  // DEBUG
  cout << inputPath << ":\n";
  cout << " nx: " << nx << ", ny: " << ny << "\n\n";


  // Calculate Fourier dimensions
  nyHP1 = (ny/2) + 1;
  sizeCompl = nx*nyHP1;
  sizeReal = 2*sizeCompl;
  // FFTW memory allocation
  double* surfR = (double *) fftw_malloc( sizeReal*sizeof(double) );
  Complex* surfF = (Complex *) fftw_malloc( sizeCompl*sizeof(Complex) );
  fftw_plan surfR2F = fftw_plan_dft_r2c_2d(nx, ny, surfR, (fftw_complex*) surfF, FFTW_ESTIMATE);
  
  // Read file
  if (mode == 1) readReal(surfR, inputPath);
  else if (mode == 2) readINM(surfR, inputPath, nx, ny, normFac);

  // Calculate surface properties in Fourier space
  cout << "Surface properties in Fourier space:\n";
  fftw_execute(surfR2F);
  
  ofstream output;
  output.open(outputPath);

  double msHeight = 0, msGrad = 0, msCurv = 0;
  dqx2 = TWOPI / lengthX; dqx2 *= dqx2;
  dqy2 = TWOPI / lengthY; dqy2 *= dqy2;
  Lint nxny = nx*ny;
  for (Lint k=0; k<sizeCompl; ++k){
    q = getQ(k);
    double q2 = q*q;
    double specLoc = surfF[k].real()*surfF[k].real() + surfF[k].imag()*surfF[k].imag();
    double weight = 2;
    if ( (iqy==0) || (iqy==ny/2) ) weight = 1; // would otherwise be counted twice

    output << q << "\t" << specLoc/nxny << "\n";
    
    specLoc *= weight;
    msHeight += specLoc;
    specLoc *= q2; msGrad += specLoc;
    specLoc *= q2; msCurv += specLoc;
  }
  output.close();
  msHeight /= nxny; msHeight /= nxny;
  msGrad /= nxny; msGrad /= nxny;
  msCurv /= nxny; msCurv /= nxny;
  cout << " rmsHeight (Fourier): " << sqrt(msHeight) << "\n";
  cout << " rmsGrad   (Fourier): " << sqrt(msGrad) << "\n";
  cout << " rmsCurv   (Fourier): " << sqrt(msCurv/2) << "\n";
  cout << "(Caution: Grad/Curv and PSD are only exact for periodic surfaces!)\n";//*/

  fftw_free(surfR2F);
}


void readReal(double* array, const string& fileName){
  int elaDim = 0; // Only random rough solid surfaces expected for this code

  ifstream configIn(fileName);
  if ( (!configIn.is_open()) ) return; // CM-ReadInput: Now usable regardless of elaDim
  cout << "Importing file:\n";

  Lint nxOld, nyOld;
  double dummy, dummyY, dispZ;
  char dummyChar; //used to hold the "#" character in the input file
  configIn >> dummyChar >> nxOld >> nyOld; // First line contains dimensions
  // Determine lengthX and lengthY
  double minX = 1.e+38, maxX = -1.e+38;
  double minY = 1.e+38, maxY = -1.e+38;

  // BEGINCM-ReadInput
  // Usable for files exported with dumpReal() [elaDim=0] or with dumpReal2() [elaDim=1]
  vector <double> arrayOld(nxOld*nyOld);
  Lint streamPos; // Current position in the stream
  double xVal1, xVal2; // Additional "dummies" to detect extra lines
  double maxHeight_Loc = -1.e+38, minHeight_Loc = 1.e+38;
  for (Lint ixOld = 0; ixOld < nxOld; ixOld++){
    for (Lint iyOld = 0; iyOld < nyOld; iyOld++){
      if (elaDim) configIn >> xVal1 >> dummyY >> dispZ >> dummy; // dumpReal2() has 4th column
      else configIn >> xVal1 >> dummyY >> dispZ; // dumpReal() hasn't
      Lint ii = ixOld*nyOld+iyOld;
      array[ii] = dispZ;
      if (array[ii] > maxHeight_Loc) maxHeight_Loc = array[ii];
      else if (array[ii] < minHeight_Loc) minHeight_Loc = array[ii];
      // Determine lengthX and lengthY
      if (xVal1 < minX) minX = xVal1;
      else if (xVal1 > maxX) maxX = xVal1;
      if (dummyY < minY) minY = dummyY;
      else if (dummyY > maxY) maxY = dummyY;
    }
    // Check for additional line
    streamPos = configIn.tellg();
    configIn >> xVal2;
    configIn.clear();
    configIn.seekg(streamPos, configIn.beg); // Return to where xVal2 was read
    // skip extra line with old x-value
    if (xVal2 == xVal1) {
      if (elaDim) configIn >> dummy >> dummy >> dispZ >> dummy;
      else configIn >> dummy >> dummy >> dispZ;
    }
  }
  // Determine lengthX and lengthY
  lengthX = nx*(maxX - minX)/(nx-1);
  lengthY = ny*(maxY - minY)/(ny-1);
  dx = lengthX/nx;
  dy = lengthY/ny;
  cout << " lengthX: " << lengthX << ", lengthY: " << lengthY << "\n\n";

  configIn.close();
  // ENDCM-ReadInput
}
/*
void readReal(double* array, const string& fileName) {

  ifstream configIn(fileName);
  if (!configIn.is_open())  return;
  cout << "Importing file:\n";

  int nxOld, nyOld;
  double dummy, dummyX, dummyY, dispZ;
  char dummyChar; //used to hold the "#" character in dispZR.old
  configIn >> dummyChar >> nxOld >> nyOld; // First line contains dimensions
  // Determine lengthX and lengthY:
  double minX = 1.e+38, maxX = -1.e+38;
  double minY = 1.e+38, maxY = -1.e+38;

  // Read Z values from file
  for (Lint ixOld = 0; ixOld < nxOld; ixOld++){
    for (Lint iyOld = 0; iyOld < nyOld; iyOld++){
      //configIn >> dummyX >> dummyY >> dispZ >> dummy; // elaDim=1
      configIn >> dummyX >> dummyY >> dispZ; // elaDim=0
      Lint ii = ixOld*nyOld+iyOld;
      array[ii] = dispZ;
      if (dummyX < minX) minX = dummyX;
      else if (dummyX > maxX) maxX = dummyX;
      if (dummyY < minY) minY = dummyY;
      else if (dummyY > maxY) maxY = dummyY;
    }
    //configIn >> dummyX >> dummyY >> dispZ >> dummy; // skip extra line, elaDim=1
    configIn >> dummyX >> dummyY >> dispZ; // skip extra line, elaDim=0
  }
  lengthX = nx*(maxX - minX)/(nx-1);
  lengthY = ny*(maxY - minY)/(ny-1);
  dx = lengthX/nx;
  dy = lengthY/ny;
  cout << " lengthX: " << lengthX << ", lengthY: " << lengthY << "\n\n";
  configIn.close();
}*/



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
      if (dispZstr != "***") {
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
