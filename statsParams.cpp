#include "header.h"
string version = "V1.01 (2020-09-18)";
// Default inputs
int nx = 512, ny = 512;
double Hurst = 0.8;
// Global
double lengthX = 1., lengthY = 1.;
int fRollOff = 2;
double lambdaS = 0.005, lambdaR = 0.5;
//double lambdaS = 0.00000001, lambdaR = lengthX+lengthY;
double rRoughNorm = 1.;

// Derived
int nxH = nx/2, nyH = ny/2, nyHP1 = nyH+1;
double dqx = TWOPI/lengthX, dqy = TWOPI/lengthY;
double dqx2 = dqx*dqx, dqy2 = dqy*dqy;

double msStat(double exponent) {
  double q2Max = TWOPI/lambdaS; q2Max *= q2Max;
  double q2Roll = TWOPI/lambdaR; q2Roll *= q2Roll;
  double normRoll = pow(q2Roll, 1+Hurst);
  double expWeight = exponent+1+Hurst;
  double expPSD = -2-2*Hurst;

  double result = 0;
  for (int iqx = 0; iqx < nx; ++iqx) {
    for (int iqy = 0; iqy < ny; ++iqy) {
      double qx = (nxH-iqx)*dqx;//TWOPI/nx;
      double qy = (nyH-iqy)*dqy;//TWOPI/ny;
      double q2 = qx*qx + qy*qy;
      if (q2 > q2Max || q2==0) continue;
      else if (fRollOff==1) result += pow(q2,expWeight)*pow(sqrt(1. + (q2/q2Roll)), expPSD); // smooth roll-off
      else if (q2 < q2Roll) {
        if (fRollOff==2) result += pow(q2, expWeight);
        // else if (fRollOff==0): result += 0;
      }
      else result += normRoll*pow(q2, exponent);
    }
  }
  return(result);
};

void dumpHelp() {
  cout << "Usage: ./stats.exe [nx ny Hurst lambdaR lambdaS fRollOff rRoughNorm]\n\n";
  cout << "Calculates the surface statistics for the given parameters. These \n";
  cout << "include RMS height, gradient and curvature. Parameters have the same \n";
  cout << "meanings and defaults as in the contMech code. Specifying one as '-'\n";
  cout << "(as well as not specifying it at all) sets it to default\n";
};

void statsF() {
  double q2Max = TWOPI/lambdaS; q2Max *= q2Max;
  double q2Roll = TWOPI/lambdaR; q2Roll *= q2Roll;

  double msHeight = 0, msGrad = 0, msCurv = 0;
  for (Lint k = 0; k<(nx)*nyHP1; ++k) {

    int iqx = k/nyHP1;
    int iqy = k%nyHP1;
    int jqx = abs(iqx-nx); // NEW-FixSelfAffine???
    iqx = (iqx<jqx) ? iqx : jqx; // NEW-FixSelfAffine???
    double q2 = pow(iqx,2)*dqx2 + pow(iqy,2)*dqy2;

    if ( (q2>q2Max) || (q2==0) ) continue;

    double weight = 2; // weightQ not usable, unless elaDim>0 !!!
    if ( (iqy==0) || (iqy==ny/2) ) {weight = 1;}

    double specLoc = weight;
    if (q2>q2Roll) specLoc *= pow(q2/q2Roll, -1. -Hurst);
    msHeight += specLoc;
    msGrad += q2*specLoc;
    msCurv += q2*q2*specLoc;
  }
  cout << "rmsHeight (unnormalized as in CM): " << msHeight << "\n";
  cout << "rmsGrad   (unnormalized as in CM): " << msGrad << "\n";
  cout << "rmsCurv   (unnormalized as in CM): " << msCurv << "\n\n";
  cout << "rmsHeight (calculated as in CM):   " << sqrt(msHeight/msGrad) << "\n";
  cout << "rmsGrad   (calculated as in CM):   " << sqrt(msGrad/msGrad) << "\n";
  cout << "rmsCurv   (calculated as in CM):   " << sqrt(msCurv/msGrad) << "\n\n";
};

int main(int argc, char **argv) {
  cout << "----- Surface Stats from Parameters " << version << " -----\n";

  // Process arguments
  vector<string> args(argc);
  for (int iArg=0; iArg<argc; ++iArg) args[iArg] = argv[iArg];
  if (argc > 1) {
    if (args[1]=="-h" || args[1]=="--help") {dumpHelp(); return(0);}
    else if (args[1]!="-") nx = stoi(args[1]);
  }
  if (argc > 2 && args[2]!="-") ny = stoi(args[2]);
  if (argc > 3 && args[3]!="-") Hurst = stod(args[3]);
  if (argc > 4 && args[4]!="-") lambdaR = stod(args[4]);
  if (argc > 5 && args[5]!="-") lambdaS = stod(args[5]);
  if (argc > 6 && args[6]!="-") fRollOff = stoi(args[6]);
  if (argc > 7 && args[7]!="-") rRoughNorm = stod(args[7]);

  cout << "nx = " << nx << ", ny = " << ny << ", Hurst = " << Hurst << "\n";
  cout << "lambdaR = " << lambdaR << ", lambdaS = " << lambdaS << ", fRollOff = " << fRollOff << "\n\n";
  int nxH = nx/2, nyH = ny/2, nyHP1 = nyH+1;

  // Calculate stats
  double msHeight = msStat(-1.-Hurst);
  double msGrad = msStat(-Hurst);
  double msCurv = msStat(1.-Hurst);

  ofstream output;
  output.open("randRough.dat",ofstream::app);
  output << "# Hurst\trmsHeight\trmsGrad\trmsCurv\n";
  output << Hurst << "\t" << sqrt(msHeight) << "\t" << sqrt(msGrad) << "\t" << sqrt(msCurv) << "\n";
  output.close();

  cout << "msHeight  (unnormalized):          " << msHeight << "\n";
  cout << "msGrad    (unnormalized):          " << msGrad << "\n";
  cout << "msCurv    (unnormalized):          " << msCurv << "\n\n";

  cout << "rmsHeight (normalized by Height):    " << rRoughNorm << "\n";
  cout << "rmsGrad   (normalized by Height):    " << rRoughNorm*sqrt(msGrad/msHeight) << "\n";
  cout << "rmsCurv   (normalized by Height):    " << rRoughNorm*sqrt(msCurv/msHeight) << "\n\n";

  cout << "rmsHeight (normalized by Grad):    " << rRoughNorm*sqrt(msHeight/msGrad) << "\n";
  cout << "rmsGrad   (normalized by Grad):    " << rRoughNorm << "\n";
  cout << "rmsCurv   (normalized by Grad):    " << rRoughNorm*sqrt(msCurv/msGrad) << "\n\n";


  //statsF();
};

