#include "header.h"
string version = "V1.00 (2021-11-09)";


// How to use: 
// - set parameters in section "//PARAMETERS"
// - refer to "//CHOOSE" to choose which GF you want to want.
// - note: GreenZX = -GreenXZ


//PARAMETERS
double lengthX = 1., lengthY = 1.;
Lint nx = 256, ny = 256;

int fThickness0 = 2;
const double POISSON = 0.5;
double elastExpnt0 = 1, stiffness0 = 0.5, thickness0 = 0.05;//0.159155;

//TEST: This is how dumpReal() should work!
vector<double*> fields;

// derived variables
Lint nxny = nx*ny;
Lint nxH = nx/2, nyHP1 = ny/2 + 1;
Lint nFour = nx*nyHP1;
Lint nReal = 2*nFour;
Lint irx, iry, iqx, iqy;
double dx = lengthX / nx, dqx = 2*M_PI / lengthX, dqx2 = dqx*dqx;
double dy = lengthY / ny, dqy = 2*M_PI / lengthY, dqy2 = dqy*dqy;


// functions
inline void get_iqxiqy(Lint k){ iqx = k/nyHP1; iqy = k%nyHP1; };
inline void get_irxiry(Lint k){ irx = k/ny; iry = k%ny; };
inline Lint getLinRIndex(Lint ix, Lint iy) {return(ix*ny+iy);};

double getQ(Lint k){
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx) ? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
}


// fftw declarations
Complex* dispF_XX    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
Complex* dispF_XY    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
Complex* dispF_XZ    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
Complex* dispF_ZZ    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
Complex* fieldFFFT = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
double* dispR     = (double *) fftw_malloc( nReal*sizeof(double) );
fftw_plan dispF2R = fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, dispR, FFTW_ESTIMATE);
//fftw_plan dispF2R = fftw_plan_dft_c2r(ny, (fftw_complex*) fieldFFFT, dispR, FFTW_ESTIMATE);

void dumpDisp(string filename){
  ofstream output;
  output.open(filename);
  for (Lint k = 0; k < nxny; ++k) {
    get_irxiry(k);
    //output << iry*dy << "\t" << dispR[k] << "\n";
    double x = (irx-nx/2)*dx, y = (iry-ny/2)*dy; 
    Lint ixShift = (irx+nx/2)%nx;
    Lint iyShift = (iry+ny/2)%ny;
    Lint kShift = getLinRIndex(ixShift,iyShift);
    if (k%ny==0) output << endl; // for gnuplot
    output << x << "\t" << y << "\t" << dispR[kShift] << "\n";
  }
  output.close();
}


int main() {
  cout << "----- 3D Green's functions " << version << " -----\n\n";
  cout << " Refer to comments in the source code for options.\n";

  // Construct Green's tensor
  double GreenXX, GreenXY, GreenZZ, ImGreenXZ;
  //Complex GreenXZ;//TODO-remove
  for (Lint k = 1; k < nFour; ++k) {
    double q = getQ(k);
    get_iqxiqy(k);
    Lint jqx = abs(iqx-nx);
    iqx = (iqx<jqx) ? iqx : -jqx;
    double q_x = iqx*dqx, q_y = iqy*dqy;
    double cosFac = q_x/q, sinFac = q_y/q;

    // from Menga. Int. J. Solids Struct. 164, 212-220 (2019)

    // semi-infinite
    double preFacXX = 1;
    double preFacXY = cosFac*sinFac/(POISSON-1);//TODO: check
    double preFacXZ = 0.5/(1.-POISSON);
    double preFacZZ = 1;
    if (fThickness0==0) {
      preFacXY *= POISSON;
      preFacXZ *= (2*POISSON-1);
    }

    // finite thickness, bottom boundary displacement = 0
    else if (fThickness0==2) {
      double width = thickness0*q;
      double width2 = width*width;

      if (width < 40) {
        double normG = 5+2*width2-4*POISSON*(3-2*POISSON)+(3-4*POISSON)*cosh(2*width);

        preFacXX *= +2*width + (3-4*POISSON)*sinh(2*width); //Menga2019 (A.9)
        preFacXX /= normG;

        //TODO: check
        preFacXY *= 2*width*(POISSON-1) + (3-4*POISSON)*POISSON*sinh(2*width) + 2*(width2 + (1-2*POISSON)*(1-2*POISSON))*tanh(width); //Menga/Scaraggi?
        preFacXY /= normG;

        preFacXZ *= +3 + 2*width2 - 2*POISSON*(5-4*POISSON) - (3-2*POISSON*(5-4*POISSON))*cosh(2*width); //Menga2019 (B.3)
        preFacXZ /= normG;

        preFacZZ *= -2*width + (3-4*POISSON)*sinh(2*width); //Carbone2008 (B.5)
        preFacZZ /= normG;
      }
      else {
        preFacXY *= POISSON;
        preFacXZ *= (2*POISSON-1);
      }

    }

    // Green's function with width-dependet q-Factor
    GreenXX = preFacXX * pow(q,-elastExpnt0) / stiffness0;
    GreenXY = preFacXY * pow(q,-elastExpnt0) / stiffness0;//TODO: check
    //GreenXZ = Complex(0, preFacXZ * pow(q,-elastExpnt0) / stiffness0);//TODO-remove
    ImGreenXZ = preFacXZ * pow(q,-elastExpnt0) / stiffness0;
    GreenZZ = preFacZZ * pow(q,-elastExpnt0) / stiffness0;

    //CHOOSE
    dispF_XX[k] = GreenXX;
    dispF_XY[k] = GreenXY;
    dispF_XZ[k] = Complex(0,1)*ImGreenXZ;
    dispF_ZZ[k] = GreenZZ;
  }

  // dump Fourier coefficients
  ofstream output;
  output.open("Fourier3D.dat");
  output << "# q\tFTGreenXX.real()\tFTGreenXY.real()\tFTGreenXZ.imag()\tFTGreenZZ.real()\n";
  for (Lint k = 0; k < nFour; ++k) {
    double q = getQ(k);
    output << q << "\t" << dispF_XX[k].real() << "\t" << dispF_XY[k].real() << "\t";
    output << dispF_XZ[k].imag() << "\t" << dispF_ZZ[k].real() << "\n";
  }
  output.close();

  // transform XX to real space and dump result
  memcpy(fieldFFFT, dispF_XX, nFour*sizeof(Complex));
  fftw_execute(dispF2R);
  dumpDisp("GreenXX.dat");

  // transform XY to real space and dump result
  memcpy(fieldFFFT, dispF_XY, nFour*sizeof(Complex));
  fftw_execute(dispF2R);
  dumpDisp("GreenXY.dat");

  // transform XZ to real space and dump result
  memcpy(fieldFFFT, dispF_XZ, nFour*sizeof(Complex));
  fftw_execute(dispF2R);
  dumpDisp("GreenXZ.dat");

  // transform XZ to real space and dump result
  memcpy(fieldFFFT, dispF_ZZ, nFour*sizeof(Complex));
  fftw_execute(dispF2R);
  dumpDisp("GreenZZ.dat");

  /*//TEST-Start: This is how dumpReal() should work!
  double* dispDummy = (double *) fftw_malloc( nReal*sizeof(double) );
  for (int k = 0; k < nxny; ++k) dispDummy[k] = 6.666e+66;
  fields.push_back(dispR);
  fields.push_back(dispDummy);
  string header = "# dispZR";
  cout << header + "\tdummy" << "\n";
  for (int k = 0; k < 20; ++k) {
    cout << fields[0][k];
    for (int iField = 1; iField < fields.size(); ++iField) {
      cout << "\t" << fields[iField][k];
    }
    cout << "\n";
  }
  cout << "pointer sizes: " << sizeof(dispR) << " " << sizeof(dispDummy) << " " << sizeof(dispF_XZ) << "\n";
  cout << "vector of pointers size: " << sizeof(fields) << "\n";
  cout << "data array sizes: " << nReal*sizeof(double) << " " << nFour*sizeof(Complex) << "\n";//TEST-End*/
}
