#include "header.h"
string version = "V1.00 (2021-11-02)";


// How to use: 
// - set parameters in section "//PARAMETERS"
// - refer to "//CHOOSE" to choose which GF you want to want.
// - note: GreenZX = -GreenXZ


//PARAMETERS
double lengthX = 1., lengthY = 1.;
Lint nx = 1, ny = 1024;

int fThickness0 = 2;
const double POISSON = 0.5;
double elastExpnt0 = 1, stiffness0 = 0.5, thickness0 = 0.159155;


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

double getQ(Lint k){
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx) ? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
}


// fftw declarations
Complex* dispF    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
Complex* fieldFFFT = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
double* dispR     = (double *) fftw_malloc( nReal*sizeof(double) );
fftw_plan dispF2R = fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, dispR, FFTW_ESTIMATE);
//fftw_plan dispF2R = fftw_plan_dft_c2r(ny, (fftw_complex*) fieldFFFT, dispR, FFTW_ESTIMATE);

void dispIFFT(){
  memcpy(fieldFFFT, dispF, nFour*sizeof(Complex));
  fftw_execute(dispF2R);
}


int main() {
  cout << "----- Green's functions " << version << " -----\n\n";
  cout << " Refer to comments in the source code for options.\n";

  // calculate Green's function 
  double GreenXX, GreenZZ;
  Complex GreenXZ;
  for (Lint k = 1; k < nFour; ++k) {
    double q = getQ(k);
    get_iqxiqy(k);
    if(iqx > nxH) iqx -= nx;
    double q_x = iqx*dqx, q_y = iqy*dqy;

    // from Menga. Int. J. Solids Struct. 164, 212-220 (2019)

    // semi-infinite
    double preFacXX = 1;
    double preFacXZ = 0.5/(1.-POISSON);
    double preFacZZ = 1;
    if (fThickness0==0) preFacXZ *= (2*POISSON-1);

    // finite thickness, bottom boundary displacement = 0
    else if (fThickness0==2) {
      double width = thickness0*q;
      double width2 = width*width;

      if (width < 40) {
        double normG = 5+2*width2-4*POISSON*(3-2*POISSON)+(3-4*POISSON)*cosh(2*width);

        preFacXX *= +2*width + (3-4*POISSON)*sinh(2*width); //Menga2019 (A.9)
        preFacXX /= normG;

        preFacXZ *= +3 + 2*width2 - 2*POISSON*(5-4*POISSON) - (3-2*POISSON*(5-4*POISSON))*cosh(2*width); //Menga2019 (B.3)
        preFacXZ /= normG;

        preFacZZ *= -2*width + (3-4*POISSON)*sinh(2*width); //Carbone2008 (B.5)
        preFacZZ /= normG;
      }
      else preFacXZ *= (2*POISSON-1);

    }

    // Green's function with width-dependet q-Factor
    GreenXX = preFacXX * pow(q,-elastExpnt0) / stiffness0;
    GreenXZ = Complex(0,preFacXZ * pow(q,-elastExpnt0) / stiffness0);
    GreenZZ = preFacZZ * pow(q,-elastExpnt0) / stiffness0;

    dispF[k] = -GreenXZ; //CHOOSE
  }

  // dump Fourier coefficients
  ofstream output;
  output.open("Fourier.dat");
  for (Lint k = 0; k < nFour; ++k) {
    double q = getQ(k);
    output << q << "\t" << dispF[k].real() << "\t" << dispF[k].imag() << "\n";
  }
  output.close();

  // transform to real space and dump result
  dispIFFT();
  output.open("Green.dat");
  for (Lint k = 0; k < nxny; ++k) {
    get_irxiry(k);
    //output << iry*dy << "\t" << dispR[k] << "\n";
    output << (iry-ny/2)*dy << "\t" << dispR[(k+ny/2)%ny] << "\n";
  }
  output.close();

}
