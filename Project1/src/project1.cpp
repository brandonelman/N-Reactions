/*  This program solves the radial equation for
 *  Be-10 + n scattering. */ 
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
//#include "lib.h"
using namespace  std;
// output file as global variable
ofstream ofile;
// function declarations
void derivatives(double, double*,  double*,
                    int, double);
void initialize ( double&, int&, int&);
void output( double, double *);
void runge_kutta_4(double *, double *, int, double, double, 
    double *, int, double, void (*)(double, double *, double *, int, double));

//This function calculates the value of the woods-saxon potential
//at some given R in fm
double calculateWS(double R){
  //Potential Constants
  const double a_ws   = 0.65; //fm 
  const double V0     = -61.1; //MeV
  const double R_WS   = 1.2*pow(10,1./3.);//fm

  return  (V0/(1+exp((R-R_WS)/a_ws)));
}

//Pass in initial wave_function and derivative as u[0], u[1]
//respectively and A the matching radius. Make sure the u
//values correspond to u(a)
double calculateRMatrix(double *u, double A){
  return (1/A)*u[0]/u[1];
}

//H[0] will be set to H+, while H[1] will be set to H-
void  calculateH(std::complex<double>* H, double r, int L, double E){
  std::complex<double> i(0,1);
  const double MU2H        = 0.0478450; // 1/(MeV*fm^2)
  double k = sqrt(MU2H*E);//1/fm
  double rho = k*r;

  H[0] = pow(i,-L)*exp(i*rho); 
  H[1] = pow(i,L)*exp(-i*rho); 
}

void  calculateH_Prime(std::complex<double>* H_prime, double r, int L, double E){
  const double MU2H        = 0.0478450; // 1/(MeV*fm^2)
  std::complex<double> i(0,1);
  double k = sqrt(MU2H*E);//1/fm
  double rho = k*r;
  H_prime[0] = k*pow(i,-L+1)*exp(i*rho); 
  H_prime[1] = -k*pow(i,L+1)*exp(-i*rho); 
}

std::complex<double> calculateSMatrix(double *u, const double R_N, double r, int L, double E){
  std::complex<double> H[2];
  std::complex<double> H_prime[2];
  calculateH(H, r,L,E); 
  calculateH_Prime(H_prime, r,L,E);
  double R_L = calculateRMatrix(u,R_N);
  //Calculated from equation 3.1.30 in the book. Note R_N = a
  std::complex<double> S = (H[1]-R_N*R_L*H_prime[1]) / (H[0]-R_N*R_L*H_prime[0]);
  return S;
}

std::complex<double> calculatePhaseShift(std::complex<double> S_L){
  std::complex<double> i(0,1);
  std::complex<double> phase_shift = (1./(2.*i))*log(S_L);
  return phase_shift;
}

int main(int argc, char* argv[])
{
  //  declarations of variables
  double *u, *dudr, *yout, r, step_size, r_max, E0;
  double initial_u, initial_u_prime;
  int i, number_of_steps, num_diff_eq;
  char *outfilename;

  const double R_N = 100;      //Range of potential
  const int NUM_L_VALUES = 1;
  //const double   A = R_N; //Matching point outside potential

  int L; //angular momentum quantum number
  double E; //Energy (MeV)

  double r_matrix[NUM_L_VALUES];
  std::complex<double> s_matrix[NUM_L_VALUES];
  std::complex<double> phase_shifts[NUM_L_VALUES];

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    return 5;
    //    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename);

  //  this is the number of differential equations  
  num_diff_eq = 2;     
  //  allocate space in memory for the arrays containing the derivatives 
  dudr = new double[num_diff_eq];
  u = new double[num_diff_eq];
  yout = new double[num_diff_eq];
  // read in the initial position, velocity and number of steps 
  initialize (E, L, number_of_steps);
  //  setting initial values, step size and max position r_max  
  step_size = R_N / ( (double) number_of_steps);   // the step size     
  r_max = step_size*number_of_steps;       // the final time    

  r = R_N/number_of_steps; // initial position      
  //Based on equation 3.1.17 in text
  u[0] = pow(r,L+1);       // initial wavefunction  
  u[1] = (L+1)*pow(r,L);   // initial wf derivative  

  // now we start solving the differential equations using the RK4 method 
  while (r <= R_N){
    derivatives(r, u, dudr, L, E);   // initial derivatives              
    runge_kutta_4(u, dudr, num_diff_eq, r, step_size, 
                  yout, L, E, derivatives); 
    for (i = 0; i < num_diff_eq; i++) {
      u[i] = yout[i];  
    }
    r += step_size;
    output(r, u);   // write to file 
  }

  r_matrix[L] = calculateRMatrix(u, R_N);
  s_matrix[L] = calculateSMatrix(u, R_N, r,L,E);
  phase_shifts[L] = calculatePhaseShift(s_matrix[L]);


  std::cout << "For L = " << L << " and E = " << E  << ": " << std::endl
            << "R Matrix is: " << r_matrix[L] << std::endl
            << "S Matrix is: " << s_matrix[L] << std::endl
            << "Phase Shift is: " << phase_shifts[L] << std::endl;


  delete [] u; delete [] dudr; delete [] yout; 
  ofile.close();  // close output file
  return 0;
}   //  End of main function 

//     Read in from screen the number of steps,
//     initial position and initial speed 
void initialize (double& E, int& L, int& number_of_steps)
{
  cout << "Number of steps: ";
  cin >> number_of_steps;

  cout << "Energy (MeV): ";
  cin >> E;
  cout << "L-Value: ";
  cin >> L;
}  // end of function initialize  

//   this function sets up the derivatives for this special case  
void derivatives(double r, double *u, double *dudr,
                 int    L, double E){
  //du^2/dr^2 = [L(L+1)/R^2+(2mu/h^2)*(V(R)-E)]u
  //Radial equation constants
  const double MU2H        = 0.0478450; // 1/(MeV*fm^2)
  dudr[0]=  u[1];    // derivative of u 
  dudr[1]= (L*(L+1)/pow(r,2) + MU2H*(calculateWS(r) - E))*u[0];    // derivative of   u' 
} // end of function derivatives  

//    function to write out the final results
void output(double r, double *u)
{
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << r;
  ofile << setw(15) << setprecision(8) << u[0];
  ofile << setw(15) << setprecision(8) << u[1] << std::endl;
}  // end of function output

/*   This function upgrades a function y (input as a pointer)
     and returns the result yout, also as a pointer. Note that
     these variables are declared as arrays.  It also receives as
     input the starting value for the derivatives in the pointer
     dydx. It receives also the variable n which represents the 
     number of differential equations, the step size h and 
     the initial value of x. It receives also the name of the
     function *derivs where the given derivative is computed
 */
void runge_kutta_4(double *y, double *dydx, int n, double x, double h, 
    double *yout, int L, double E, void (*derivs)(double, double *, double *, int, double))
{
  int i;
  double      xh,hh,h6; 
  double *dym, *dyt, *yt;
  //   allocate space for local vectors   
  dym = new double [n];
  dyt =  new double [n];
  yt =  new double [n];
  hh = h*0.5;
  h6 = h/6.;
  xh = x+hh;
  for (i = 0; i < n; i++) {
    yt[i] = y[i]+hh*dydx[i];
  }
  (*derivs)(xh,yt,dyt, L, E);     // computation of k2, eq. 3.60   
  for (i = 0; i < n; i++) {
    yt[i] = y[i]+hh*dyt[i];
  }
  (*derivs)(xh,yt,dym, L, E); //  computation of k3, eq. 3.61   
  for (i=0; i < n; i++) {
    yt[i] = y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x+h,yt,dyt, L, E);    // computation of k4, eq. 3.62   
  //      now we upgrade y in the array yout  
  for (i = 0; i < n; i++){
    yout[i] = y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  }
  delete []dym;
  delete [] dyt;
  delete [] yt;
}       //  end of function Runge-kutta 4
