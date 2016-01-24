/*  This program solves the radial equation for
 *  Be-10 + n scattering. */ 
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#define NUM_E_VALUES 5000
//#include "lib.h"
using namespace  std;
// output file as global variable
// function declarations
void derivatives(double, double*,  double*,
                    int, double);
void initialize (int&);
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

void output_radial(std::ofstream& ofile, double r, double *u){
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << r;
  ofile << setw(15) << setprecision(8) << u[0];
  ofile << setw(15) << setprecision(8) << u[1] << std::endl;
}  // end of function output


void output_phase_shifts(std::ofstream& ofile, std::complex<double> phase_shift[3][NUM_E_VALUES], double E[NUM_E_VALUES]){
//  const int NUM_E_VALUES = 500;
  const int NUM_L_VALUES = 3;
  const int PI = 3.14159265358979323846;

    for (int i = 0; i < NUM_E_VALUES; i++){
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      ofile << setw(15) << setprecision(8) << E[i];
      if (phase_shift[0][i].real() < 0){
        phase_shift[0][i] += PI;
      }
      ofile << setw(15) << setprecision(8) << phase_shift[0][i].real();
      ofile << setw(15) << setprecision(8) << phase_shift[1][i].real();

      if (phase_shift[2][i].real() > 0){
        phase_shift[2][i] -= PI;
      }
      ofile << setw(15) << setprecision(8) << phase_shift[2][i].real() << std::endl;
    }
}  // end of function output


int main(int argc, char* argv[])
{
  const int NUM_DIFF_EQ = 2;
  const double R_N = 100;      //Range of potential
  const int NUM_L_VALUES = 3;
  const double ENERGY_SPACING = 0.001;
//  const int NUM_E_VALUES = 500;
  const int MIN_E = ENERGY_SPACING;
  //  declarations of variables
  double *u, *dudr, *uout, r, step_size;
  double initial_u, initial_u_prime;
  int i, number_of_steps;
  ofstream ofile_energy_0_1_files[NUM_L_VALUES];//0.1 MeV Energy 
  ofstream ofile_energy_3_0_files[NUM_L_VALUES];//3.0 MeV Energy
  ofstream phase_shift_file;

  ofile_energy_0_1_files[0].open("../output/output_L_0_energy_0_1_radial.dat");
  ofile_energy_0_1_files[1].open("../output/output_L_1_energy_0_1_radial.dat");
  ofile_energy_0_1_files[2].open("../output/output_L_2_energy_0_1_radial.dat");

  ofile_energy_3_0_files[0].open("../output/output_L_0_energy_3_0_radial.dat");
  ofile_energy_3_0_files[1].open("../output/output_L_1_energy_3_0_radial.dat");
  ofile_energy_3_0_files[2].open("../output/output_L_2_energy_3_0_radial.dat");

  phase_shift_file.open("../output/phase_shifts.dat");

//  double E_values[NUM_E_VALUES]; //Energy (MeV)

  double E_values[NUM_E_VALUES];
            
  
  int e_index_0_1;
  int e_index_3_0;
               
  for (int i = 0; i < NUM_E_VALUES; i++){
    E_values[i] = MIN_E + ENERGY_SPACING*i; 
    
    if (abs(E_values[i] - 0.1) < 0.00001){
      std::cout << "0.1 is at position" << i << std::endl;
      e_index_0_1 = i;
    }
    else if (abs(E_values[i] - 3.0) < 0.00001){
      std::cout << "3.0 is at position" << i << std::endl;
      e_index_3_0 = i;
    }
  }
  double r_matrix[NUM_L_VALUES][NUM_E_VALUES];
  std::complex<double> s_matrix[NUM_L_VALUES][NUM_E_VALUES];
  std::complex<double> phase_shifts[NUM_L_VALUES][NUM_E_VALUES];
  for (int i = 0; i < NUM_L_VALUES; i++){
    for (int j = 0; j < NUM_E_VALUES; j++){
      s_matrix[i][j] = 0.0;
      phase_shifts[i][j] = 0.0;
      r_matrix[i][j] = 0.0;
    }
  }

  // Read in output file, abort if there are too few command-line arguments
//if( argc <= 1 ){
//  cout << "Bad Usage: " << argv[0] <<
//    " read also output file on same line" << endl;
//  return 5;
//  //    exit(1);
//}
//else{
//  outfilename=argv[1];
//}
  
  initialize(number_of_steps);

  for (int L = 0; L < NUM_L_VALUES; L++){ 
    for (int e_index = 0; e_index < NUM_E_VALUES; e_index++){
      //  this is the number of differential equations  
      //  allocate space in memory for the arrays containing the derivatives 
      double E = E_values[e_index];
      dudr = new double[NUM_DIFF_EQ];
      u = new double[NUM_DIFF_EQ];
      uout = new double[NUM_DIFF_EQ];


      //  setting initial values, step size 
      step_size = R_N / ( (double) number_of_steps);   // the step size     

      r = R_N/number_of_steps; // initial position      
      //Based on equation 3.1.17 in text
      u[0] = pow(r,L+1);       // initial wavefunction  
      u[1] = (L+1)*pow(r,L);   // initial wf derivative  

      // now we start solving the differential equations using the RK4 method 
      while (r <= R_N){
        derivatives(r, u, dudr, L, E);   // initial derivatives              
        runge_kutta_4(u, dudr, NUM_DIFF_EQ, r, step_size, 
            uout, L, E, derivatives); 
        for (i = 0; i < NUM_DIFF_EQ; i++) {
          u[i] = uout[i];  
        }
        r += step_size;
        if (e_index == e_index_0_1){
          output_radial(ofile_energy_0_1_files[L], r,u);
        }
        else if (e_index == e_index_3_0){
          output_radial(ofile_energy_3_0_files[L], r,u);
        }
      }

      r_matrix[L][e_index] = calculateRMatrix(u, R_N);
      s_matrix[L][e_index] = calculateSMatrix(u, R_N, r,L,E);
      phase_shifts[L][e_index] = calculatePhaseShift(s_matrix[L][e_index]);
      
//    std::cout << "For L = " << L << " and E = " << E  << ": " << std::endl
//              << "R Matrix: " << r_matrix[L][e_index]
//              << "S Matrix: " << s_matrix[L][e_index]
//              << "Phase Shift is: " << phase_shifts[L][e_index] << std::endl;
      delete [] u; delete [] dudr; delete [] uout; 
    }
  }

  output_phase_shifts(phase_shift_file , phase_shifts, E_values);
  phase_shift_file.close();
  for (int L = 0; L < NUM_L_VALUES;L++){
    ofile_energy_0_1_files[L].close();
    ofile_energy_3_0_files[L].close();
  }
  return 0;
}   //  End of main function 

//     Read in from screen the number of steps,
//     initial position and initial speed 
void initialize (int& number_of_steps)
{
  cout << "Number of steps: ";
  cin >> number_of_steps;
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


/*   This function upgrades a function y (input as a pointer)
     and returns the result uout, also as a pointer. Note that
     these variables are declared as arrays.  It also receives as
     input the starting value for the derivatives in the pointer
     dydx. It receives also the variable n which represents the 
     number of differential equations, the step size h and 
     the initial value of x. It receives also the name of the
     function *derivs where the given derivative is computed
 */
void runge_kutta_4(double *y, double *dydx, int n, double x, double h, 
    double *uout, int L, double E, void (*derivs)(double, double *, double *, int, double))
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
  //      now we upgrade y in the array uout  
  for (i = 0; i < n; i++){
    uout[i] = y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  }
  delete []dym;
  delete [] dyt;
  delete [] yt;
}       //  end of function Runge-kutta 4
