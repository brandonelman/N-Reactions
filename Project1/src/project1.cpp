/*    This program solves Newton's equation for a block
      sliding on a horizontal frictionless surface. The block
      is tied  to a wall with a spring, and Newton's equation
      takes the form
      m d^2x/dt^2 =-kx
      with k the spring tension and m the mass of the block.
      The angular frequency is omega^2 = k/m and we set it equal
      1 in this example program. 

      Newton's equation is rewritten as two coupled differential
      equations, one for the position x  and one for the velocity v
      dx/dt = v    and
      dv/dt = -x   when we set k/m=1

      We use therefore a two-dimensional array to represent x and v
      as functions of t
      y[0] == x
      y[1] == v
      dy[0]/dt = v
      dy[1]/dt = -x

      The derivatives are calculated by the user defined function 
      derivatives.

      The user has to specify the initial velocity (usually v_0=0)
      the number of steps and the initial position. In the programme
      below we fix the time interval [a,b] to [0,2*pi].

 */ 
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include "lib.h"
using namespace  std;
// output file as global variable
ofstream ofile;
// function declarations
void derivatives(double, double*,  double*,
                    int, double);
void initialise ( double&, double&, int&);
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

int main(int argc, char* argv[])
{
  //  declarations of variables
  double *u, *dudr, *yout, r, step_size, r_max, E0;
  double initial_u, initial_u_prime;
  int i, number_of_steps, num_diff_eq;
  char *outfilename;

  const double R_N = 100;      //Range of potential
  const double   A = R_N; //Matching point outside potential

  int L = 0; //angular momentum quantum number
  double E = 0.1; //Energy (MeV)
  
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
  initialise (initial_u, initial_u_prime, number_of_steps);
  //  setting initial values, step size and max position r_max  
  step_size = A / ( (double) number_of_steps);   // the step size     
  r_max = step_size*number_of_steps;       // the final time    
  u[0] = initial_u;                       // initial wavefunction  
  u[1] = initial_u_prime;                 // initial wf derivative  
  r=0.1;                                  // initial position      
  // now we start solving the differential equations using the RK4 method 
  while (r <= A){
    derivatives(r, u, dudr, L, E);   // initial derivatives              
    runge_kutta_4(u, dudr, num_diff_eq, r, step_size, 
                  yout, L, E, derivatives); 
    for (i = 0; i < num_diff_eq; i++) {
      u[i] = yout[i];  
    }
    r += step_size;
    output(r, u);   // write to file 
  }
  delete [] u; delete [] dudr; delete [] yout;
  ofile.close();  // close output file
  return 0;
}   //  End of main function 

//     Read in from screen the number of steps,
//     initial position and initial speed 
void initialise (double& initial_u, double& initial_u_prime, int& number_of_steps)
{
  cout << "Wavefunction at R = 0: ";
  cin >> initial_u;
  cout << "Derivative of WF at R = 0: ";
  cin >> initial_u_prime;
  cout << "Number of steps = ";
  cin >> number_of_steps;
}  // end of function initialise  

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
