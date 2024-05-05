#include<iostream>
#include<fstream>
#include<math.h>
#include<sstream>
#include <iomanip>
using namespace std;

//Spherical case: the mass should be 25.9148, and the energy should be -0.692229.

/*

This code solves the solution of eq 1 and eq 2 with given phi_s.
Psi, chi, and the minimum eigenvalue gamma are calculated.

1) \nabla^2 \phi = \chi^2
2) \nabla^2 \chi = 2 * ( \phi + \phi_s - \gamma ) + \chi

\phi : gravitational potential
\chi : scalar field wave function
\phi_s : star potential (external field)
\gamma : energy eigenvalue

The setup is
 a) cylindrical symmetry and parity symmetry (z -> -z) are assumed.
 b) The boundary conditions are
    * phi(rmax, z) = phi(r, zmax) = 0
    * chi(rmax, z) = chi(r, zmax) = 0
    * chi(0,0) is fixed to some nonzero value
*/

//Reminder: Maximum lam_power = 6.0

//std::ostringstream galaxy_name <<  "UGC01281";

//galaxy_name ;

const double lam_power = 4.7;

const double lambda = 1.0*pow(10,-lam_power);

const int num = 400 + 1; // number of points in each dimension
const double rmax = 10./lambda; // maximum value of r = sqrt(x^2+y^2)
const double zmax = 10./lambda; // maximum value of z
const double dr = rmax/(num-1); // spacing of r
const double dz = zmax/(num-1); // spacing of z

inline double r(int i){ return rmax * i / (num-1); }
inline double z(int i){ return zmax * i / (num-1); }


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// field defined in cylindrical coordinate


class field{
  public:
	double val[num][num];
	// appling Laplacian
	double laplacian(int i, int j) {
		double result = 0.;

		// d^2phi / drho^2 + 1/rho dphi/drho
		if(i==0) {
			result += (val[i+1][j] -   val[i][j]              )/dr/dr * 4.;
		} else if (i==(num-1)) {
			result += (            -1.*val[i][j] + val[i-1][j])/dr/dr +(val[i][j] - val[i-1][j])/dr/r(i);
		} else {
			//result += (val[i+1][j] -2.*val[i][j] + val[i-1][j])/dr/dr +(val[i+1][j] - val[i][j])/dr/r(i); // old derivative
            result += (val[i+1][j] -2.*val[i][j] + val[i-1][j])/dr/dr +(val[i+1][j] - val[i-1][j])/2./dr/r(i); // new derivative
		}
		// d^2phi / dz^2
		if(j==0) {
			result += (val[i][j+1] -   val[i][j]              )/dz/dz * 2.;
		} else if (j==(num-1)) {
			result += (            -1.*val[i][j] + val[i][j-1])/dz/dz;
		} else {
			result += (val[i][j+1] -2.*val[i][j] + val[i][j-1])/dz/dz;
		}
	
		return result;
	}
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*

Solver for Schrodinger-Poisson equation.

1) Poission equation is solved by SOR method.

2) Schrodinger equation is solved by calculating "time" evolution of
   (d/dt) chi = -\nabla^2 chi + 2(\phi + \phi_s)\chi
In large t limit, chi satisfy Schrodinger eq.

chi(0,0) is always set to chi0.

*/
class SchrodingerPoissonSolver{
  public:
	field chi, phi, phistar;
	double chi0; // chi at the origin.
	SchrodingerPoissonSolver(){
		chi0 = 1.;

		// initial condition of phi and chi
		for (int i=0; i< num; i++){ for (int j=0; j< num; j++){
			if( (i==(num-1))||(j==(num-1)) ){
				chi.val[i][j] = 0.;
			} else {
				chi.val[i][j] = exp( -0.2*(r(i)*r(i)+z(j)*z(j)));
			}
			phi.val[i][j]=0.;
		}}
		for (int i=0; i< num; i++){ for (int j=0; j< num; j++){
			phistar.val[i][j]=0.;
		}}
	}

	// SOR method
	double improvePhi(double dt){
		double maxdiff=0.;
		double phinew[num][num];
		for (int i=0; i<(num-1); i++){ for (int j=0; j<(num-1); j++){
			double phinew;
			phinew = phi.val[i][j] + dt * (phi.laplacian(i,j) - pow(chi.val[i][j],2));
			maxdiff = max(maxdiff, fabs(phinew-phi.val[i][j]));
			phi.val[i][j] = phinew;
		}}
		return maxdiff;
	}

	// "time" evolution for (d/dt) chi = -H chi
	// wave function is always renormalized to satisfy chi(0) = chi0
	// in large t limit, chi will be the lowest energy state
	double improveChi(double dt){
		double chinew[num][num];
		// "time" evolution
		for (int i=0; i<(num-1); i++){ for (int j=0; j<(num-1); j++){
			chinew[i][j] = chi.val[i][j] + dt*(chi.laplacian(i,j) - 2.*(phi.val[i][j] + phistar.val[i][j])*chi.val[i][j]);
		}}

		// renormalization to satisfy the same condition at the origin
		// return the maximum value of the difference between old phi and new phi
		double maxdiff = 0.;
		double chi0new = chinew[0][0];
		for (int i=0; i<(num-1); i++){ for (int j=0; j<(num-1); j++){
			chinew[i][j] = chinew[i][j] / chi0new * chi0;
			maxdiff = max(maxdiff, fabs(chinew[i][j]-chi.val[i][j]));
			chi.val[i][j] = chinew[i][j];
		}}

		return maxdiff;
	}

	// evaluate (H psi)(x) / psi(x) at x = (r[i], z[j])
	double energy(int i, int j){
		return ( -0.5 * chi.laplacian(i,j) + (phi.val[i][j] + phistar.val[i][j]) * chi.val[i][j] ) / chi.val[i][j];
	}

	// \int d^3x chi^2(x)
	double mass(){
		double result = 0.;
		for (int i=0; i<(num-1)-2; i++){ for (int j=0; j<(num-1)-2; j++){
			//double rho = pow(chi.val[i][j], 2.);
            double chiavg= (chi.val[i][j] + chi.val[i+1][j] + chi.val[i][j+1] + chi.val[i+1][j+1]) / 4.0;
            double rho = pow(chiavg, 2.);
			result += 2. * dz * M_PI*( r(i+1)*r(i+1)-r(i)*r(i) )*rho;
		}}
		return result;
	}

	// update the boundary conditions for phi on the edge (phi(r) = -M/4/Pi/r) 
	void updatePhiBoundary(){
		double totalmass = mass();
		for(int i=0; i<(num-1); i++){
			phi.val[i][num-1] = -1.*totalmass/4./M_PI/sqrt( r(i)*r(i) + z(num-1)*z(num-1) );
			phi.val[num-1][i] = -1.*totalmass/4./M_PI/sqrt( r(num-1)*r(num-1) + z(i)*z(i) );
		}
		phi.val[num-1][num-1] = -1.*totalmass/4./M_PI/sqrt( r(num-1)*r(num-1) + z(num-1)*z(num-1) );
	}

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[] ){
	SchrodingerPoissonSolver spsolver;
    
	// set boundary condition at the origin
    spsolver.chi0 = 1.*(pow(lambda,2));
    
    //std::ostringstream galaxy_name;
    
    std::stringstream galaxy_name(argv[1]);
        
    cout << "Galaxy is " << galaxy_name.str() << endl;
    cout << "chi(0) = " << pow(lambda,2) << "\n";
    
    std::ostringstream infilename;
    infilename << galaxy_name.str() << "_potential_lampower=" << std::setprecision(5) << lam_power << "_fine.txt";
    //infilename << galaxy_name.str() << "_potential_lampower=5_fine.txt";
    
    std::ifstream in(infilename.str());
    std::string line;
    
    if(!in) cout << "=============== No File! =================\n";
    
    double gal_data[401][401] = {0};
    int i = 0, k = 0;
    
    while (std::getline(in, line))
    {
        double value;
        k = 0;
        std::stringstream ss(line);
        
        while (ss >> value)
        {
            gal_data[i][k] = value;
            //cout << gal_data[i][k] << "\n";
            ++k;
        }
        ++i;
    }

	// set baryonic potential
	for (int i=0; i< num; i++){ for (int j=0; j< num; j++){     
		//spsolver.phistar.val[i][j] = 0.;
		//spsolver.phistar.val[i][j] = -1.0*exp( - 0.1*r(i)*r(i)- 1.*z(j)*z(j));
        spsolver.phistar.val[i][j] = gal_data[i][j];
	}}

	// calculate phi & chi to satisfy Schrodinger-Poisson eqs.
	cout << "nsteps\tdeltaPhi\tdeltaChi\tenergy" << endl;
	double dt = dr*dr/4.;
	double diffPhi;
	double diffChi;
	for(int n=0; n<300000; n++){
		diffPhi = spsolver.improvePhi(dt*1.5);
		diffChi = spsolver.improveChi(dt*0.8);
		spsolver.updatePhiBoundary();
		if(((n+1) % 10000)==0) {
			cout << n+1 << "\t";
			cout << diffPhi << "\t";
			cout << diffChi << "\t";
			cout << spsolver.energy(0,0) << endl;
		}
	}
	cout << endl;
	cout << "total mass is \t" << spsolver.mass() << endl;
	cout << "energy eigenvalue is \t" << spsolver.energy(0,0) << endl;
    cout << "{-log(lam), M, gamma} = {" << lam_power << ", " << spsolver.mass() << ", " << spsolver.energy(0,0) << "}" << endl;

	// save the result of phi & chi
	ofstream outputfile("solution.dat");
	outputfile << "# r\tz\tphi\tchi" << endl;
	for (int i=0; i< num; i++){ for (int j=0; j< num; j++){
		outputfile << r(i) << "\t";
		outputfile << z(j) << "\t";
		outputfile << spsolver.phi.val[i][j] << "\t";
		outputfile << spsolver.chi.val[i][j] << "\t";
		outputfile << endl;
	}
		outputfile << endl;
	}

	return 0;
}

