#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;
int n=51;
vector<double> u_new, u_old;
float error_tolerance = 1.e-3;
extern double RUN_NewtonStep();
ofstream out_file;
ofstream out_python;

int main(int argc, char *argv[])
{
	// Datenstrukturen (dynamnisch)
	u_new.resize(n);
    u_old.resize(n);
	out_file.open("out.txt");
    //out_file.setf(ios::scientific);
    out_file.precision(3);
	out_python.open("out.csv");
    // Anfangsbedingungen
    for(int x=0;x<n;x++)
    {
      u_old[x] = 0.25;
    }
    // Randbedingungen
    u_old[n-1] = 0.35; // Wasserstand flußabwärts [m]
    u_new[n-1] = 0.35; // Wasserstand flußabwärts [m]
	int i = 0;
	float error = 1.1*error_tolerance;
	while(error>error_tolerance)
	{
		error = RUN_NewtonStep();
		cout << error << endl;
		i++;
	}
	cout << i;	
    // Ausgabe der Ergebnisse
      // File 
	out_file << "Water depth (new):\t\t";
    for(int i=0;i<n;i++)
    {
      out_file << "\t" << u_new[i] << " ";
    }
    out_file << endl;
    out_python.close();
}

double RUN_NewtonStep()
{
	// Geometrie
	double x[n];
    for(int i=0;i<n;i++)
      x[i] = -100. + i*10.;
	double bottom_elevation[n];
    for(int i=0;i<n;i++)
      bottom_elevation[i] = 0.04 - i*0.004;
	// Parameter (konstant)
	double discharge = 0.07; // Volumenfließrate [m3/s]
    double gravity = 9.81; // [m/s2]
    double friction_law_exponent = 0.5; // Chezy, Manning-Strickler [-]
    double error_tolerance = 1e-3; // [m]
    double bed_slope = 0.0004; // [m/m]
    double bottom_width = 1.; // [m]
    double m = 1.; //
    double friction_coefficient = 10.; //
	// Parameter (variabel)
	double wetted_cross_section[n];
    double water_level_elevation[n];
    double flow_velocity[n];
    double Froude_number[n];
    double wetted_perimeter[n];
    double hydraulic_radius[n];
    double friction_slope[n];
	ofstream out_file2;
	out_file2.open("out2.txt");
    //out_file.setf(ios::scientific);
    out_file2.precision(3);
  //local variables
  double N,N1,N2,N3,D,D1,D2,D21,D22;
  double error = 0;
  {
    //start values
    for(int i=0;i<n;i++)
    {
      wetted_perimeter[i] = bottom_width + 2.*sqrt(1.+m*m)*u_old[i];
      wetted_cross_section[i] = (bottom_width + m*u_old[i])*u_old[i];
      hydraulic_radius[i] = wetted_cross_section[i] / wetted_perimeter[i];
      water_level_elevation[i] = bottom_elevation[i] + u_old[i];
      flow_velocity[i] = discharge/wetted_cross_section[i];
      Froude_number[i] = flow_velocity[i]/(sqrt(gravity*wetted_cross_section[i]\
                       /sqrt(bottom_width*bottom_width+4.*m*wetted_cross_section[i])));
      friction_slope[i] = pow(flow_velocity[i]/(friction_coefficient*pow(hydraulic_radius[i],friction_law_exponent)),2);
    }
	// Test output
    out_file2 << "Water depth (old):\t\t";
    for(int i=0;i<n;i++)
    {
      out_file2 << "\t" << u_old[i] << " ";
    }
    out_file2 << endl;
    out_file2 << "Wetted perimeter:\t\t";
    for(int i=0;i<n;i++)
    {
      out_file2 << "\t" << wetted_perimeter[i] << " ";
    }
    out_file2 << endl;	
    //Newton step
    for(int i=0;i<n-1;i++)
    {
      N1 = pow(discharge,2)/pow(wetted_cross_section[i+1],2) + gravity*u_old[i+1];
      N2 = pow(discharge,2)/pow(wetted_cross_section[i],2) + gravity*u_old[i];
      N3 = gravity*(bed_slope - (friction_slope[i+1]+friction_slope[i])/2.)*(x[i+1]-x[i]);
      N = N1 - N2 - N3;
      D1 = pow(discharge,2)/pow(wetted_cross_section[i],3) * (bottom_width+2.*m*u_old[i]) - gravity;
      D21 = friction_law_exponent*2.*(sqrt(1+m*m))/wetted_perimeter[i];
      D22 = (1.+friction_law_exponent)/wetted_cross_section[i] * (bottom_width+2.*m*u_old[i]);
      D2 = gravity*friction_slope[i]*(D21-D22)*(x[i+1]-x[i]);
      D = D1 + D2;
      u_new[i] = u_old[i] - N/D;
    }
    //calc Newton error
    for(int i=0;i<n-1;i++)
    {
      error += u_old[i] -u_new[i];
    }
    error = sqrt(error*error);
    //error = abs(error);
    //save Newton step
    for(int i=0;i<n-1;i++)
    {
      u_old[i] = u_new[i];
    }
	// für Python plot
    for(int i=0;i<n;i++)
    {
      out_python << x[i] << "," << u_new[i] << endl;
    }
  }
  return error;
}