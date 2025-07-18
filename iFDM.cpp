#include <vector>
//#include <cmath>
#include <fstream>
using namespace std;
//
extern void Gauss(double *matrix, double *vecb, double *vecx, int g);
void AssembleEquationSystem();
//
int nj = 101;
int n=nj;
double* matrix;
double* vecb;
double* vecx;
double Ne;
vector<double> u_new, u_old;
std::ofstream out_file;
std::ofstream out_file_solver;
  
int main(int argc, char *argv[])
{
  //1-Definitionen
  //1.1 Diskretisierung
  double dt,dx,x;
  int nt = 100;
  dx = 1./(nj-1);
  //1.2 Lösungsvektoren
  u_new.resize(nj);
  u_old.resize(nj); 
  matrix = new double[nj*nj];
  vecb = new double[nj];
  vecx = new double[nj];
  //1.3 Kennzahlen, Parameter
  double Ne = 0.5;
  double alpha=1.;
  dt = Ne * dx*dx / alpha ; // Zeitschrittsteuerung
  //Ne = alpha * dt / (dx*dx);
  //1.4 Ausgabe
  out_file.open("out.csv");
  out_file_solver.open("matrix.txt");
  //2-Anfangsbedingungen
  for(int ix=0;ix<nj;ix++)
  {
    u_old[ix] = 0.;
  }
  //3-Randbedingungen
  double u_bc_l = 3.;
  double u_bc_r = -1.;
  u_new[0] = u_old[0] = u_bc_l;
  u_new[nj-1] = u_old[nj-1] = u_bc_r;
  //4-Diskretisierung
  //dt = 1.;
  dx = 1./(nj-1);
  //5-Formel
  //dt = 100. * dx*dx / alpha; // Zeitschrittsteuerung
  dt = 0.001;
  Ne = alpha * dt / (dx*dx);  
  //------------------------------------------------------
  //6-Berechnung
  int i,j;
  for(int t=0;t<nt;t++)
  {
//	AssembleEquationSystem();
//
  // Matrix entries
  for(i=0;i<n;i++)
  {
    vecb[i] = u_old[i];
    for(j=0;j<n;j++)
    {
      matrix[i*n+j] = 0.0;
      if(i==j)
        matrix[i*n+j] = 1. + 2.*Ne;
      else if(abs((i-j))==1)
        matrix[i*n+j] = - Ne;
    }
  }
  // Treat boundary conditions
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      if(i==0||i==n-1)
        matrix[i*n+j] = 0.0;
    }
  for(i=0;i<n;i++)
  {
    if(i!=0&&i!=n-1)
      continue;
    for(j=0;j<n;j++)
    {
       if(i==j)
         matrix[i*n+j] = 1.0;
       else
         matrix[i*n+j] = 0.0;
    }
  }
//
    Gauss(matrix,vecb,vecx,nj);	
	for(int i=0;i<nj;i++)
    {
      u_new[i] = vecx[i];
    }
    for(int ix=0;ix<nj;ix++)
    {
      x = dx*ix;
      out_file << x << "," << u_new[ix] << endl;
    }	
    // Daten speichern
    for(int i=1;i<nj-1;i++)
    {
      u_old[i] = u_new[i];
    }	
  }  
}

void AssembleEquationSystem()
{
  int i,j;
  // Matrix entries
  for(i=0;i<n;i++)
  {
    vecb[i] = u_old[i];
    for(j=0;j<n;j++)
    {
      matrix[i*n+j] = 0.0;
      if(i==j)
        matrix[i*n+j] = 1. + 2.*Ne;
      else if(abs((i-j))==1)
        matrix[i*n+j] = - Ne;
    }
  }
  // Treat boundary conditions
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      if(i==0||i==n-1)
        matrix[i*n+j] = 0.0;
    }
  for(i=0;i<n;i++)
  {
    if(i!=0&&i!=n-1)
      continue;
    for(j=0;j<n;j++)
    {
       if(i==j)
         matrix[i*n+j] = 1.0;
       else
         matrix[i*n+j] = 0.0;
    }
  }
  // Matrix output
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      out_file_solver << matrix[i*n+j] << " ";
    }
    out_file_solver << "b:" << vecb[i] << " ";
    out_file_solver << endl;
  }
}

