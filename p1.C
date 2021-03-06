// Calculate the Voltage for a 2D system with wire at fixed voltage immediately above
// a 3-sided box at 0V

#include "TGraph2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBox.h"
#include "TApplication.h"

#include <getopt.h>
#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;

const int L = 100;   // Length of any side
const double V0 = 100;  // Voltage at top of box
const int dist = 45; // distance from plate to wall, should be < Npts / 2 for best results 
const int maxgraphlines = 200;  // Max lines to draw in each direction

// Generic code to do one iteration of finite difference method
// Jacobi Method
double iterateJ(vector<vector<double>> &V) {
  auto Vtmp = V;
  double dVmax = 1e-50;
  int nx = V.size();
  int ny = V[0].size();
  for(int i = 1; i < nx - 1; i++) {
    for(int j = 1; j < ny - 1; j++) {
      double Vnew = 0.25 * (Vtmp[i + 1][j] + Vtmp[i - 1][j] +
			    Vtmp[i][j + 1] + Vtmp[i][j - 1]);
      double dV = fabs(Vnew - V[i][j]);
      dVmax = std::max(dVmax, dV);  // Keep track of max change in this sweep
      V[i][j] = Vnew;
    }
  }
  return dVmax;
}


// Gauss-Seidel Method
double iterateGS(vector<vector<double>> &V, int Npts) {
  double dVmax = 1e-50;
  int nx = V.size();
  int ny = V[0].size();
  for(int i = 1; i < nx - 1; i++) {
    for(int j = 1; j < ny - 1; j++) {
      if(!((j == dist || j == ny - 1 - dist) && (i < L/2+Npts/2 || i > L/2-Npts/2))) {  // Preserve potential at plate
	double Vnew = 0.25 * (V[i + 1][j] + V[i - 1][j] +
			      V[i][j + 1] + V[i][j - 1]);
	double dV = fabs(Vnew - V[i][j]);
	dVmax=std::max(dVmax, dV);  // Keep track of max change in this sweep
	V[i][j] = Vnew;
      }
    }
  }
  return dVmax;
}

// Fill a TGraph2D object from a vector of voltages
// delta: grid spacing
// The optional range parameter defines the subregion to plot
void fillGraph(TGraph2D* tg, const vector<vector<double>> &V, double delta,
	       TBox *range = 0) {
  int nx = V.size();
  int ny = V[0].size();
  tg->Clear();               // Reset the graph
  for(int i = 0; i < nx; i++) {
    double x = i * delta;
    for(int j = 0; j < ny; j++) {
      double y = j * delta;
      if(range && range->IsInside(x, y))
	tg->SetPoint(tg->GetN(), x, y, V[i][j]);
    }
  }
}

// Define box 0 < x < L, 0 < y < L
// Potential on top edge at y = L
// eps: convergence criteria (max size of change at any grid point in an iteration)
// maxIter: max iterations in case of non-converence
// Npts : smoothness parameter, number of grid points in x,y
// pass a tcanvas for an animated solution, with specified max rate of frames/second
TGraph2D* LaplaceLine(int maxIter = 100, double eps = 0.001, int Npts = 100,
		      TCanvas *tc = 0, int rate = 10) {
   // create N x N vector, init to 
  vector<vector<double>> V(L, vector<double> (L, 0));
  double delta = 1;  // Grid spacing

  
  for(int i = L/2-Npts/2; i < L/2+Npts/2; i++) {
    V[i][L - 1 - dist] = -V0;  // Set voltage at plates
    V[i][dist]        =  V0;
  }
  
  

  int msec = 1000 / rate;   // Milliseconds sleep between frames
  TBox *plotRange = new TBox(0, 0, 1.1 * L, 1.1 * L);

  TGraph2D* tgV = new TGraph2D();  // Graph to store result               
  tgV->SetLineColor(kRed);
  tgV->SetNpx(std::min(maxgraphlines, L));
  tgV->SetNpy(std::min(maxgraphlines, L)); 
  tgV->SetTitle("Voltage;x;y;V");
  
  double dV;
  int niter = 0;

  do {
    //dV = iterateJ(V);   // iterate using Jacobi method
    dV = iterateGS(V, Npts);   // iterate using Gauss-Seidel method
    ++niter;
    if (tc) {
      tc->cd();
      fillGraph(tgV, V, delta, plotRange);
      //tgV->Draw("surf");
      //tc->Update();
      //gSystem->Sleep(msec);
    }
  } while(dV > eps && niter < maxIter);
  

  cout << "Ended calculation with " << niter << " iterations, dVmax = " << dV << endl;

  fillGraph(tgV, V, delta, plotRange);
  tgV->Draw("surf");
  //tc->Update();
  return tgV;
}

void usage(char *prog) {
  std::cerr << "Usage: " << prog << " <option(s)> SOURCES"
	    << "Options:\n"
	    << "\t-h\t\tShow this help message\n"
	    << "\t-a\t\tDisplay animation of solution\n"
    	    << "\t-I\t\t(max) Number of iterations [100]\n"
	    << "\t-e\t\tconvergence criteria [0.001]\n"
    	    << "\t-N\t\tNumber of points in x,y [100]\n"
	    << "\t-R\t\tmax frames/second with animation [10]"
	    << std::endl;
  exit(0);
}


int main(int argc, char *argv[]) {
  TApplication theApp("App", &argc, argv, NULL, -1);
  // -1 disables ROOT arg processing

  // Defaults for LaplaceLine
  int maxIter = 100;
  double eps  = 0.001;
  int Npts    = 20;
  TCanvas *tc = 0;
  int rate    = 10;
  
  int opt;
  while((opt = getopt(argc, argv, "haI:e:N:r:")) != -1) {
    switch(opt) {
    case 'h':
      usage(argv[0]);
      break;
    case 'a':
      tc = new TCanvas();
      break;
    case 'I':
      maxIter = atoi(optarg);
      break;
     case 'e':
      eps = atof(optarg);
      break; 
    case 'N':
      Npts = atoi(optarg);
      break;
    case 'r':
      rate = atoi(optarg);
      break;
    }
  }
 
  auto tg = LaplaceLine(maxIter, eps, Npts, tc, rate);

  // Display final result
  if(!tc) tc = new TCanvas();
  tg->Draw("surf");  // explore other drawing options!
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30, ".q");
  // set up a failsafe timer to end the program  
  theApp.Run();
}
