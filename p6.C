#include "TGraph2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBox.h"
#include "TApplication.h"
#include "TMath.h"

#include <getopt.h>
#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;


const int L  = 200; // Radius of grounded container
const double V0 = 100; // Voltage of outer cylinder
const int maxgraphlines = 200; // Max lines to draw in each direction
const int a = 40;
const int b = 70;
const int crack_size = 3;


// Generic code to do one iteration of finite difference method
// Jacobi Method
double iterateJ(vector<vector<double>> &V, vector<vector<double>> &rho,
		double delta) {
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
double iterateGS(vector<vector<double>> &V, vector<vector<double>> &rho,
		 double delta) {
  double dVmax = 1e-50;
  int nx = V.size();
  int ny = V[0].size();
  for(int i = 1; i < nx - 1; i++) {
    for(int j = 1; j < ny - 1; j++) {
      if((i != a && i != b) ||
	 (j >= ny / 2 - crack_size && j <= ny / 2 + crack_size)) {
	double Vnew = 0.5 * (V[i + 1][j] * (1 + 1 / (2 * i)) / (delta * delta) +
			     V[i - 1][j] * (1 - 1 / (2 * i)) / (delta * delta) +
			     (V[i][j + 1] + V[i][j - 1])
			     / (i * i * pow(delta, 4))) * delta * delta /
	  (1 + 1 / (i * i * delta * delta));
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
      double y = j * delta * 2 * TMath::Pi() / ny;
      if(range && range->IsInside(x, y)) {
	tg->SetPoint(tg->GetN(), x * TMath::Cos(y), x * TMath::Sin(y), V[i][j]);
      }
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
  vector<vector<double>> V(Npts, vector<double> (Npts, 0));
  vector<vector<double>> rho(Npts, vector<double> (Npts, 0));
  double delta =1;  // Grid spacing
  
  for(int i = 0; i < Npts; i++) {
    V[a][i] = -V0; // Inner cylinder
    V[b][i] =  V0; // Outer cylinder
    V[Npts - 1][i] = 0; // Grounded container
  }
  
  for(int j = Npts / 2 - crack_size; j < Npts / 2 + crack_size; j++) {
    V[a][j] = 0;
    V[b][j] = 0;
  }
  

  int msec = 1000 / rate;   // Milliseconds sleep between frames
  TBox *plotRange = new TBox(0, 0, 1.1 * L, 1.1 * L);

  TGraph2D* tgV = new TGraph2D();  // Graph to store result                      
  tgV->SetLineColor(kRed);
  tgV->SetNpx(std::min(maxgraphlines, Npts));
  tgV->SetNpy(std::min(maxgraphlines, Npts)); 
  tgV->SetTitle("Voltage;x;y;V");
  
  double dV;
  int niter = 0;
  
  do {
    //dV = iterateJ(V);   // iterate using Jacobi method
    dV=iterateGS(V, rho, delta);   // iterate using Gauss-Seidel method
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
  for(int i = 0; i < Npts; i++) {
    V[i][Npts - 1] = V[i][Npts - 2]; // Hacky way to associate 0 and 2pi in V
    V[i][0] = V[i][1];
  }
  
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
  int Npts    = 200;
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
  TCanvas *tc2 = new TCanvas();
  TCanvas *tc3 = new TCanvas();
  tc->cd();
  tg->Draw("surf");  // explore other drawing options!
  tc2->cd();
  tg->Draw("CONT");
  tc3->cd();
  TGraph2D *tg2 = ((TGraph2D*)tg->Clone());
  tg2->SetNpx(100);
  tg2->SetNpy(100);
  tg2->Draw("ARR COLZ");

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30, ".q");
  // set up a failsafe timer to end the program  
  theApp.Run();
}
