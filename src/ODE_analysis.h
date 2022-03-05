#ifndef __pendulum_analysis__
#define __pendulum_analysis__

// ROOT INCLUDES
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include <iostream>
#include <vector>
#include <functional>
#include <map>
#include <cmath>

#include "ODEpoint.h"

class ODE_analysis
{
public:
    ODE_analysis() = default;
    ODE_analysis(int dim, const std::initializer_list<double> &);                              // pendulum_parametric(10, {80, 0})
    ODE_analysis(int dim, const std::vector<double> &, const std::initializer_list<double> &); // pendulum_parametric(10, {80, 0})
    // length in meters, angles in degrees, velocity degrees/sec
    ~ODE_analysis() = default;

    void SetFunction(int, std::function<double(ODEpoint)>);

    // Solvers
    const std::vector<ODEpoint> &EulerSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &RungeKutta4(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &LeapFrogImprovedSolver(double Time = 0, double step = 1E-3); // método usado para confirmar valores do RK4


private:
    std::vector<double> obj; // for
    Xvar X0;                 // initial conditions: angle (rad), angular velocity (rad/s)
    std::vector<ODEpoint> resultado;
    std::function<double(ODEpoint)> *f; // Recebe uma odepoint e retorna um double
};

#endif
