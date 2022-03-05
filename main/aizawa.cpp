#include <iostream>
#include <vector>
#include <fstream>

#include "TGraph2D.h"
#include "ODE_analysis.h"

void DrawImage(double tmax, double *xf, double *yf, double *zf);
void OpenApp(double tmax, double *xf, double *yf, double *zf);

int main()
{
    double a = 0.95, b = 0.7, c = 0.6, d = 3.5, e = 0.25, f = 0.1;
    double x0 = 0., y0 = 1., z0 = 1.05;

    double tmax = 1000.;

    ODE_analysis lorenz(3, {x0, y0, z0});

    lorenz.SetFunction(0, [&](ODEpoint p)
                       { return (p.X()[2] - b) * p.X()[0] - d * p.X()[1]; });

    lorenz.SetFunction(1, [&](ODEpoint p)
                       { return d * p.X()[0] + (p.X()[2] - b) * p.X()[1]; });

    lorenz.SetFunction(2, [&](ODEpoint p)
                       { return c + a * p.X()[2] - (p.X()[2] * p.X()[2] * p.X()[2]) / 3 - (p.X()[0] * p.X()[0] + p.X()[1] * p.X()[1]) * (1 + e * p.X()[2]) + f * p.X()[2] * (p.X()[0] * p.X()[0] * p.X()[0]); });

    auto result = lorenz.RungeKutta4(tmax, 1e-3);

    double *tempo = new double[result.size()];
    double *xf = new double[result.size()];
    double *yf = new double[result.size()];
    double *zf = new double[result.size()];

    for (int i = 0; i < result.size(); i++)
    {
        tempo[i] = result[i].T();
        xf[i] = result[i].X()[0];
        yf[i] = result[i].X()[1];
        zf[i] = result[i].X()[2];
    }

    DrawImage(tmax, xf, yf, zf);

    OpenApp(tmax, xf, yf, zf);

    delete[] tempo;
    delete[] xf;
    delete[] yf;
    delete[] zf;

    return 0;
}

void DrawImage(double tmax, double *xf, double *yf, double *zf)
{
    TCanvas *c = new TCanvas("canvas", "Pendulum", 0, 0, 1280, 720);
    TGraph2D *gr = new TGraph2D(tmax / 1e-2, xf, yf, zf);
    gr->SetTitle("Aizawa Attractor");

    gStyle->SetPalette(1);
    gr->Draw("pcol LC");

    c->SaveAs("aizawa.png");

    delete gr;
    delete c;
};

void OpenApp(double tmax, double *xf, double *yf, double *zf)
{
    TApplication *app = new TApplication("app", nullptr, nullptr);
    TCanvas *c = new TCanvas("canvas", "attractors", 0, 0, 1280, 720);

    TRootCanvas *r = (TRootCanvas *)c->GetCanvasImp();
    r->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TGraph2D *gr = new TGraph2D(tmax / 1e-2, xf, yf, zf);
    gr->SetTitle("Aizawa Attractor");

    gStyle->SetPalette(1);
    gr->Draw("pcol LC");

    c->Update();
    app->Run();
    delete gr;
    delete c;
};