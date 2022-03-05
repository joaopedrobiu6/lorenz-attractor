#include "ODE_analysis.h"

ODE_analysis::ODE_analysis(int dim, const std::vector<double> &info, const std::initializer_list<double> &X)
{
    obj = info;
    X0 = X;
    f = new std::function<double(ODEpoint)>[dim]();
};

ODE_analysis::ODE_analysis(int dim, const std::initializer_list<double> &X)
{
    Xvar var(X);
    X0 = var;
    f = new std::function<double(ODEpoint)>[dim]();
};

void ODE_analysis::SetFunction(int index, std::function<double(ODEpoint)> f_)
{
    f[index] = f_;
};

const std::vector<ODEpoint> &ODE_analysis::EulerSolver(double Time, double step)
{
    // Time = tempo final
    // n  = quantos pontos vamos retirar entre tempo0 e Time
    /* basicamente em cada n (instante tempo0 < t < Time) vamos calcular um os Xvar (variaveis dependentes)
    (os instantes sao ODEpoints - a cada um deles esta associado um Xvar) */

    int n = Time / step;

    // resizes
    resultado.resize(n);    // precisamos que haja tantos ODEpoints quanto instantes de tempo!
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes (dimensão dos Xvar!!!)

    for (int i = 0; i < n; i++)
        resultado[i].X().resize(n1); // redimensionar o Xvar de cada um de n ODEpoints

    // condições iniciais - o primeiro ODEpoint que corresponde às condições do construtor!!!!
    resultado[0].T() = 0;

    for (int i = 0; i < n1; i++)
        resultado[0].X()[i] = X0.X()[i]; // colocar os dados do construtor no primeiro ODEpoint!!!

    // calculo das soluções

    for (int i = 0; i < n - 1; i++)
    {
        resultado[i + 1].T() = step + resultado[i].T();
        for (int j = 0; j < n1; j++)
        {
            resultado[i + 1].X()[j] = resultado[i].X()[j] + step * f[j](resultado[i]);
        }
    }
    return resultado;
};

const std::vector<ODEpoint> &ODE_analysis::RungeKutta4(double Time, double step)
{
    std::vector<double> K1, K2, K3, K4;

    // resizes
    int n = Time / step;
    int n1 = X0.X().size();

    resultado.resize(n);
    K1.resize(n1);
    K2.resize(n1);
    K3.resize(n1);
    K4.resize(n1);

    for (int i = 0; i < n; i++)
    {
        resultado[i].X().resize(n1);
    }

    // condições iniciais no tempo
    resultado[0].T() = 0;

    for (int i = 0; i < n1; i++)
    {
        resultado[0].X()[i] = X0.X()[i];
    }

    ODEpoint p0;
    Xvar xvar_temp;
    double tempo;
    p0.X().resize(n1);
    xvar_temp.X().resize(n1);

    for (int i = 0; i < n - 1; i++)
    {
        // Calcular o tempo em i+1
        resultado[i + 1].T() = resultado[i].T() + step;

        // Cálculo dos K1
        for (int j = 0; j < n1; j++)
        {
            K1[j] = f[j](resultado[i]);
        }
        // K1[1] = f[1](resultado[i]);
        for (int j = 0; j < n1; j++)
        {
            xvar_temp.X()[j] = resultado[i].X()[j] + (0.5 * step * K1[j]);
        }
        // xvar_temp.X()[1] = resultado[i].X()[1] + (0.5 * step * K1[1]);
        tempo = resultado[i].T() + (step * 0.5);

        p0.SetODEpoint(tempo, xvar_temp);

        // Cálculo dos K2
        for (int j = 0; j < n1; j++)
            K2[j] = f[j](p0);
        // K2[1] = f[1](p0);

        for (int j = 0; j < n1; j++)
            xvar_temp.X()[j] = resultado[i].X()[j] + (0.5 * step * K2[j]);
        // xvar_temp.X()[1] = resultado[i].X()[1] + (0.5 * step * K2[1]);

        p0.SetODEpoint(tempo, xvar_temp);

        // Cálculo dos K3
        for (int j = 0; j < n1; j++)
            K3[j] = f[j](p0);
        // K3[1] = f[1](p0);

        // Cálculo dos K4
        for (int j = 0; j < n1; j++)
            xvar_temp.X()[j] = resultado[i].X()[j] + (step * K3[j]);
        // xvar_temp.X()[1] = resultado[i].X()[1] + (step * K3[1]);
        tempo = resultado[i].T() + (step);

        p0.SetODEpoint(tempo, xvar_temp);

        for (int j = 0; j < n1; j++)
            K4[j] = f[j](p0);
        // K4[1] = f[1](p0);

        // Cálculo das soluções
        for (int j = 0; j < n1; j++)
            resultado[i + 1].X()[j] = resultado[i].X()[j] + step * (K1[j] + 2 * K2[j] + 2 * K3[j] + K4[j]) / 6;
        // resultado[i + 1].X()[1] = resultado[i].X()[1] + step * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]) / 6;

        K1.clear();
        K2.clear();
        K3.clear();
        K4.clear();
    }

    return resultado;
};

const std::vector<ODEpoint> &ODE_analysis::LeapFrogImprovedSolver(double Time, double step)
{
    // std::vector<ODEpoint> resultado;
    int n = Time / step;

    // resizes
    resultado.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        resultado[i].X().resize(n1);

    // condições iniciais
    resultado[0].T() = 0;

    for (int i = 0; i < n1; i++)
        resultado[0].X()[i] = X0.X()[i];

    // primeira iteração
    resultado[1].T() = step + resultado[0].T();
    resultado[1].X()[0] = resultado[0].X()[0] + step * f[0](resultado[0]);
    resultado[1].X()[1] = resultado[0].X()[1] + step * f[1](resultado[0]);

    for (int i = 1; i < n - 1; i++)
    {
        resultado[i + 1].T() = step + resultado[i].T();
        resultado[i + 1].X()[0] = resultado[i].X()[0] + step * f[0](resultado[i]) + 0.5 * step * step * f[1](resultado[i]);
        resultado[i + 1].X()[1] = resultado[i].X()[1] + 0.5 * step * (f[1](resultado[i]) + f[1](resultado[i + 1]));
    }

    return resultado;
};

