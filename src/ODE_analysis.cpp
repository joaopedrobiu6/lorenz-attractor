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

void ODE_analysis::ODE_Draw(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<ODEpoint> resultado)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o gráfico da soluções da ODE..." << std::endl;
    // FAZER GRÁFICO DA EDO
    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 3200, 1800);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    int n = resultado.size();

    double *tempo = new double[n];
    double *posicao = new double[n];

    for (int i = 0; i < n; i++)
    {
        tempo[i] = resultado[i].T();
        posicao[i] = resultado[i].X()[0];
    }

    TGraph *gr = new TGraph(n, tempo, posicao);
    gStyle->SetTitleFontSize(0.05);
    gr->SetLineColor(color);
    gr->SetLineWidth(10);

    gr->SetTitle("Pendulo");
    gr->GetXaxis()->SetTitle(xaxis);
    gr->GetYaxis()->SetTitle(yaxis);
    gr->Draw("AC");

    c->SaveAs(filename);

    delete[] tempo;
    delete[] posicao;

    delete c;
};

std::vector<double> ODE_analysis::Moving_Average(std::vector<ODEpoint> vetor, double time_step, double tw)
{
    double sum_sq;
    std::vector<double> amp_mean;
    amp_mean.resize(vetor.size() - tw);
    int loop;

    for (int i = 0; i < vetor.size() - tw; i++)
    {
        amp_mean[i] = 0;
    }

    loop = 0;
    while (loop < tw)
    {
        sum_sq += vetor[loop].X()[0] * vetor[loop].X()[0];
        loop++;
    }
    double media = sqrt(sum_sq / tw);
    amp_mean[0] = media;

    for (int i = 1; i < vetor.size() - tw; i++)
    {
        sum_sq += vetor[i + tw].X()[0] * vetor[i + tw].X()[0];
        sum_sq -= vetor[i - 1].X()[0] * vetor[i - 1].X()[0];
        media = sqrt(sum_sq / tw);
        amp_mean[i] = media;
    }

    return amp_mean;
};

void ODE_analysis::Moving_Average_Draw(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<ODEpoint> resultado, std::vector<double> media, double epsilon, double time_step, double tw)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo Moving Average..." << std::endl;
    // FAZER GRÁFICO DA EDO
    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 3200, 1800);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    int n = resultado.size();

    double *tempo = new double[n];
    double *posicao = new double[n];
    double *velocidade = new double[n];

    for (int i = 0; i < n; i++)
    {
        tempo[i] = resultado[i].T();
        posicao[i] = resultado[i].X()[0];
        velocidade[i] = resultado[i].X()[1];
    }

    TGraph *gr = new TGraph(n, tempo, posicao);
    gStyle->SetTitleFontSize(0.05);
    gr->SetLineColor(50);
    gr->SetLineWidth(10);

    gr->SetTitle("Grafico Media Deslizante");
    gr->GetXaxis()->SetTitle(xaxis);
    gr->GetYaxis()->SetTitle(yaxis);
    gr->Draw("AC");

    double *media_arr = &media[0];

    double *tempo_novo = new double[n - 5000]; // Array com novo tamanho para descartar as medias que não foram calculadas nos extremos

    for (int i = 0; i < n - tw; i++)
    {
        tempo_novo[i] = tempo[(int)(i + tw / 2)];
    }

    TGraph *gr2 = new TGraph(n - tw, tempo_novo, media_arr);
    gr2->SetLineColor(62);
    gr2->SetLineWidth(10);
    gr2->Draw("same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillStyle(0);

    std::string s_epsilon;
    if (epsilon == 0.36)
        s_epsilon = "0.36";
    else
        s_epsilon = "0.2";

    std::string legenda_texto = "#delta = 0.4 #varepsilon = " + s_epsilon; // Criar string para a legenda com o valor da particula
    legenda->AddEntry(gr, legenda_texto.c_str(), "L");
    legenda->AddEntry(gr2, "Media Deslizante", "L");
    legenda->Draw();
    c->SaveAs(filename);

    delete[] tempo;
    delete[] posicao;
    delete[] velocidade;
    delete[] tempo_novo;

    delete c;
};

void ODE_analysis::Phase_Draw(const char *title, const char *s, int color, std::vector<ODEpoint> resultado)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o Retrato de Fase..." << std::endl;
    // FAZER GRÁFICO DA EDO
    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 6100, 3600);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    int n = resultado.size();

    double posicao[n];
    double velocidade[n];

    for (int i = 0; i < n; i++)
    {
        posicao[i] = resultado[i].X()[0];
        velocidade[i] = resultado[i].X()[1];
    }

    TGraph *gr = new TGraph(n, posicao, velocidade);
    gStyle->SetTitleFontSize(0.05);
    gr->SetLineColor(color);
    gr->SetLineWidth(10);

    gr->SetTitle(title);
    gr->GetYaxis()->SetTitle("#dot{#hat{#theta}}(#hat{t})");
    gr->GetXaxis()->SetTitle("#hat{#theta}(#hat{t})");

    gr->Draw("AC");

    c->SaveAs(s);

    delete c;
};

int ODE_analysis::Is_Stable(double *media, int n)
{
    int flag_estab = 1;
    for (int i = 0; i < n; i++)
    {
        if (media[i] > M_PI)
        {
            flag_estab = 0;
            break;
        }
    }
    return flag_estab;
};

void ODE_analysis::Draw_Stability(const char *filename, int numero, double step)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o Diagrama de Estabilidade..." << std::endl;

    int n;
    int count = 0;
    double step_delta = (10 - 0) / numero;
    double step_epsilon = (10 - 0) / numero;

    auto h = new TH2D("h", "Stability Regions", numero + 1, 0.02, 1.2, numero + 1, 0, 0.5);

    double est = 1;

    for (double delta = 0.00; delta <= 10; delta = delta + step_delta)
    {
        for (double epsilon = 0; epsilon <= 10 + step_epsilon; epsilon = epsilon + step_epsilon)
        {
            ODE_analysis pendulum(2, {0.1, 0});
            pendulum.SetFunction(0, [](ODEpoint p)
                                 { return p.X()[1]; });

            pendulum.SetFunction(1, [&](ODEpoint p)
                                 { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

            std::vector<ODEpoint> resultado = pendulum.LeapFrogImprovedSolver(120, step);
            std::vector<double> mov_avg = Moving_Average(resultado, 5e-3, 5 / 5e-3);

            n = mov_avg.size();

            double *media = &mov_avg[0];

            int value = Is_Stable(media, n);
            if (value == 0)
                est = 0.1;
            else
                est = 1;

            h->Fill(delta, epsilon, est);
        }
        count++;
        if (count % 5 == 0)
            std::cout << 100 * delta / 1.2 << "% " << std::endl;
    }
    std::cout << std::endl;

    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 1920, 1080);

    TStyle *mcStyle = new TStyle("mcStyle", "FC's Root Styles");
    mcStyle->SetPalette(kGreyYellow);

    h->SetTitle("Diagrama de Estabilidade para (#delta, #varepsilon)");
    h->GetXaxis()->SetTitle("#delta");
    h->GetYaxis()->SetTitle("#varepsilon");

    gStyle->SetOptStat(0);

    h->SetContour(1000);
    h->SetStats(0);
    h->SetLineColor(1);
    h->SetBarWidth(2);
    h->Draw("COL");

    // TGraph para as cores
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    gr1->SetLineColor(92);
    gr2->SetLineColor(112);
    gr1->SetLineWidth(20);
    gr2->SetLineWidth(20);
    // gr1->Draw("same");
    // gr2->Draw("same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillColor(0);
    legenda->AddEntry(gr1, "Solucoes Estaveis", "LP");
    legenda->AddEntry(gr2, "Solucoes Instaveis", "LP");
    legenda->Draw();

    c->Update();
    c->SaveAs(filename);

    delete c;
    delete h;
};

void ODE_analysis::Draw_Stability_Zoom(const char *filename, int numero, double step)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o Diagrama de Estabilidade Aumentado (Zoom In)..." << std::endl;

    int count = 0;
    int n;
    double step_delta = (0.35 - 0.2) / numero;
    double step_epsilon = (0.25 - 0.1) / numero;

    auto h = new TH2D("h", "Stability Regions", numero + 1, 0.2, 0.35, numero + 1, 0.1, 0.25);

    double est = 1;

    for (double delta = 0.2; delta <= 0.35; delta = delta + step_delta)
    {
        for (double epsilon = 0.1; epsilon <= 0.25 + step_epsilon; epsilon = epsilon + step_epsilon)
        {
            ODE_analysis pendulum(2, {0.1, 0});
            pendulum.SetFunction(0, [](ODEpoint p)
                                 { return p.X()[1]; });

            pendulum.SetFunction(1, [&](ODEpoint p)
                                 { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

            std::vector<ODEpoint> resultado = pendulum.LeapFrogImprovedSolver(100, step);
            std::vector<double> mov_avg = Moving_Average(resultado, step, 5 / step);

            n = mov_avg.size();

            double *media = &mov_avg[0];

            int value = Is_Stable(media, n);
            if (value == 0)
                est = 0.1;
            else
                est = 1;

            h->Fill(delta, epsilon, est);
        }
        count++;
        if (count % 5 == 0)
            std::cout << 100 * delta / 0.35 << "% " << std::endl;
    }

    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 6400, 3600);

    TStyle *mcStyle = new TStyle("mcStyle", "FC's Root Styles");
    mcStyle->SetPalette(kGreyYellow);

    h->SetTitle("Diagrama de Estabilidade para (#delta, #varepsilon)");
    h->GetXaxis()->SetTitle("#delta");
    h->GetYaxis()->SetTitle("#varepsilon");

    gStyle->SetOptStat(0);

    h->SetContour(1000);
    h->SetStats(0);
    h->SetLineColor(1);
    h->SetBarWidth(2);
    h->Draw("COL");

    // TGraph para as cores
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    gr1->SetLineColor(92);
    gr2->SetLineColor(112);
    gr1->SetLineWidth(20);
    gr2->SetLineWidth(20);
    // gr1->Draw("same");
    // gr2->Draw("same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillColor(0);
    legenda->AddEntry(gr1, "Solucoes Estaveis", "LP");
    legenda->AddEntry(gr2, "Solucoes Instaveis", "LP");
    legenda->Draw();

    c->Update();
    c->SaveAs(filename);

    delete c;
};