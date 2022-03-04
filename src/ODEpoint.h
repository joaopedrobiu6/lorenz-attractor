#ifndef __ODEPOINT__
#define __ODEPOINT__

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <map> 

class Xvar // classe para armazenar todas as variávies dependentes - posição, velocidade, aceleração, etc
{
public:
    // Xvar() = default;
    Xvar(int num = 1);                                        // number of dependent variables
    Xvar(const std::vector<double> &v) : x(v) { ; }           // passing vector //ESTES JA ESTAO DECLARADOS AQUI!!!!
    Xvar(const std::initializer_list<double> &l) : x(l) { ; } // using initializer list to build object: X({1,2}) //ESTES JA ESTAO DECLARADOS AQUI!!!!
    Xvar(const Xvar &);                                       // copy constructor
    ~Xvar();

    Xvar &operator=(const Xvar &); // assignment operator
    Xvar &operator+(const Xvar &); // operator+

    double &operator[](int); // X[i]

    friend Xvar operator*(double, const Xvar &); // scalar*X
    // Xvar operator*()
    friend std::ostream &operator<<(std::ostream &, const Xvar &);

    std::vector<double> &X(); // accessor to x, para o pendulo x[0] = posição, x[1] = velocidade, x[2] = aceleração

protected:
    std::vector<double> x;
};

class ODEpoint : public Xvar // classe para armazenar a variável independente - tempo
{
public:
    ODEpoint();
    ODEpoint(double t_, Xvar a_) : t(t_), Xvar(a_) { ; }
    ODEpoint(double t_, const std::vector<double> &v) : t(t_), Xvar(v) { ; }
    ODEpoint(double t_, const std::initializer_list<double> &v) : t(t_), Xvar(v) { ; }

    void SetODEpoint(double t_, Xvar &p);
    void SetODEpoint(double t_, const std::initializer_list<double> &v);
    void SetODEpoint(double t_, std::vector<double> v);

    double &T() { return t; } // accessor to time

private:
    double t; // time
};

#endif