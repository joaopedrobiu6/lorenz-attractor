#include "ODEpoint.h"

Xvar::Xvar(int num)
{
    x.resize(num);
};

Xvar::Xvar(const Xvar &v)
{
    x = v.x;
};

Xvar::~Xvar()
{
    x.clear();
};

Xvar &Xvar::operator=(const Xvar &v)
{
    if (this != &v)
    {
        x = v.x;
    }
    return *this;
};

Xvar &Xvar::operator+(const Xvar &v)
{
    for (int i = 0; i < x.size(); i++)
    {
        x[i] += v.x[i];
    }
    return *this;
};

double &Xvar::operator[](int a)
{
    return x[a];
};

Xvar operator*(double k, const Xvar &v)
{
    Xvar c(v);
    for (int i = 0; i < (c.X()).size(); i++)
    {
        c.X()[i] *= k;
    }
    return c;
};

std::ostream &operator<<(std::ostream &s, const Xvar &v)
{
    Xvar c(v);
    s << "Number of Dependent Variables:" << (c.X()).size() << std::endl;

    for (int i = 0; i < (c.X()).size(); ++i)
    {
        s << (c.X())[i] << std::endl;
    }
    return s;
};

std::vector<double> &Xvar::X()
{
    return x;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////ODE POINTS//////////////////////////////////
////////////////////////////////////////////////////////////////////////

ODEpoint::ODEpoint()
{
    t = -1;
};

void ODEpoint::SetODEpoint(double t_, Xvar &p)
{
    t = t_;
    x.resize(p.X().size());
    for (int i = 0; i < p.X().size(); i++)
    {
        x[i] = p.X()[i];
    }
};
void ODEpoint::SetODEpoint(double t_, const std::initializer_list<double> &v)
{
    t = t_;
    x.resize(v.size());
    x = v;
};
void ODEpoint::SetODEpoint(double t_, std::vector<double> v)
{
    t = t_;
    x.resize(v.size());

    // x = v;

    for (int i = 0; i < v.size(); i++)
    {
        x[i] = v[i];
    }
};
