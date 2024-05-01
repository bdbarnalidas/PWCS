#include "Simulate.h"

Simulate::Simulate() // default constructor
{
 timeInterval=0.005;
 nsteps=200;
 duration=nsteps*timeInterval;
 type=-1.0; // -1 for attractive and +1 for repulsive
 nrow=1;
 ncol=1;
}

Simulate::Simulate(const double a,const int b,const int c,const int d) // overloading constructor
{
 timeInterval=a;
 nsteps=b;
 duration=timeInterval*nsteps;
 nrow=c;
 ncol=d;
}

double Simulate::get_timeInterval() const // read the value of timeInterval
{
 return this->timeInterval;
}

int Simulate::get_nsteps() const // get the no of steps
{
 return this->nsteps;
}
 
double Simulate::get_type() const // get the type of force-field
{
 return this->type;
} 

int Simulate::get_nrow() const // get the no of rows
{
 return this->nrow;
}

int Simulate::get_ncol() const // get the no of cols
{
 return this->ncol;
}


