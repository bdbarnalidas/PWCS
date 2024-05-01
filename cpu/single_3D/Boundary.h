#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>

#include "Protein.h"
#include "Vector.h"

class Boundary
{
 private:
 Vector a; // convert box dimensions to vectors
 Vector b;
 Vector c;
 public:
 Boundary(); // default constructor
 Boundary(const double,const double,const double); // overloading constructor
 void rigid(std::vector<Protein>&) const; // apply rigid boundary conditions
 Vector geta() const; // return a
 Vector getb() const; // return b
 Vector getc() const; // return c
 void periodic(std::vector<Protein>&) const; // apply periodic boundary conditions
 void seta(const double,const double,const double); // set a
 void setb(const double,const double,const double); // set b
 void setc(const double,const double,const double); // set c
};

double sign(const double); // extract the sign of a number

#endif
