#ifndef FORCE_H
#define FORCE_H

#define G 1 // the gravitational constant

#include <vector>

#include "Vector.h"
#include "Protein.h"
#include "Simulate.h"

class Force
{
 private:
 public:
 Force(); // default constructor
 void GravitationalForces(std::vector<Protein>&,const Simulate&) const; // calculate force on all the proteins
};

Vector GravitationalForce(const Protein&,const Protein&); // calculate force acting between 2 proteins

#endif
