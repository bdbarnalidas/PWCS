#ifndef PROTEIN_H
#define PROTEIN_H

#include "Particle.h"

class Protein : public Particle
{
 public:
 Protein(); // default constructor 
 Protein(const std::string,const double,const double,const Vector,const Vector,const Vector,const Vector); // overloading constructor 
};
 
#endif
