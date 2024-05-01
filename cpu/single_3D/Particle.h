#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <string>

#include "Vector.h"

class Particle
{
 protected:
 std::string symbol; // symbol of particle
 double mass; // mass of particle
 double radius; // radius of particle
 Vector position; // position of particle
 Vector velocity; // velocity of particle
 Vector halfvelocity; // halfvelocity of particle
 Vector force; // force on particle
 public:
 Particle(); // default constructor
 Particle(const std::string,const double,const double,const Vector,const Vector,const Vector,const Vector); // overloading constructor
 double getradius() const; // to get the value of radius of the particle
 double getmass() const; // to get the value of mass of the particle
 std::string getsymbol() const; // to get the value of the symbol of the particle
 Vector getposition() const; // to get the position of the particle
 void print() const; // to print all the parameters of the object of the Particle class
 void setforce(const Vector&); // to set the force of a particle
 Vector getforce() const; // to get the force of the particle
 void setposition(const Vector&); // to set the position of a particle
 void sethv(const Vector&); // to set the halfvelocity of a particle
 void setvelocity(const Vector&); // to set the velocity of a particle
 Vector getvelocity() const; // to get the velocity of the particle
 Vector gethv() const; // to get the halfvelocity of the particle
};

#endif

