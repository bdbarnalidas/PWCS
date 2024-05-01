#ifndef VERLET_H
#define VERLET_H

#include <vector>

#include "Protein.h"
#include "Cell.h"
#include "Boundary.h"
#include "Collision.h"
#include "Simulate.h"

class Verlet
{
 public:
 Verlet(); // default constructor
 void VelocityVerlet(Cell&,std::vector<Protein>&,const Simulate&,const Boundary&,const Collision&) const; // run Velocity-Verlet algorithm
};

void update_position(std::vector<Protein>&,const double); // update position of the proteins
void update_halfvelocity(std::vector<Protein>&,const double); // update halfvelocity of the proteins
void update_velocity(std::vector<Protein>&,const double); // update velocity of the proteins

#endif
