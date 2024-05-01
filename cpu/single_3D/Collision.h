#ifndef COLLISION_H
#define COLLISION_H

#include <vector>

#include "Protein.h"
#include "Simulate.h"

class Collision
{
 private:
 int COR; // coefficient of restitution
 public:
 Collision(); // default constructor
 Collision(const int); // overloading constructor
 void ResolveCollision(std::vector<Protein>&) const; // detect and resolve collision (both single-body collision and multi-body collision)
 void PredictResolveCollision(std::vector<Protein>&,const std::vector<double>&) const; // detect and resolve collision (both single-body collision and multi-body collision) for pseudoparticles
 void PredictCollision(std::vector<Protein>&,const Simulate&) const; // predict future collisions to avoid overlapping of particles
 int get_COR() const; // return the value of COR
};

void SingleRestitution(Protein&,Protein&,const int); // COR equation for single-body collision
void ElasticMultiRestitution(Protein&,Protein&,const int,const Vector&); // COR equation for multi-body collision
void InelasticMultiRestitution(std::vector<Protein>&); // COR equation for multi-body collision

#endif
