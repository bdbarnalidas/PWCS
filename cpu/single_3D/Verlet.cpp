#include "Verlet.h"
#include "Force.h"

Verlet::Verlet() // default constructor
{
 
}

void Verlet::VelocityVerlet(Cell& cell,std::vector<Protein>& proteins,const Simulate& sim,const Boundary& boundary,const Collision& collision) const // run Velocity-Verlet algorithm
{
 collision.ResolveCollision(proteins); // detect and resolve collisions 
 collision.PredictCollision(proteins,sim);
 boundary.rigid(proteins); // apply rigid boundary conditions
 update_position(proteins,sim.get_timeInterval()); // update position of the proteins
 update_halfvelocity(proteins,sim.get_timeInterval()); // update halfvelocity of the proteins
 Force().GravitationalForces(proteins,sim); // update force of the proteins
 update_velocity(proteins,sim.get_timeInterval()); // update velocity of the proteins
}

void update_position(std::vector<Protein>& proteins,const double delta_t) // update position of the proteins
{
 int size=proteins.size(); // get the size of the proteins vector
 for(int i=0;i<size;i++)
 {
  Vector v; // dummy vector for storing the result
  double temp; 
  v=vector_add(proteins[i].getposition(),scalar_mult(delta_t,proteins[i].getvelocity()));
  temp=(delta_t*delta_t)/(2*proteins[i].getmass());
  v=vector_add(v,scalar_mult(temp,proteins[i].getforce()));
  proteins[i].setposition(v); // new position of the i-th protein
 }
}
 
void update_halfvelocity(std::vector<Protein>& proteins,const double delta_t) // update halfvelocity of the proteins 
{
 int size=proteins.size(); // get the size of the proteins vector
 for(int i=0;i<size;i++)
 {
  Vector v; // dummy vector for storing the result
  double temp; 
  temp=delta_t/(2*proteins[i].getmass());
  v=vector_add(proteins[i].getvelocity(),scalar_mult(temp,proteins[i].getforce()));
  proteins[i].sethv(v);
 }
}

void update_velocity(std::vector<Protein>& proteins,const double delta_t) // update velocity of the proteins 
{
 int size=proteins.size(); // get the size of the proteins vector
 for(int i=0;i<size;i++)
 {
  Vector v; // dummy vector for storing the result
  double temp; 
  temp=delta_t/(2*proteins[i].getmass());
  v=vector_add(proteins[i].gethv(),scalar_mult(temp,proteins[i].getforce()));
  proteins[i].setvelocity(v);
 }
}

 
