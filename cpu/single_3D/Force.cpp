#include <cmath>

#include "Force.h"

Force::Force() // default constructor
{

}

Vector GravitationalForce(const Protein& a,const Protein& b) // calculate force between 2 proteins
{
 Vector f; // force between 2 proteins
 double num;
 num=G*a.getmass()*b.getmass(); // numerator of force-field
 Vector diff; 
 diff=getdistance(a.getposition(),b.getposition()); // difference between the positional vectors of a and b
 double magn,cube_magn,front;
 magn=getmagnitude(diff); // magnitude of the difference vector
 cube_magn=pow(magn,3); // cube of the magnitude of the difference vector
 front=num/cube_magn; // scalar of the product
 f=scalar_mult(front,diff); // the force acting between 2 proteins a and b
 return f;
}

void Force::GravitationalForces(std::vector<Protein>& proteins,const Simulate& sim) const // calculate force on each protein
{
 int size=proteins.size(); // get the no of proteins in the cell
 Vector vec(0,0,0); //  dummy vector
 for(int t=0;t<size;t++) // flush the forces of the proteins
 proteins[t].setforce(vec);
 for(int i=0;i<size-1;i++)
 {
  for(int j=i+1;j<size;j++)
  {
   Vector vf=GravitationalForce(proteins[i],proteins[j]);
   vf=scalar_mult(sim.get_type(),vf); // attractive force field
   proteins[i].setforce(vector_add(proteins[i].getforce(),vf));
   proteins[j].setforce(vector_sub(proteins[j].getforce(),vf));
  }
 }
}



