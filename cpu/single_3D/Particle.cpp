#include "Particle.h"

Particle::Particle() //default constructor
{
 symbol="Atom";
 mass=1;
 radius=0.5;
 position=Vector(0,0,0);
 velocity=Vector(0,0,0);
 halfvelocity=Vector(0,0,0);
 force=Vector(0,0,0);
}

Particle::Particle(const std::string c,const double a,const double b,const Vector pos,const Vector vel,const Vector halfvel,const Vector forc) // overloading constructor
{
 symbol=c;
 mass=a;
 radius=b;
 position=pos;
 velocity=vel;
 halfvelocity=halfvel;
 force=forc;
}

double Particle::getradius() const // read the value of radius
{
 return radius;
}

double Particle::getmass() const // read the value of mass
{
 return mass;
}

std::string Particle::getsymbol() const // read the value of symbol
{
 return symbol;
}

Vector Particle::getposition() const // read the position of the particle
{
 return position;
}

Vector Particle::getforce() const // read the force of the particle
{
 return force;
}

void Particle::setforce(const Vector& vec) // set the force of particle
{
 this->force=vec;
}

void Particle::print() const // print the values of all the parameters of the objects of the Particle class
{
 std::cout<<"\nThe symbol of the protein is "<<symbol;
 std::cout<<"\nThe mass of the protein is "<<mass;
 std::cout<<"\nThe radius of the protein is "<<radius;
 std::cout<<"\nThe position of protein is "<<position.getx()<<"i+"<<position.gety()<<"j+"<<position.getz()<<"k";
 std::cout<<"\nThe velocity of protein is "<<velocity.getx()<<"i+"<<velocity.gety()<<"j+"<<velocity.getz()<<"k";
 std::cout<<"\nThe halfvelocity of protein is "<<halfvelocity.getx()<<"i+"<<halfvelocity.gety()<<"j+"<<halfvelocity.getz()<<"k";
 std::cout<<"\nThe force on protein is "<<force.getx()<<"i+"<<force.gety()<<"j+"<<force.getz()<<"k\n\n";
}

void Particle::setposition(const Vector& vec) // set the position of particle
{
 this->position=vec;
}

void Particle::sethv(const Vector& vec) // set the halfvelocity of particle
{
 this->halfvelocity=vec;
}

void Particle::setvelocity(const Vector& vec) // set the velocity of particle
{
 this->velocity=vec;
}

Vector Particle::getvelocity() const // read the velocity of the particle
{
 return velocity;
}

Vector Particle::gethv() const // read the halfvelocity of the particle
{
 return halfvelocity;
}

