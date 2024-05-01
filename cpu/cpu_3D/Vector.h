#ifndef VECTOR_H
#define VECTOR_H

class Vector
{
 private:
 double xcomp; // component along X-axis
 double ycomp; // component along Y-axis
 double zcomp; // component along Z-axis
 public:
 Vector(); // default constructor
 Vector(const double,const double,const double); // overloading constructor
 double getx () const; // to get the value of the x component
 double gety () const; // to get the value of the y component
 double getz () const; // to get the value of the z component
 void setx(const double); // to set the value of xcomp
 void sety(const double); // to set the value of ycomp
 void setz(const double); // to set the value of zcomp
};

Vector getdistance(const Vector&,const Vector&); // calculate distance vector between 2 vectors
double getmagnitude(const Vector&); // calculate the magnitude of a vector
Vector scalar_mult(const double,const Vector&); // multiply vector by a scalar
Vector vector_add(const Vector&,const Vector&); // add 2 vectors
Vector vector_sub(const Vector&,const Vector&); // subtract 2 vectors
Vector cross(const Vector&,const Vector&); // cross product of 2 vectors
double dot(const Vector&,const Vector&); // dot product of 2 vectors

#endif
