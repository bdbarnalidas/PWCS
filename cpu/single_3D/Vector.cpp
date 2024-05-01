#include <cmath>

#include "Vector.h"

Vector::Vector() // default constructor
{
 xcomp=0;
 ycomp=0;
 zcomp=0;
}

Vector::Vector(const double a,const double b,const double c) // overloading constructor
{
 xcomp=a;
 ycomp=b;
 zcomp=c;
}

double Vector::getx () const // read the value of x component
{
 return xcomp;
}

double Vector::gety () const // read the value of y component
{
 return ycomp;
}

double Vector::getz () const // read the value of z component
{
 return zcomp;
}

Vector getdistance(const Vector& obja,const Vector& objb) // calculate distance between 2 vectors
{
 double x,y,z;
 x=obja.getx()-objb.getx(); // difference between the components along X-axis
 y=obja.gety()-objb.gety(); // difference between the components along Y-axis
 z=obja.getz()-objb.getz(); // difference between the components along Z-axis
 Vector vec(x,y,z); // create the difference vector
 return vec;
}

double getmagnitude(const Vector& obj) // calculate the magnitude of a vector
{
 double x,y,z,length;
 x=pow(obj.getx(),2); // square the X-component
 y=pow(obj.gety(),2); // square the Y-component
 z=pow(obj.getz(),2); // square the Z-component
 length=sqrt(x+y+z); // magnitude of the vector
 return length;
}

Vector scalar_mult(const double a,const Vector& obj) // multiply vector by a scalar
{
 double x,y,z;
 x=obj.getx()*a; // multiply the X-component of the vector by a
 y=obj.gety()*a; // multiply the Y-component of the vector by a
 z=obj.getz()*a; // multiply the Z-component of the vector by a
 Vector vec(x,y,z); // create the result vector
 return vec;
}

Vector vector_add(const Vector& obja,const Vector& objb) // add 2 vectors
{
 double x,y,z;
 x=obja.getx()+objb.getx(); // add the X-componenets of the 2 vectors
 y=obja.gety()+objb.gety(); // add the Y-components of the 2 vectors
 z=obja.getz()+objb.getz(); // add the Z-componenets of the 2 vectors
 Vector vec(x,y,z); // create the result vector
 return vec;
}

Vector vector_sub(const Vector& obja,const Vector& objb) // subtract 2 vectors
{
 double x,y,z;
 x=obja.getx()-objb.getx(); // subtract the X-componenets of the 2 vectors
 y=obja.gety()-objb.gety(); // subtract the Y-components of the 2 vectors
 z=obja.getz()-objb.getz(); // subtract the Z-componenets of the 2 vectors
 Vector vec(x,y,z); // create the result vector
 return vec;
}
 
Vector cross(const Vector& obja,const Vector& objb) // cross product of 2 vectors
{
 double x,y,z;
 x=obja.gety()*objb.getz()-obja.getz()*objb.gety(); // a_y*b_z-a_z*b_y
 y=obja.getz()*objb.getx()-obja.getx()*objb.getz(); // a_z*b_x-a_x*b_z
 z=obja.getx()*objb.gety()-obja.gety()*objb.getx(); // a_x*b_y-a_y*b_x
 Vector vec(x,y,z); // result vector
 return vec;
}

double dot(const Vector& obja,const Vector& objb) // dot product of 2 vectors
{
 double res;
 res=obja.getx()*objb.getx()+obja.gety()*objb.gety()+obja.getz()*objb.getz(); // a_x*b_x+a_y*b_y+a_z*b_z
 return res;
}


 
