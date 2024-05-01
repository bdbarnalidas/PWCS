#include "Protein.h"

Protein::Protein() // default constructor
{
 symbol="Prot";
 mass=1;
 radius=0.5;
 position=Vector(0,0,0);
 velocity=Vector(0,0,0);
 halfvelocity=Vector(0,0,0);
 force=Vector(0,0,0);
}

Protein::Protein(const std::string g,const double a,const double b,const Vector c,const Vector d,const Vector e,const Vector f) // overloading constructor
{
 symbol=g;
 mass=a;
 radius=b;
 position=c;
 velocity=d;
 halfvelocity=e;
 force=f;
}

