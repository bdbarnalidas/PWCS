#include "Boundary.h"

Boundary::Boundary() // default constructor
{

}

Boundary::Boundary(const double x,const double y,const double z) // overloading constructor
{
 a=Vector(x,0.0,0.0);
 b=Vector(0.0,y,0.0);
 c=Vector(0.0,0.0,z); 
} 

double sign(const double a) // extract the sign of a number
{
 if (a > 0) return 1.0;
 if (a < 0) return -1.0;
 return 0;
}
 
void Boundary::rigid(std::vector<Protein>& proteins) const // apply rigid boundary conditions
{
 int size=proteins.size(); // get the size of the proteins vector
 Vector temp_pos; // temporary vector for storing the position 
 Vector temp_vel; // temporary vector for storing the velocity
 
 Vector nbc=scalar_mult(sign(dot(a,cross(b,c)))*(1.0/getmagnitude(cross(b,c))),cross(b,c)); // normal to plane bc
 Vector nca=scalar_mult(sign(dot(b,cross(a,c)))*(1.0/getmagnitude(cross(a,c))),cross(a,c)); // normal to plane ac
 Vector nab=scalar_mult(sign(dot(c,cross(a,b)))*(1.0/getmagnitude(cross(a,b))),cross(a,b)); // normal to plane ab

 std::vector<int>connected(size); //  create a connected vector to keep a track of the particles connected to each other
 for(int i=0;i<size;i++)
 connected[i]=0;
 for(int i=0;i<size;i++) // initialize the connected vector by the particle numbers
 connected.push_back(i);
 for(int i=0;i<size-1;i++)
 {
  for(int j=i+1;j<size;j++)
  {
   double dist=getmagnitude(vector_sub(proteins[i].getposition(),proteins[j].getposition())); // distance between the i-th protein and j-th protein
   if(dist<=(proteins[i].getradius()+proteins[j].getradius())) // i-th and j-th proteins are connected to each other
   connected[j]=connected[i];
  }
 }

 
 for(int i=0;i<size;i++)
 { 
  Vector temp=proteins[i].getposition();
  if((dot(nbc,a))<(dot(nbc,proteins[i].getposition()))) // check if particle is further than a's component in direction parallel to nbc
  {
   temp_pos=vector_add(scalar_mult((dot(nbc,proteins[i].getposition())-dot(nbc,a))*(-2.0),nbc),proteins[i].getposition());
   proteins[i].setposition(temp_pos);
   temp_vel=vector_add(scalar_mult((dot(nbc,proteins[i].getvelocity()))*(-2.0),nbc),proteins[i].getvelocity());
   proteins[i].setvelocity(temp_vel);

   for(int j=0;j<size;j++)
   {
    if((connected[j]==connected[i]) && (i!=j))
    {
     temp_pos=vector_add(scalar_mult((dot(nbc,temp)-dot(nbc,a))*(-2.0),nbc),proteins[j].getposition()); // give the displacement of i to the connected j's
     proteins[j].setposition(temp_pos);
     proteins[j].setvelocity(temp_vel); // set the velocity of j as the same velocity as i
    }
   }
  }
  else if(dot(nbc,proteins[i].getposition())<0)
  {
   temp_pos=vector_add(scalar_mult(dot(nbc,proteins[i].getposition())*(-2.0),nbc),proteins[i].getposition());
   proteins[i].setposition(temp_pos);
   temp_vel=vector_add(scalar_mult(dot(nbc,proteins[i].getvelocity())*(-2.0),nbc),proteins[i].getvelocity());
   proteins[i].setvelocity(temp_vel);
 
   for(int j=0;j<size;j++)
   {
    if((connected[j]==connected[i]) && (i!=j))
    {
     temp_pos=vector_add(scalar_mult(dot(nbc,temp)*(-2.0),nbc),proteins[j].getposition());
     proteins[j].setposition(temp_pos);
     proteins[j].setvelocity(temp_vel);
    }
   }
  }
  else if((dot(nca,b))<(dot(nca,proteins[i].getposition())))
  {
   temp_pos=vector_add(scalar_mult((dot(nca,proteins[i].getposition())-dot(nca,b))*(-2.0),nca),proteins[i].getposition());
   proteins[i].setposition(temp_pos);
   temp_vel=vector_add(scalar_mult((dot(nca,proteins[i].getvelocity()))*(-2.0),nca),proteins[i].getvelocity());
   proteins[i].setvelocity(temp_vel);
   
   for(int j=0;j<size;j++)
   {
    if((connected[j]==connected[i]) && (i!=j))
    {
     temp_pos=vector_add(scalar_mult((dot(nca,temp)-dot(nca,b))*(-2.0),nca),proteins[j].getposition());
     proteins[j].setposition(temp_pos);
     proteins[j].setvelocity(temp_vel);
    }
   }
  }
  else if(dot(nca,proteins[i].getposition())<0)
  {
   temp_pos=vector_add(scalar_mult(dot(nca,proteins[i].getposition())*(-2.0),nca),proteins[i].getposition());
   proteins[i].setposition(temp_pos);
   temp_vel=vector_add(scalar_mult(dot(nca,proteins[i].getvelocity())*(-2.0),nca),proteins[i].getvelocity());
   proteins[i].setvelocity(temp_vel);

   for(int j=0;j<size;j++)
   {
    if((connected[j]==connected[i]) && (i!=j))
    {
     temp_pos=vector_add(scalar_mult(dot(nca,temp)*(-2.0),nca),proteins[j].getposition());
     proteins[j].setposition(temp_pos);
     proteins[j].setvelocity(temp_vel);
    }
   }
  }
  else if((dot(nab,c))<(dot(nab,proteins[i].getposition())))
  {
   temp_pos=vector_add(scalar_mult((dot(nab,proteins[i].getposition())-dot(nab,c))*(-2.0),nab),proteins[i].getposition());
   proteins[i].setposition(temp_pos);
   temp_vel=vector_add(scalar_mult((dot(nab,proteins[i].getvelocity()))*(-2.0),nab),proteins[i].getvelocity());
   proteins[i].setvelocity(temp_vel);

   for(int j=0;j<size;j++)
   {
    if((connected[j]==connected[i]) && (i!=j))
    {
     temp_pos=vector_add(scalar_mult((dot(nab,temp)-dot(nab,c))*(-2.0),nab),proteins[j].getposition());
     proteins[j].setposition(temp_pos);
     proteins[j].setvelocity(temp_vel);
    }
   }
  }
  else if(dot(nab,proteins[i].getposition())<0)
  {
   temp_pos=vector_add(scalar_mult(dot(nab,proteins[i].getposition())*(-2.0),nab),proteins[i].getposition());
   proteins[i].setposition(temp_pos);
   temp_vel=vector_add(scalar_mult(dot(nab,proteins[i].getvelocity())*(-2.0),nab),proteins[i].getvelocity());
   proteins[i].setvelocity(temp_vel);

   for(int j=0;j<size;j++)
   {
    if((connected[j]==connected[i]) && (i!=j))
    {
     temp_pos=vector_add(scalar_mult(dot(nab,temp)*(-2.0),nab),proteins[j].getposition());
     //double calx=temp_pos.getx()+((wholecell.getbreadth()+1)*colno);
     //temp_pos.setx(calx);
     //double caly=temp_pos.gety()+((wholecell.getlength()+1)*rowno);
     //temp_pos.sety(caly);
     proteins[j].setposition(temp_pos);
     proteins[j].setvelocity(temp_vel);
    }
   }
  }
 }
}


Vector Boundary::geta() const // return a
{
 return a;
}

Vector Boundary::getb() const
{
 return b;
}

Vector Boundary::getc() const
{
 return c;
}

void Boundary::periodic(std::vector<Protein>& proteins) const // apply periodic boundary conditions
{
 int size=proteins.size(); // get the size of the proteins vector
 Vector temp_pos; // temporary vector for storing the position 
 Vector temp_vel; // temporary vector for storing the velocity
 Vector nbc=scalar_mult(sign(dot(a,cross(b,c)))*(1.0/getmagnitude(cross(b,c))),cross(b,c)); // normal to plane bc
 Vector nca=scalar_mult(sign(dot(b,cross(a,c)))*(1.0/getmagnitude(cross(a,c))),cross(a,c)); // normal to plane ac
 Vector nab=scalar_mult(sign(dot(c,cross(a,b)))*(1.0/getmagnitude(cross(a,b))),cross(a,b)); // normal to plane ab
 for(int i=0;i<size;i++)
 { 
  if((dot(nbc,a))<(dot(nbc,proteins[i].getposition()))) // check if particle is further than a's component in direction parallel to nbc
  {
   temp_pos=vector_sub(proteins[i].getposition(),a);
   proteins[i].setposition(temp_pos);
  }
  else if(dot(nbc,proteins[i].getposition())<0)
  {
   temp_pos=vector_add(proteins[i].getposition(),a);
   proteins[i].setposition(temp_pos);
  }
  else if((dot(nca,b))<(dot(nca,proteins[i].getposition())))
  {
   temp_pos=vector_sub(proteins[i].getposition(),b);
   proteins[i].setposition(temp_pos);
  }
  else if(dot(nca,proteins[i].getposition())<0)
  {
   temp_pos=vector_add(proteins[i].getposition(),b);
   proteins[i].setposition(temp_pos);
  }
  else if((dot(nab,c))<(dot(nab,proteins[i].getposition())))
  {
   temp_pos=vector_sub(proteins[i].getposition(),c);
   proteins[i].setposition(temp_pos);
  }
  else if(dot(nab,proteins[i].getposition())<0)
  {
   temp_pos=vector_add(proteins[i].getposition(),c);
   proteins[i].setposition(temp_pos);
  }
 }
}

void Boundary::seta(const double d,const double e,const double f) // set vector a of the boundary
{
 a=Vector(d,e,f);
}

void Boundary::setb(const double d,const double e,const double f) // set vector b of the boundary
{
 b=Vector(d,e,f);
}

void Boundary::setc(const double d,const double e,const double f) // set vector c of the boundary
{
 c=Vector(d,e,f);
}




