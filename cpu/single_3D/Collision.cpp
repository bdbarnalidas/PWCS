#include "Collision.h"
#include "Vector.h"

Collision::Collision() //  default constructor
{
 COR=0; // 1 for elastic and 0 for inelastic
}

Collision::Collision(const int a) // overloading constructor
{
 COR=a;
}

void SingleRestitution(Protein& obja,Protein& objb,const int a) // COR equation for single-body collision
{
 Vector first;
 Vector second;
 Vector third;
 Vector fourth;
 Vector fifth;
 first=scalar_mult((obja.getmass()/(obja.getmass()+objb.getmass())),obja.getvelocity());
 second=scalar_mult((objb.getmass()/(obja.getmass()+objb.getmass())),objb.getvelocity());
 third=scalar_mult(((objb.getmass()/(obja.getmass()+objb.getmass()))*a),vector_sub(objb.getvelocity(),obja.getvelocity()));
 fourth=vector_add(vector_add(first,second),third);
 third=scalar_mult(((obja.getmass()/(obja.getmass()+objb.getmass()))*a),vector_sub(obja.getvelocity(),objb.getvelocity()));
 fifth=vector_add(vector_add(first,second),third);
 obja.setvelocity(fourth); // v_a final velocity of a is updated
 objb.setvelocity(fifth); // v_b final velocity of b is updated
}

void ElasticMultiRestitution(Protein& obja,Protein& objb,const int a,const Vector& vec) // COR equation for multi-body collision
{
 Vector first;
 Vector second;
 Vector third;
 Vector fourth;
 Vector fifth;
 first=scalar_mult((obja.getmass()/(obja.getmass()+objb.getmass())),vec);
 second=scalar_mult((objb.getmass()/(obja.getmass()+objb.getmass())),objb.getvelocity());
 third=scalar_mult(((objb.getmass()/(obja.getmass()+objb.getmass()))*a),vector_sub(objb.getvelocity(),vec));
 fourth=vector_add(vector_add(first,second),third);
 third=scalar_mult(((obja.getmass()/(obja.getmass()+objb.getmass()))*a),vector_sub(vec,objb.getvelocity()));
 fifth=vector_add(vector_add(first,second),third);
 obja.setvelocity(vector_add(fourth,obja.getvelocity())); // resultant of the previous vectors v_a final velocity of a is updated
 objb.setvelocity(fifth); // v_b final velocity of b is updated
}

void InelasticMultiRestitution(std::vector<Protein>& proteins) // COR equation for multi-body collision
{
 int size=proteins.size(); // get the size of the proteins vector
 Vector vec(0,0,0); // zero vector
 Vector temp1; // temporary vector 1
 Vector temp2; // temporary vector 2
 Vector temp;
 for(int i=0;i<size-1;i++)
 {
  int flag=0;
  Vector temp_vel=proteins[i].getvelocity(); // temp_vel temporarily storing the velocity of i-th protein
  for(int j=i+1;j<size;j++)
  {
   double dist=getmagnitude(vector_sub(proteins[i].getposition(),proteins[j].getposition())); // distance between the i-th protein and j-th protein
   if(dist<=(proteins[i].getradius()+proteins[j].getradius())) // collision between i-th protein and j-th protein
   {
    if(!flag)
    {
     temp1=scalar_mult((proteins[i].getmass()/(proteins[i].getmass()+proteins[j].getmass())),temp_vel);
     temp2=scalar_mult((proteins[j].getmass()/(proteins[i].getmass()+proteins[j].getmass())),proteins[j].getvelocity());
     temp=vector_add(temp1,temp2);
     proteins[i].setvelocity(temp); // set i-th protein's velocity 
     proteins[j].setvelocity(temp); // set j-th protein's velocity 
     flag=1;
    }
    else if(flag)
    {
     proteins[j].setvelocity(temp);
    }
   }
  }
 }
}


   
  

void Collision::ResolveCollision(std::vector<Protein>& proteins) const //detect and resolve collision (both single-body collision and multi-body collision)
{
 int size=proteins.size(); // get the size of the proteins vector
 for(int i=0;i<size-1;i++) // take the i-th protein
 {
  int flag=0; // flag for differentiating between single-body and multi-body collision
  Vector temp_vel=proteins[i].getvelocity(); // temp_vel temporarily storing the velocity of i-th protein
  for(int j=i+1;j<size;j++) // take the j-th protein
  {
   double dist=getmagnitude(vector_sub(proteins[i].getposition(),proteins[j].getposition())); // distance between the i-th protein and j-th protein
   if(dist<=(proteins[i].getradius()+proteins[j].getradius())) // collision between i-th protein and j-th protein
   {
    if(COR==1) // elastic collision
    {
     if(flag) // multi-body elastic collision
     ElasticMultiRestitution(proteins[i],proteins[j],COR,temp_vel); // update velocitites using the COR equation for multi-body collision
     else // single-body collision
     SingleRestitution(proteins[i],proteins[j],COR); // update velocitites using the COR equation for single-body collision
     flag=1; // already collided with a single body, from now on collisions will be multi-body
    }
    else // inelastic collision
    InelasticMultiRestitution(proteins); // both single-body and multi-body collisions
   }
  }
 }
}

void Collision::PredictResolveCollision(std::vector<Protein>& proteins,const std::vector<double>& dummy) const //detect and resolve collision (both single-body collision and multi-body collision)
{
 int size=proteins.size(); // get the size of the proteins vector
 for(int i=0;i<size-1;i++) // take the i-th protein
 {
  int flag=0; // flag for differentiating between single-body and multi-body collision
  Vector temp_vel=proteins[i].getvelocity(); // temp_vel temporarily storing the velocity of i-th protein
  for(int j=i+1;j<size;j++) // take the j-th protein
  {
   double dist=getmagnitude(vector_sub(proteins[i].getposition(),proteins[j].getposition())); // distance between the i-th protein and j-th protein
   if(dist<=(dummy[i]+dummy[j])) // collision between i-th protein and j-th protein
   {
    if(COR==1) // elastic collision
    {
     if(flag) // multi-body collision
     ElasticMultiRestitution(proteins[i],proteins[j],COR,temp_vel); // update velocitites using the COR equation for multi-body collision
     else // single-body collision
     SingleRestitution(proteins[i],proteins[j],COR); // update velocitites using the COR equation for single-body collision
     flag=1; // already collided with a single body, from now on collisions will be multi-body
    }
    else // inelastic collision
    InelasticMultiRestitution(proteins); // both single-body and multi-body collisions
   }
  }
 }
}

void Collision::PredictCollision(std::vector<Protein>& proteins,const Simulate& sim) const // predict future collisions to avoid overlapping of particles
{
 int size=proteins.size(); // get the size of the proteins vector
 double temp; // temporary variable to store the radii of the pseudospheres
 std::vector<double>dummy(size); // create a dummy vector of the same size as the proteins vector which will store the radii of the future predicted particles
 for(int i=0;i<size;i++)
 {
  temp=(getmagnitude(proteins[i].getvelocity())*sim.get_timeInterval())+proteins[i].getradius();
  dummy.push_back(temp);
 }
 PredictResolveCollision(proteins,dummy);
}
 
int Collision::get_COR() const // get the value of COR
{
 return this->COR;
}   
