__host__ __device__ double** vanderwalforce(double **force,double **pos,double *mass,int npart)
{
 

 int i=0,j=0,t=0;
 double s,temp;
 double *temp1 = (double *)malloc(NDIM*sizeof(double)); // temporary array 1 to store the coordinates of i-th particle as vector
 double *temp2 = (double *)malloc(NDIM*sizeof(double)); // temporary array 2 to store the coordinates of j-th particle as vector
 double *temp3 = (double *)malloc(NDIM*sizeof(double)); // temporary array 3 to store the vector (b-a)
 double *temp4 = (double *)malloc(NDIM*sizeof(double)); // temporary array 3 to store the vector temp*temp3
 for(i=0;i<npart-1;i++)
 {
  temp1[0]=pos[i][2];temp1[1]=pos[i][3];temp1[2]=pos[i][4]; // coordinates of i-th particle
  for(j=i+1;j<npart;j++)
  {
   

   temp2[0]=pos[j][2];temp2[1]=pos[j][3];temp2[2]=pos[j][4]; // coordinates of j-th particle
   s=DISTV(s,temp1,temp2); // s=((b-a).length)
   
   temp=(G*mass[i]*mass[j])/(pow(s,3)); // temp=(G*m1*m2)/(s^3)
   temp3=SUBV(temp3,temp1,temp2); // temp3=(b-a)
   temp4=MULVS(temp4,temp3,temp); // temp4=temp*temp3
  

   for(t=0;t<NDIM;t++)
   {
    

    force[i][t+2]+=temp4[t]; // to calculate net force due to many body interaction (i->j so +ve)
    force[j][t+2]-=temp4[t]; // to calculate net force due to many body interaction (j->i so -ve)
   }
  }
 }

 return force;
}
