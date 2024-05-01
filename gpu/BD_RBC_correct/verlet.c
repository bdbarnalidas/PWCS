__global__ void velocityverlet(double *dev_halfvelocity,double *dev_velocity,double *dev_force,double *dev_position,double *mass,double *radius,int *startrow,int *lastrow,int npart,int *startx,int *starty,int *startz,int *width,int *height,int *depth,double *organelle,int dev_pupu)
{
 int i=0,j=0,m=0,n=0,t=0;
 int flag=0; // flag is 1 if any particle has collided with the wall of the box
 double *temp1 = (double *)malloc(NDIM*sizeof(double)); // temporary array 1 to store the coordinates of i-th particle as vector
 double *temp2 = (double *)malloc(NDIM*sizeof(double)); // temporary array 2 to store the coordinates of j-th particle as vector
 double *temp3 = (double *)malloc(NDIM*sizeof(double)); // temporary array 3 to store the vector (b-a)
 double temp[N];
 double *temp5 = (double *)malloc(NDIM*sizeof(double)); // temporary array 1 to store the coordinates of i-th particle as vector
 double *temp6 = (double *)malloc(NDIM*sizeof(double)); // temporary array 2 to store the coordinates of j-th particle as vector
 double *temp7 = (double *)malloc(NDIM*sizeof(double)); // temporary array 3 to store the vector (b-a)
 double *temp8 = (double *)malloc(NDIM*sizeof(double)); // temporary array 3 to store the vector temp*temp3
 double tempx[N],tempy[N],tempz[N];
 double tempa[N],tempb[N],tempc[N];
 double huha[N],s[N],forc[N],so[N];


/**************************************************** Elastic collision (equal masses) ************************************************************/
 for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5-5;i+=5)
 {
  temp1[0]=dev_position[i+2];temp1[1]=dev_position[i+3];temp1[2]=dev_position[i+4]; // coordinates of i-th particle
  for(j=i+5;j<lastrow[blockIdx.x]*5;j+=5)
  {
   temp2[0]=dev_position[j+2];temp2[1]=dev_position[j+3];temp2[2]=dev_position[j+4]; // coordinates of j-th particle
   huha[blockIdx.x]=(temp1[0]-temp2[0])*(temp1[0]-temp2[0])+(temp1[1]-temp2[1])*(temp1[1]-temp2[1])+(temp1[2]-temp2[2])*(temp1[2]-temp2[2]);
   s[blockIdx.x]=sqrt(huha[blockIdx.x]);
   if((s[blockIdx.x]<=(radius[i/5]+radius[j/5])) && (!organelle[i/5]) && (!organelle[j/5]))// collision, i particle, j particle 
   {
    if((int)dev_position[i+1]==15) 
    printf("collision\n");

    tempa[blockIdx.x]=dev_velocity[i+2];tempb[blockIdx.x]=dev_velocity[i+3];tempc[blockIdx.x]=dev_velocity[i+4];
    dev_velocity[i+2]=(CR*mass[j/5]*(dev_velocity[j+2]-dev_velocity[i+2])+mass[i/5]*dev_velocity[i+2]+mass[j/5]*dev_velocity[j+2])/(mass[i/5]+mass[j/5]);
    dev_velocity[i+3]=(CR*mass[j/5]*(dev_velocity[j+3]-dev_velocity[i+3])+mass[i/5]*dev_velocity[i+3]+mass[j/5]*dev_velocity[j+3])/(mass[i/5]+mass[j/5]);
    dev_velocity[i+4]=(CR*mass[j/5]*(dev_velocity[j+4]-dev_velocity[i+4])+mass[i/5]*dev_velocity[i+4]+mass[j/5]*dev_velocity[j+4])/(mass[i/5]+mass[j/5]);
    dev_velocity[j+2]=(CR*mass[i/5]*(tempa[blockIdx.x]-dev_velocity[j+2])+mass[i/5]*tempa[blockIdx.x]+mass[j/5]*dev_velocity[j+2])/(mass[i/5]+mass[j/5]);
    dev_velocity[j+3]=(CR*mass[i/5]*(tempb[blockIdx.x]-dev_velocity[j+3])+mass[i/5]*tempb[blockIdx.x]+mass[j/5]*dev_velocity[j+3])/(mass[i/5]+mass[j/5]);
    dev_velocity[j+4]=(CR*mass[i/5]*(tempc[blockIdx.x]-dev_velocity[j+4])+mass[i/5]*tempc[blockIdx.x]+mass[j/5]*dev_velocity[j+4])/(mass[i/5]+mass[j/5]);
   }
   else if((s[blockIdx.x]<=(radius[i/5]+radius[j/5])) && (organelle[i/5]) && (!organelle[j/5]))// collision, i organelle, j particle
   {
    dev_velocity[j+2]=0-dev_velocity[j+2];
    dev_velocity[j+3]=0-dev_velocity[j+3];
    dev_velocity[j+4]=0-dev_velocity[j+4];
   }
   else if((s[blockIdx.x]<=(radius[i/5]+radius[j/5])) && (!organelle[i/5]) && (organelle[j/5]))// collision, j organelle, i particle
   {
    dev_velocity[i+2]=0-dev_velocity[i+2];
    dev_velocity[i+3]=0-dev_velocity[i+3];
    dev_velocity[i+4]=0-dev_velocity[i+4];
   }
  }
 }

/*******************************************************************/
 for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
 {
  dev_halfvelocity[i+2]=dev_velocity[i+2]+(dev_force[i+2]*delta_t)/(2*mass[i/5]); // update half-velocity
  dev_position[i+2]+=dev_halfvelocity[i+2]*delta_t; // update position
  dev_halfvelocity[i+3]=dev_velocity[i+3]+(dev_force[i+3]*delta_t)/(2*mass[i/5]); // update half-velocity
  dev_position[i+3]+=dev_halfvelocity[i+3]*delta_t; // update position
  dev_halfvelocity[i+4]=dev_velocity[i+4]+(dev_force[i+4]*delta_t)/(2*mass[i/5]); // update half-velocity
  dev_position[i+4]+=dev_halfvelocity[i+4]*delta_t; // update position
 }

/************************************************ flush force ************************************************************************/
 for(m=startrow[blockIdx.x]*5;m<lastrow[blockIdx.x]*5;m+=5)
 {
  dev_force[m+2]=0;
  dev_force[m+3]=0;
  dev_force[m+4]=0;
 }

/************************************************ update dev_force ****************************************************************************/
 
 for(m=startrow[blockIdx.x]*5;m<lastrow[blockIdx.x]*5-5;m+=5)
 {
  temp5[0]=dev_position[m+2];temp5[1]=dev_position[m+3];temp5[2]=dev_position[m+4]; // coordinates of i-th particle
  for(n=m+5;n<lastrow[blockIdx.x]*5;n+=5)
  {
   temp6[0]=dev_position[n+2];temp6[1]=dev_position[n+3];temp6[2]=dev_position[n+4]; // coordinates of j-th particle
   forc[blockIdx.x]=(temp5[0]-temp6[0])*(temp5[0]-temp6[0])+(temp5[1]-temp6[1])*(temp5[1]-temp6[1])+(temp5[2]-temp6[2])*(temp5[2]-temp6[2]);
   so[blockIdx.x]=sqrt(forc[blockIdx.x]);
   temp[blockIdx.x]=(G*mass[m/5]*mass[n/5])/(so[blockIdx.x]*so[blockIdx.x]*so[blockIdx.x]);
   for(i=0;i<NDIM;i++)
   temp7[i]=temp5[i]-temp6[i];
   for(i=0;i<NDIM;i++)
   temp8[i]=temp7[i]*temp[blockIdx.x];

    dev_force[m+2]+=temp8[0];dev_force[m+3]+=temp8[1];dev_force[m+4]+=temp8[2];
    dev_force[n+2]-=temp8[0];dev_force[n+3]-=temp8[1];dev_force[n+4]-=temp8[2];
  }
 }

/*********************************************************************************************************************************************/ 
 for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
 {
  dev_velocity[i+2]=dev_halfvelocity[i+2]+(dev_force[i+2]*delta_t)/(2*mass[i/5]); // update dev_velocity
  dev_velocity[i+3]=dev_halfvelocity[i+3]+(dev_force[i+3]*delta_t)/(2*mass[i/5]); // update dev_velocity
  dev_velocity[i+4]=dev_halfvelocity[i+4]+(dev_force[i+4]*delta_t)/(2*mass[i/5]); // update dev_velocity
 }
 
 
/**************************************************** Rigid boundary conditions *****************************************************************/
/****************************** if the particles are outside the box, put them inside the box ***************************************************/


 if(BC==1) // if we specify rigid boundary conditions
 {
  for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
  {
   if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 1 
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 2 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 3 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 4 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 5 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 6 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 7 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 8 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 12 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 23 
   { 
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 34 
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 14 
   { 
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 56
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 67
   {
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 78
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 58
   { 
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 15
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 26
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 37
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 48
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 1265
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 4378
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 3267
   { 
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x]; 
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) && (dev_position[i+3]<=height[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 4158
   { 
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]>height[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 1234
   { 
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
   }
   else if((dev_position[i+2]<=width[blockIdx.x]) && (dev_position[i+3]<starty[blockIdx.x]) && (dev_position[i+4]<=depth[blockIdx.x])) // Quadrant - 5678
   { 
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
   }
  } 
  if(flag) // update force if the particles have collided with the wall
  {
   for(m=startrow[blockIdx.x]*5;m<lastrow[blockIdx.x]*5;m+=5) // flush force
   {
    dev_force[m+2]=0;
    dev_force[m+3]=0;
    dev_force[m+4]=0;
   }

   for(m=startrow[blockIdx.x]*5;m<lastrow[blockIdx.x]*5-5;m+=5) // update force
   {
    temp5[0]=dev_position[m+2];temp5[1]=dev_position[m+3];temp5[2]=dev_position[m+4]; // coordinates of i-th particle
    for(n=m+5;n<lastrow[blockIdx.x]*5;n+=5)
    {
     temp6[0]=dev_position[n+2];temp6[1]=dev_position[n+3];temp6[2]=dev_position[n+4]; // coordinates of j-th particle
     forc[blockIdx.x]=(temp5[0]-temp6[0])*(temp5[0]-temp6[0])+(temp5[1]-temp6[1])*(temp5[1]-temp6[1])+(temp5[2]-temp6[2])*(temp5[2]-temp6[2]);
     so[blockIdx.x]=sqrt(forc[blockIdx.x]);
     temp[blockIdx.x]=(G*mass[m/5]*mass[n/5])/(so[blockIdx.x]*so[blockIdx.x]*so[blockIdx.x]);
     for(t=0;t<NDIM;t++)
     temp7[t]=temp5[t]-temp6[t];
     for(t=0;t<NDIM;t++)
     temp8[t]=temp7[t]*temp[blockIdx.x];

      dev_force[m+2]+=temp8[0];dev_force[m+3]+=temp8[1];dev_force[m+4]+=temp8[2];
      dev_force[n+2]-=temp8[0];dev_force[n+3]-=temp8[1];dev_force[n+4]-=temp8[2];
    }
   }
   flag=0;
  }
 }
 dev_pupu=funpupi(5);
/***************************************************************************************************************************/

 free(temp1);free(temp2);free(temp3);free(temp5);free(temp6);free(temp7);free(temp8);  
}

int funpupi(int x)
{
 return x*2;
}
