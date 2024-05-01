__global__ void velocityverlet(double *x_interval,double *y_interval,double *dev_halfvelocity,double *dev_velocity,double *dev_force,double *dev_position,double *mass,double *radius,int *startrow,int *lastrow,int npart,double *startx,double *starty,double *startz,double *width,double *height,double *depth)
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
 int particle_i,particle_j,particle_ia,particle_ja,particle_jb,particle_jc,part_i,part_j;


/**************************************************** Elastic collision (equal masses) ************************************************************/
 for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5-5;i+=5)
 {
  temp1[0]=dev_position[i+2];temp1[1]=dev_position[i+3];temp1[2]=dev_position[i+4]; // coordinates of i-th particle
  particle_i=(int)dev_position[i+1];

  for(j=i+5;j<lastrow[blockIdx.x]*5;j+=5)
  {
   temp2[0]=dev_position[j+2];temp2[1]=dev_position[j+3];temp2[2]=dev_position[j+4]; // coordinates of j-th particle
   particle_j=(int)dev_position[j+1];

   huha[blockIdx.x]=(temp1[0]-temp2[0])*(temp1[0]-temp2[0])+(temp1[1]-temp2[1])*(temp1[1]-temp2[1])+(temp1[2]-temp2[2])*(temp1[2]-temp2[2]);
   s[blockIdx.x]=sqrt(huha[blockIdx.x]);
   if(s[blockIdx.x]<=(radius[particle_i]+radius[particle_j])) // collision so swap the velocities
   {
    printf("collision\n");
    tempa[blockIdx.x]=dev_velocity[i+2];tempb[blockIdx.x]=dev_velocity[i+3];tempc[blockIdx.x]=dev_velocity[i+4];
    dev_velocity[i+2]=(CR*mass[particle_j]*(dev_velocity[j+2]-dev_velocity[i+2])+mass[particle_i]*dev_velocity[i+2]+mass[particle_j]*dev_velocity[j+2])/(mass[particle_i]+mass[particle_j]);
    dev_velocity[i+3]=(CR*mass[particle_j]*(dev_velocity[j+3]-dev_velocity[i+3])+mass[particle_i]*dev_velocity[i+3]+mass[particle_j]*dev_velocity[j+3])/(mass[particle_i]+mass[particle_j]);
    dev_velocity[i+4]=(CR*mass[particle_j]*(dev_velocity[j+4]-dev_velocity[i+4])+mass[particle_i]*dev_velocity[i+4]+mass[particle_j]*dev_velocity[j+4])/(mass[particle_i]+mass[particle_j]);
    dev_velocity[j+2]=(CR*mass[particle_i]*(tempa[blockIdx.x]-dev_velocity[j+2])+mass[particle_i]*tempa[blockIdx.x]+mass[particle_j]*dev_velocity[j+2])/(mass[particle_i]+mass[particle_j]);
    dev_velocity[j+3]=(CR*mass[particle_i]*(tempb[blockIdx.x]-dev_velocity[j+3])+mass[particle_i]*tempb[blockIdx.x]+mass[particle_j]*dev_velocity[j+3])/(mass[particle_i]+mass[particle_j]);
    dev_velocity[j+4]=(CR*mass[particle_i]*(tempc[blockIdx.x]-dev_velocity[j+4])+mass[particle_i]*tempc[blockIdx.x]+mass[particle_j]*dev_velocity[j+4])/(mass[particle_i]+mass[particle_j]);
   }
  }
 }

/*******************************************************************/
 for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
 {
  particle_ia=(int)dev_position[i+1];
  dev_halfvelocity[i+2]=dev_velocity[i+2]+(dev_force[i+2]*delta_t)/(2*mass[particle_ia]); // update half-velocity
  dev_position[i+2]+=dev_halfvelocity[i+2]*delta_t; // update position
  dev_halfvelocity[i+3]=dev_velocity[i+3]+(dev_force[i+3]*delta_t)/(2*mass[particle_ia]); // update half-velocity
  dev_position[i+3]+=dev_halfvelocity[i+3]*delta_t; // update position
  dev_halfvelocity[i+4]=dev_velocity[i+4]+(dev_force[i+4]*delta_t)/(2*mass[particle_ia]); // update half-velocity
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
  particle_ja=(int)dev_position[m+1];

  for(n=m+5;n<lastrow[blockIdx.x]*5;n+=5)
  {
   temp6[0]=dev_position[n+2];temp6[1]=dev_position[n+3];temp6[2]=dev_position[n+4]; // coordinates of j-th particle
   particle_jb=(int)dev_position[n+1];

   forc[blockIdx.x]=(temp5[0]-temp6[0])*(temp5[0]-temp6[0])+(temp5[1]-temp6[1])*(temp5[1]-temp6[1])+(temp5[2]-temp6[2])*(temp5[2]-temp6[2]);
   so[blockIdx.x]=sqrt(forc[blockIdx.x]);
   temp[blockIdx.x]=(G*mass[particle_ja]*mass[particle_jb])/(so[blockIdx.x]*so[blockIdx.x]*so[blockIdx.x]);
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
  particle_jc=(int)dev_position[i+1];
  dev_velocity[i+2]=dev_halfvelocity[i+2]+(dev_force[i+2]*delta_t)/(2*mass[particle_jc]); // update dev_velocity
  dev_velocity[i+3]=dev_halfvelocity[i+3]+(dev_force[i+3]*delta_t)/(2*mass[particle_jc]); // update dev_velocity
  dev_velocity[i+4]=dev_halfvelocity[i+4]+(dev_force[i+4]*delta_t)/(2*mass[particle_jc]); // update dev_velocity
 }
 
 
/**************************************************** Rigid boundary conditions *****************************************************************/
/****************************** if the particles are outside the box, put them inside the box ***************************************************/


 if(BC==1) // if we specify rigid boundary conditions
 {
  for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
  {
   if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 1 
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
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 2 
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
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 3 
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
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 4 
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
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 5 
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
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 6 
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
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 7 
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
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 8 
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
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 12 
   {// only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 23 
   { // reflection for boundary boxes - UP & RIGHT 
     // inter-core communication among rest
    if(((startx[blockIdx.x]==x_interval[2]) & (width[blockIdx.x]==x_interval[3])) | ((starty[blockIdx.x]==y_interval[2]) & (height[blockIdx.x]==y_interval[3]))) // reflect
    { 
     flag=1;
     dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
     dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
     tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
     dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
     tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
     dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x+(col+1);
     dev_force[i]=blockIdx.x+(col+1);
     dev_velocity[i]=blockIdx.x+(col+1);
     dev_halfvelocity[i]=blockIdx.x+(col+1);
    }
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 34 
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
    dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 14 
   { // reflection for boundary boxes - UP & LEFT 
     // inter-core communication among rest
    //printf("\nQUADRANT 14\n");
    //printf("dev_position[i+3]=%f\n",dev_position[i+3]);
    //printf("height=%f\n",height[blockIdx.x]);
    if(((startx[blockIdx.x]==x_interval[0]) & (width[blockIdx.x]==x_interval[1])) | ((starty[blockIdx.x]==y_interval[2]) & (height[blockIdx.x]==y_interval[3]))) // reflect
    {
     flag=1;
     dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
     dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
     tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
     dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
     tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
     dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    } 
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x+(col-1);
     dev_force[i]=blockIdx.x+(col-1);
     dev_velocity[i]=blockIdx.x+(col-1);
     dev_halfvelocity[i]=blockIdx.x+(col-1);
    }
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 56
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 67
   { // reflection for boundary boxes - RIGHT & DOWN 
     // inter-core communication among rest
    if(((startx[blockIdx.x]==x_interval[2]) & (width[blockIdx.x]==x_interval[3])) | ((starty[blockIdx.x]==y_interval[0]) & (height[blockIdx.x]==y_interval[1]))) // reflect
    {
     flag=1;
     dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
     dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
     tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
     dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
     tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
     dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x-(col-1);
     dev_force[i]=blockIdx.x-(col-1);
     dev_velocity[i]=blockIdx.x-(col-1);
     dev_halfvelocity[i]=blockIdx.x-(col-1);
    }
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 78
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
    dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 58
   { // reflection for boundary boxes - LEFT & DOWN 
     // inter-core communication among rest
    if(((startx[blockIdx.x]==x_interval[0]) & (width[blockIdx.x]==x_interval[1])) | ((starty[blockIdx.x]==y_interval[0]) & (height[blockIdx.x]==y_interval[1]))) // reflect
    {
     flag=1;
     dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
     dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
     tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
     dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
     tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
     dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x-(col+1);
     dev_force[i]=blockIdx.x-(col+1);
     dev_velocity[i]=blockIdx.x-(col+1);
     dev_halfvelocity[i]=blockIdx.x-(col+1);
    }
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 15
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 26
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 37
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
    dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 48
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
    dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]>depth[blockIdx.x])) // Quadrant - 1265
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempz[blockIdx.x]=dev_position[i+4]-depth[blockIdx.x];
    dev_position[i+4]=depth[blockIdx.x]-tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]<startz[blockIdx.x])) // Quadrant - 4378
   { // only reflection happens...no inter-core communication needed
    flag=1;
    dev_velocity[i+4]=0-dev_velocity[i+4]; // reverse the k-th component of the d_velocity vector
    tempz[blockIdx.x]=startz[blockIdx.x]-dev_position[i+4];
    dev_position[i+4]=startz[blockIdx.x]+tempz[blockIdx.x];
   }
   else if((dev_position[i+2]>width[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 3267
   { // reflection for boundary boxes - RIGHT
     // inter-core communication among rest
    if((startx[blockIdx.x]==x_interval[2]) & (width[blockIdx.x]==x_interval[3])) // reflection
    {
     flag=1;
     dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
     tempx[blockIdx.x]=dev_position[i+2]-width[blockIdx.x];
     dev_position[i+2]=width[blockIdx.x]-tempx[blockIdx.x]; 
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x+1;
     dev_force[i]=blockIdx.x+1;
     dev_velocity[i]=blockIdx.x+1;
     dev_halfvelocity[i]=blockIdx.x+1;
    }
   }
   else if((dev_position[i+2]<startx[blockIdx.x]) & (dev_position[i+3]>starty[blockIdx.x] & dev_position[i+3]<height[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 4158
   { // reflection for boundary boxes - LEFT
     // inter-core communication among rest
    if((startx[blockIdx.x]==x_interval[0]) & (width[blockIdx.x]==x_interval[1])) // reflection
    {
     flag=1;
     dev_velocity[i+2]=0-dev_velocity[i+2]; // reverse the i-th component of the d_velocity vector
     tempx[blockIdx.x]=startx[blockIdx.x]-dev_position[i+2];
     dev_position[i+2]=startx[blockIdx.x]+tempx[blockIdx.x];
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x-1;
     dev_force[i]=blockIdx.x-1;
     dev_velocity[i]=blockIdx.x-1;
     dev_halfvelocity[i]=blockIdx.x-1;
    }
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]>height[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 1234
   { // reflection for boundary boxes - UP
     // inter-core communication among rest
    if((starty[blockIdx.x]==y_interval[2]) & (height[blockIdx.x]==y_interval[3])) // reflection
    {
     flag=1;
     dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
     tempy[blockIdx.x]=dev_position[i+3]-height[blockIdx.x];
     dev_position[i+3]=height[blockIdx.x]-tempy[blockIdx.x];
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x+col;
     dev_force[i]=blockIdx.x+col;
     dev_velocity[i]=blockIdx.x+col;
     dev_halfvelocity[i]=blockIdx.x+col;
    }
   }
   else if((dev_position[i+2]>startx[blockIdx.x] & dev_position[i+2]<width[blockIdx.x]) & (dev_position[i+3]<starty[blockIdx.x]) & (dev_position[i+4]>startz[blockIdx.x] & dev_position[i+4]<depth[blockIdx.x])) // Quadrant - 5678
   { // reflection for boundary boxes - DOWN
     // inter-core communication among rest
    if((starty[blockIdx.x]==y_interval[0]) & (height[blockIdx.x]==y_interval[1])) // reflection
    {
     flag=1;
     dev_velocity[i+3]=0-dev_velocity[i+3]; // reverse the j-th component of the d_velocity vector
     tempy[blockIdx.x]=starty[blockIdx.x]-dev_position[i+3];
     dev_position[i+3]=starty[blockIdx.x]+tempy[blockIdx.x];
    }
    else
    {
     flag=0;
     dev_position[i]=blockIdx.x-col;
     dev_force[i]=blockIdx.x-col;
     dev_velocity[i]=blockIdx.x-col;
     dev_halfvelocity[i]=blockIdx.x-col;
    }
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
    part_i=(int)dev_position[m+1];

    for(n=m+5;n<lastrow[blockIdx.x]*5;n+=5)
    {
     temp6[0]=dev_position[n+2];temp6[1]=dev_position[n+3];temp6[2]=dev_position[n+4]; // coordinates of j-th particle
     part_j=(int)dev_position[n+1];

     forc[blockIdx.x]=(temp5[0]-temp6[0])*(temp5[0]-temp6[0])+(temp5[1]-temp6[1])*(temp5[1]-temp6[1])+(temp5[2]-temp6[2])*(temp5[2]-temp6[2]);
     so[blockIdx.x]=sqrt(forc[blockIdx.x]);
     temp[blockIdx.x]=(G*mass[part_i]*mass[part_j])/(so[blockIdx.x]*so[blockIdx.x]*so[blockIdx.x]);
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


/*******************************************************************/
 /*for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
 {
  dev_halfvelocity[i+2]=dev_velocity[i+2]+(dev_force[i+2]*delta_t)/(2*mass[i/5]); // update half-velocity
  dev_position[i+2]+=dev_halfvelocity[i+2]*delta_t; // update position
  dev_halfvelocity[i+3]=dev_velocity[i+3]+(dev_force[i+3]*delta_t)/(2*mass[i/5]); // update half-velocity
  dev_position[i+3]+=dev_halfvelocity[i+3]*delta_t; // update position
  dev_halfvelocity[i+4]=dev_velocity[i+4]+(dev_force[i+4]*delta_t)/(2*mass[i/5]); // update half-velocity
  dev_position[i+4]+=dev_halfvelocity[i+4]*delta_t; // update position
 }*/

/************************************************ flush force ************************************************************************/
 /*for(m=startrow[blockIdx.x]*5;m<lastrow[blockIdx.x]*5;m+=5)
 {
  dev_force[m+2]=0;
  dev_force[m+3]=0;
  dev_force[m+4]=0;
 }*/

/************************************************ update dev_force ****************************************************************************/

 /*for(m=startrow[blockIdx.x]*5;m<lastrow[blockIdx.x]*5-5;m+=5)
 {
  temp5[0]=dev_position[m+2];temp5[1]=dev_position[m+3];temp5[2]=dev_position[m+4]; // coordinates of i-th particle
  for(n=m+5;n<lastrow[blockIdx.x]*5;n+=5)
  {
   temp6[0]=dev_position[n+2];temp6[1]=dev_position[n+3];temp6[2]=dev_position[n+4]; // coordinates of j-th particle
   for(i=0;i<NDIM;i++)
   tmp += (temp5[i]-temp6[i]) *(temp5[i]-temp6[i]);
   so=sqrt(tmp);
   temp=(G*mass[m/5]*mass[n/5])/(so*so*so);
   for(i=0;i<NDIM;i++)
   temp7[i]=temp5[i]-temp6[i];
   for(i=0;i<NDIM;i++)
   temp8[i]=temp7[i]*temp;

    dev_force[m+2]+=temp8[0];dev_force[m+3]+=temp8[1];dev_force[m+4]+=temp8[2];
    dev_force[n+2]-=temp8[0];dev_force[n+3]-=temp8[1];dev_force[n+4]-=temp8[2];
  }
 }*/

/*********************************************************************************************************************************************/
 /*for(i=startrow[blockIdx.x]*5;i<lastrow[blockIdx.x]*5;i+=5)
 {
  dev_velocity[i+2]=dev_halfvelocity[i+2]+(dev_force[i+2]*delta_t)/(2*mass[i/5]); // update dev_velocity
  dev_velocity[i+3]=dev_halfvelocity[i+3]+(dev_force[i+3]*delta_t)/(2*mass[i/5]); // update dev_velocity
  dev_velocity[i+4]=dev_halfvelocity[i+4]+(dev_force[i+4]*delta_t)/(2*mass[i/5]); // update dev_velocity
 }*/


/***************************************************************************************************************************/

 free(temp1);free(temp2);free(temp3);free(temp5);free(temp6);free(temp7);free(temp8);  
}
