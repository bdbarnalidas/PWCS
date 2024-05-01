#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <string.h>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>


#define G 1 // gravitational constant
#define delta_t 0.1 // timestep of integration
int total_time=5000; // time of integration
double rangex=0,rangey=0,rangez=0; // ranges in x and y and z direction
#define row 6 // no of rows
#define col 6 // no of columns
int boxcount=0; // count the boxes
int boxes;
double width=10; // width of each box
double height=10; // height of each box 
double depth=60; // depth of each box
double widthc=10; // width of each box copy
double heightc=10; // height of each box copy
double depthc=60; // depth of each box copy
double rad=0.5; // modify boxes according to the radius of the particle
double startx=0;
double starty=0;
double startz=0;
int npart; // no of particles
int sumpart=0; // sum of the no of particles
#define N 36 // launch cuda cores = no of boxes in the cell
double arr_startx[N],arr_starty[N],arr_startz[N];
double arr_width[N],arr_height[N],arr_depth[N]; 
int arr_startrow[N],arr_lastrow[N];

#include "vectdefs.h"
int x[NDIM]; // boxvector in x-direction
int y[NDIM]; // boxvector in y-direction
int z[NDIM]; // boxvector in z-direction
#include "maths.c"
#include "force.c"
#include "verlet.c"

int main()
{
 boxes=row*col; // total no of boxes
 int i=0,j=0,temporary,seed=6,k_2=0,k_3=0,k_4=0,count=0,colcount=0;
 int t=0;
 //int startrow=0,lastrow=0;
 double xr; // xr will store the random number
 FILE *fp;
 fp=fopen("output.xyz","w+"); // output coordinate file
 int rearr=0,count_core=0;

/*************************************** Boundary boxes ***************************************************************************************/
 
 double x_interval[4];
 double y_interval[4];
 x_interval[0]=startx;
 y_interval[0]=starty;
 x_interval[1]=width+0.5; // Make transitions continuous
 y_interval[1]=height+0.5;
 double calx;
 calx=col*width+(col-1);
 double caly;
 caly=row*height+(row-1);
 x_interval[2]=calx-width-0.5;x_interval[3]=calx;
 y_interval[2]=caly-height-0.5;y_interval[3]=caly;

 printf("x_interval\n");
 for(i=0;i<4;i++)
 printf("%f ",x_interval[i]);
 printf("\n");
 printf("y_interval\n");
 for(i=0;i<4;i++)
 printf("%f ",y_interval[i]);
 printf("\n");

/*************************** Dynamically allocate particle matrix which stores no of particles of boxes ***************************************/
 int *particle = (int *)malloc(boxes*sizeof(int));
 for(i=0;i<boxes;i++)
 {
  srand(seed);
  temporary = (rand() % 10);
  if(temporary==0) // npart should not be 0 and try to keep the min no of particles to be 2
  temporary+=2;
  else if(temporary==1)
  temporary++;
  particle[i]=temporary; 
  seed++;
 }
/****************************************** Print particle matrix  ***************************************************************************/
 
 printf("Printing particle matrix --------------------------\n");
 for(i=0;i<boxes;i++)
 printf("%d ",particle[i]);
 printf("\n");

/****************************************** Allocate arr_startrow ************************************************************/
 arr_startrow[0]=0;
 for(i=1;i<N;i++)
 {
  arr_startrow[i]=arr_startrow[i-1]+particle[i-1];
 }

 printf("Printing arr_startrow--------\n");
 for(i=0;i<N;i++)
 printf("%d ",arr_startrow[i]);
 printf("\n");

/****************************************** Sum particle matrix  *****************************************************************************/
 for(i=0;i<boxes;i++)
 sumpart+=particle[i];
 //printf("%d\n",sumpart);
 npart=sumpart;
/****************************************** Allocate arr_lastrow ***************************************************************************/
 arr_lastrow[N-1]=npart-1;
 for(i=N-2;i>=0;i--)
 arr_lastrow[i]=arr_lastrow[i+1]-particle[i+1];
 for(i=0;i<N;i++)
 arr_lastrow[i]=arr_lastrow[i]+1;
 
 printf("Printing arr_lastrow--------\n");
 for(i=0;i<N;i++)
 printf("%d ",arr_lastrow[i]);
 printf("\n");

/****************************************** Dynamically allocate entire box position matrix **************************************************/
 double **position = (double **)malloc(sumpart*sizeof(double *));
 for(i=0;i<sumpart;i++)
 position[i]= (double *)malloc(5*sizeof(double));
 for(i=0;i<sumpart;i++)
 {
  for(j=0;j<5;j++)
  position[i][j]=0;
 }

/****************************************** Dynamically allocate entire force  matrix **************************************************/

 double **force = (double **)malloc(sumpart*sizeof(double *)); // dynamically allocate the force matrix
 for(i=0;i<sumpart;i++)
 force[i]= (double *)malloc(5*sizeof(double));
 for(i=0;i<sumpart;i++) // initialize the force matrix to 0
 {
  for(j=0;j<5;j++)
  {
   force[i][j]=0;
  }
 }

/****************************************** Dynamically allocate entire velocity matrix **************************************************/

 double **velocity = (double **)malloc(sumpart*sizeof(double *)); // dynamically allocate the velocity matrix
 for(i=0;i<sumpart;i++)
 velocity[i]= (double *)malloc(5*sizeof(double));
 for(i=0;i<sumpart;i++) // initialize the velocity matrix to 0
 {
  for(j=0;j<5;j++)
  {
   velocity[i][j]=0;
  }
 }

/****************************************** Dynamically allocate entire halfvelocity matrix **************************************************/

 double **halfvelocity = (double **)malloc(sumpart*sizeof(double *)); // dynamically allocate the half-velocity matrix
 for(i=0;i<sumpart;i++)
 halfvelocity[i]= (double *)malloc(5*sizeof(double));
 for(i=0;i<sumpart;i++) // initialize the half-velocity matrix to 0
 {
  for(j=0;j<5;j++)
  {
   halfvelocity[i][j]=0;
  }
 }

/****************************************** Fill the 1st column of the position matrix ****************************************************/
 boxcount=0;j=0;
 while(boxcount<boxes)
 {
  for(i=0;i<particle[boxcount];i++)
  {
   position[j][0]=boxcount;
   force[j][0]=boxcount;
   velocity[j][0]=boxcount;
   halfvelocity[j][0]=boxcount;
   j++;
  }
  boxcount++;
 }
/****************************************** Fill the 2nd column of the position matrix ****************************************************/
 j=0;
 for(i=0;i<sumpart;i++)
 {
  position[j][1]=i;
  force[j][1]=i;
  velocity[j][1]=i;
  halfvelocity[j][1]=i;
  j++;
 }
/************************************************* Test print the position matrix ************************************************************/ 
 /*for(i=0;i<sumpart;i++)
 {
  for(j=0;j<5;j++)
  printf("%f ",position[i][j]);
  printf("\n");
 } 
 printf("\n\n");*/  
/****************************************** Generate random coordinates for all the boxes ****************************************************/

 i=0;boxcount=0;k_2=0;k_3=0;k_4=0;
 while(i<row)
 {
  colcount=0;
  while(colcount<col)
  {
   x[0]=startx;x[1]=width;x[2]=0;
   y[0]=starty;y[1]=height;y[2]=0;
   z[0]=startz;z[1]=depth;z[2]=0;
   rangex=abs(x[1]-x[0]);
   rangey=abs(y[1]-y[0]);
   rangez=abs(z[1]-z[0]);
   for(j=0;j<particle[boxcount];j++) // generate x-coordinates
   {
    xr = (float)rand()/(float)(RAND_MAX/rangex)+x[0];
    position[k_2][2]=xr;k_2++;
   }
   for(j=0;j<particle[boxcount];j++) // generate y-coordinates
   {
    xr = (float)rand()/(float)(RAND_MAX/rangey)+y[0];
    position[k_3][3]=xr;k_3++;
   }
   for(j=0;j<particle[boxcount];j++) // generate z-coordinates
   {
    xr = (float)rand()/(float)(RAND_MAX/rangez)+z[0];
    position[k_4][4]=xr;k_4++;
   }
   boxcount++;
   colcount++;
   startx=width+1;width=startx+widthc;
  }
  i++;startx=0;width=widthc;starty=height+1;height=starty+heightc;  
 }
/************************************************* Test print the position matrix ************************************************************/ 
 
 for(i=0;i<sumpart;i++)
 {
  for(j=0;j<5;j++)
  printf("%f ",position[i][j]);
  printf("\n");
 }
 
/************************************************* Write the file ************************************************************/
 fprintf(fp,"%d\n",sumpart); 
 for(i=0;i<sumpart;i++)
 {
  fprintf(fp,"Atom%d ",i);
  fprintf(fp,"%f %f %f\n",position[i][2],position[i][3],position[i][4]);
 }  

/*************************** Dynamically allocate organelle_pos matrix which stores coordinates of organelles ***********************************/
 /*double **organelle_pos = (double **)malloc(npart*sizeof(double *));
 for(i=0;i<npart;i++)
 organelle_pos[i]= (double *)malloc(5*sizeof(double));

 for(i=0;i<npart;i++)
 {
  for(j=0;j<5;j++)
  organelle_pos[i][j]=0;
 }*/

/*************************Dynamically allocate mass matrix which stores masses of particles ***************************************/
 double *mass = (double *)malloc(sumpart*sizeof(double));
 for(i=0;i<sumpart;i++)
 mass[i]=1; // all the particles are of mass 1

 //mass[1]=999; // organelle
 /*mass[7]=999;
 mass[12]=999;
 mass[15]=999;*/

/*************************** Dynamically allocate radius matrix which stores radii of particles **************************************/
 double *radius = (double *)malloc(sumpart*sizeof(double));
 for(i=0;i<sumpart;i++)
 radius[i]=0.5; // all the particles are of radius 0.5

 //radius[1]=1.0; // organelle
 /*radius[7]=1.0;
 radius[12]=1.0;
 radius[15]=1.0;*/

/*************************** Dynamically allocate organelle matrix which stores the organelle info **************************************/
 /*double *organelle = (double *)malloc(npart*sizeof(double));
 for(i=0;i<npart;i++)
 organelle[i]=0; // initialise the organelle array to 0 i.e., no organelles are present

 organelle[1]=1;*/ // this particle is an organelle
 /*organelle[7]=1;
 organelle[12]=1;
 organelle[15]=1;*/
 //organelle[5]=1;
 //organelle[6]=1;
 //organelle[11]=1;

/******************************** store coordinates of organelles into organelle_pos matrix ***********************************/
 /*for(i=0;i<npart;i++)
 {
  if(organelle[i]==1)
  {
   organelle_pos[i][0]=position[i][0];
   organelle_pos[i][1]=position[i][1];
   organelle_pos[i][2]=position[i][2];
   organelle_pos[i][3]=position[i][3];
   organelle_pos[i][4]=position[i][4];
  }
 }*/

/*************************** Generate the force matrix *******************************************************************************/
 
 force=vanderwalforce(force,position,mass,npart); // to calculate vanderwal force
 printf("\n\n FORCE MATRIX ----------------\n");
 for(i=0;i<sumpart;i++)
 {
  for(j=0;j<5;j++)
  printf("%f ",force[i][j]);
  printf("\n");
 }

/************************* Allocate arr_startx,arr_starty,arr_startz,arr_width,arr_height,arr_depth *****************************************************/
 i=0;boxcount=0;startx=0;width=widthc;starty=0;height=heightc;startz=0;depth=depthc;
 while(i<row)
 {
  colcount=0;
  while(colcount<col)
  {
   arr_startx[boxcount]=startx;
   arr_starty[boxcount]=starty;
   arr_startz[boxcount]=startz;
   arr_width[boxcount]=width;
   arr_height[boxcount]=height;
   arr_depth[boxcount]=depth;
   boxcount++;
   colcount++;
   startx=width+1;width=startx+widthc;
  }
  i++;startx=0;width=widthc;starty=height+1;height=starty+heightc; 
 }

 for(i=0;i<N;i++) // making transitions between boxes continuous
 {
  if((i%col)==0)
  continue;
  else
  arr_startx[i]-=0.5;
 }
  
 for(i=0;i<N;i++)
 {
  if((i%col)==(col-1))
  continue;
  else
  arr_width[i]+=0.5;
 }

  
 for(i=col;i<N;i++)
 arr_starty[i]-=0.5;

 for(i=0;i<N-col;i++)
 arr_height[i]+=0.5;
  

 printf("Printing arr_startx-------------------------\n");
 for(i=0;i<N;i++)
 printf("%f ",arr_startx[i]);
 printf("\n");
 printf("Printing arr_starty-------------------------\n");
 for(i=0;i<N;i++)
 printf("%f ",arr_starty[i]);
 printf("\n");
 printf("Printing arr_startz-------------------------\n");
 for(i=0;i<N;i++)
 printf("%f ",arr_startz[i]);
 printf("\n");
 printf("Printing arr_width-------------------------\n");
 for(i=0;i<N;i++)
 printf("%f ",arr_width[i]);
 printf("\n");
 printf("Printing arr_height-------------------------\n");
 for(i=0;i<N;i++)
 printf("%f ",arr_height[i]);
 printf("\n");
 printf("Printing arr_depth-------------------------\n");
 for(i=0;i<N;i++)
 printf("%f ",arr_depth[i]);
 printf("\n");

 











/************************* Convert main 2D arrays to 1D arrays to ease cuda processing ******************************************************************/
 
 double *scale_halfvelocity;// 1D scaled version of halfvelocity        host copies
 double *scale_velocity;// 1D scaled version of velocity
 double *scale_force;// 1D scaled version of force
 double *scale_position;// 1D scaled version of position
 double *dup_scale_position,*dup_scale_halfvelocity,*dup_scale_velocity,*dup_scale_force;

 scale_halfvelocity=(double *)malloc(5*sumpart*sizeof(double));
 scale_velocity=(double *)malloc(5*sumpart*sizeof(double));
 scale_force=(double *)malloc(5*sumpart*sizeof(double));
 scale_position=(double *)malloc(5*sumpart*sizeof(double));
 dup_scale_position=(double *)malloc(5*sumpart*sizeof(double));
 dup_scale_halfvelocity=(double *)malloc(5*sumpart*sizeof(double));
 dup_scale_velocity=(double *)malloc(5*sumpart*sizeof(double));
 dup_scale_force=(double *)malloc(5*sumpart*sizeof(double));



 
 count=0;
 for(i=0;i<sumpart;i++) // Allocating halfvelocity matrix to scale_halfvelocity array
 {
  for(j=0;j<5;j++)
  {
   scale_halfvelocity[count]=halfvelocity[i][j];
   count++;
  }
 }


 
 count=0;
 for(i=0;i<sumpart;i++) // Allocating velocity matrix to scale_velocity array
 {
  for(j=0;j<5;j++)
  {
   scale_velocity[count]=velocity[i][j];
   count++;
  }
 }

 count=0;
 for(i=0;i<sumpart;i++) // Allocating force matrix to scale_force array
 {
  for(j=0;j<5;j++)
  {
   scale_force[count]=force[i][j];
   count++;
  }
 }

 count=0;
 for(i=0;i<sumpart;i++) // Allocating position matrix to scale_position array
 {
  for(j=0;j<5;j++)
  {
   scale_position[count]=position[i][j];
   count++;
  }
 }

 count=0;
 for(i=0;i<sumpart;i++) // Initialising dup_scale_position 
 {
  for(j=0;j<5;j++)
  {
   dup_scale_position[count]=0;
   count++;
  }
 }

 count=0;
 for(i=0;i<sumpart;i++) // Initialising dup_scale_halfvelocity 
 {
  for(j=0;j<5;j++)
  {
   dup_scale_halfvelocity[count]=0;
   count++;
  }
 }

 count=0;
 for(i=0;i<sumpart;i++) // Initialising dup_scale_velocity 
 {
  for(j=0;j<5;j++)
  {
   dup_scale_velocity[count]=0;
   count++;
  }
 }

 count=0;
 for(i=0;i<sumpart;i++) // Initialising dup_scale_force 
 {
  for(j=0;j<5;j++)
  {
   dup_scale_force[count]=0;
   count++;
  }
 }

 free(halfvelocity);free(velocity);free(force);free(position);


/************************* Device initialisations and allocating space **********************************************************************************/
 //cudaDeviceReset();
 double *dev_halfvelocity,*dev_velocity,*dev_force,*dev_position,*dev_mass,*dev_radius;// device copies
 int *dev_arr_startrow,*dev_arr_lastrow;
 double *dev_arr_startx,*dev_arr_starty,*dev_arr_startz,*dev_arr_width,*dev_arr_height,*dev_arr_depth;
 double *dev_x_interval,*dev_y_interval;

 cudaMalloc((void **) &dev_halfvelocity,5*sumpart*sizeof(double));
 cudaMalloc((void **) &dev_velocity,5*sumpart*sizeof(double));
 cudaMalloc((void **) &dev_force,5*sumpart*sizeof(double));
 cudaMalloc((void **) &dev_position,5*sumpart*sizeof(double));
 cudaMalloc((void **) &dev_mass,1*sumpart*sizeof(double));
 cudaMalloc((void **) &dev_radius,1*sumpart*sizeof(double));
 cudaMalloc((void **) &dev_arr_startrow,N*sizeof(int));
 cudaMalloc((void **) &dev_arr_lastrow,N*sizeof(int));
 cudaMalloc((void **) &dev_arr_startx,N*sizeof(double));
 cudaMalloc((void **) &dev_arr_starty,N*sizeof(double));
 cudaMalloc((void **) &dev_arr_startz,N*sizeof(double));
 cudaMalloc((void **) &dev_arr_width,N*sizeof(double));
 cudaMalloc((void **) &dev_arr_height,N*sizeof(double));
 cudaMalloc((void **) &dev_arr_depth,N*sizeof(double));
 cudaMalloc((void **) &dev_x_interval,1*4*sizeof(double));
 cudaMalloc((void **) &dev_y_interval,1*4*sizeof(double));



 cudaMemcpy(dev_halfvelocity,scale_halfvelocity,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_velocity,scale_velocity,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_force,scale_force,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_position,scale_position,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_mass,mass,1*sumpart*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_radius,radius,1*sumpart*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_startrow,arr_startrow,N*sizeof(int),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_lastrow,arr_lastrow,N*sizeof(int),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_startx,arr_startx,N*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_starty,arr_starty,N*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_startz,arr_startz,N*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_width,arr_width,N*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_height,arr_height,N*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_arr_depth,arr_depth,N*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_x_interval,x_interval,1*4*sizeof(double),cudaMemcpyHostToDevice);
 cudaMemcpy(dev_y_interval,y_interval,1*4*sizeof(double),cudaMemcpyHostToDevice);


 printf("\nInitial Force-------------------\n");
 for(i=0;i<5*sumpart;i+=5)
 printf("%d %d %f %f %f\n",(int)scale_force[i],(int)scale_force[i+1],scale_force[i+2],scale_force[i+3],scale_force[i+4]);


 for(t=0;t<total_time;t++)
 {
  velocityverlet<<<N,1>>>(dev_x_interval,dev_y_interval,dev_halfvelocity,dev_velocity,dev_force,dev_position,dev_mass,dev_radius,dev_arr_startrow,dev_arr_lastrow,npart,dev_arr_startx,dev_arr_starty,dev_arr_startz,dev_arr_width,dev_arr_height,dev_arr_depth);
  //gpuErrchk( cudaPeekAtLastError() );
  //gpuErrchk( cudaDeviceSynchronize() );
  cudaDeviceSynchronize();
  cudaThreadSynchronize();

  cudaMemcpy(scale_position,dev_position,sumpart*5*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(scale_force,dev_force,sumpart*5*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(scale_velocity,dev_velocity,sumpart*5*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(scale_halfvelocity,dev_halfvelocity,sumpart*5*sizeof(double),cudaMemcpyDeviceToHost);

  printf("\n%d Force be4 exchange -------------------\n",t);
  for(i=0;i<5*sumpart;i+=5)
  printf("%d %d %f %f %f\n",(int)scale_force[i],(int)scale_force[i+1],scale_force[i+2],scale_force[i+3],scale_force[i+4]);
 
  printf("\n%d Position be4 exchange -------------------\n",t);
  for(i=0;i<5*sumpart;i+=5)
  printf("%d %d %f %f %f\n",(int)scale_position[i],(int)scale_position[i+1],scale_position[i+2],scale_position[i+3],scale_position[i+4]);

  printf("\n%d Velocity be4 exchange -------------------\n",t);
  for(i=0;i<5*sumpart;i+=5)
  printf("%d %d %f %f %f\n",(int)scale_velocity[i],(int)scale_velocity[i+1],scale_velocity[i+2],scale_velocity[i+3],scale_velocity[i+4]);

  printf("\n%d HalfVelocity be4 exchange -------------------\n",t);
  for(i=0;i<5*sumpart;i+=5)
  printf("%d %d %f %f %f\n",(int)scale_halfvelocity[i],(int)scale_halfvelocity[i+1],scale_halfvelocity[i+2],scale_halfvelocity[i+3],scale_halfvelocity[i+4]);

  printf("\n%d Arr_startrow bef4 exchange--------------------\n",t);
  for(i=0;i<N;i++)
  printf("%d ",arr_startrow[i]);
 
  printf("\n%d Arr_lastrow bef4 exchange--------------------\n",t);
  for(i=0;i<N;i++)
  printf("%d ",arr_lastrow[i]);
  
  /*for(i=0;i<npart;i++) // keep organelle posiiton constant and organelle velocity to 0
  {
   if(organelle[i]==1)
   {
    for(j=0;j<npart*5;j+=5)
    {
     if(((int)scale_position[j+1])==i)
     {
      scale_position[j+0]=organelle_pos[i][0];
      //scale_position[(j*5)+1]=organelle_pos[i][1];
      scale_position[j+2]=organelle_pos[i][2];
      scale_position[j+3]=organelle_pos[i][3];
      scale_position[j+4]=organelle_pos[i][4];
      scale_halfvelocity[j+0]=organelle_pos[i][0];
      scale_halfvelocity[j+2]=0;scale_halfvelocity[j+3]=0;scale_halfvelocity[j+4]=0;
      scale_velocity[j+0]=organelle_pos[i][0];
      scale_velocity[j+2]=0;scale_velocity[j+3]=0;scale_velocity[j+4]=0;
     }
    }
   }
  }*/

/******* rearranging scale_position according to core nos i.e., using 1st column of scale_position *********************/
  count=0;count_core=0;
  for(rearr=0;rearr<N;rearr++)
  {
   for(j=0;j<5*npart;j+=5)
   {
    if(scale_position[j]==rearr)
    {
     count_core++;
     dup_scale_position[count]=scale_position[j];
     dup_scale_force[count]=scale_force[j];
     dup_scale_velocity[count]=scale_velocity[j];
     dup_scale_halfvelocity[count]=scale_halfvelocity[j];
     count++;
     dup_scale_position[count]=scale_position[j+1];
     dup_scale_force[count]=scale_force[j+1];
     dup_scale_velocity[count]=scale_velocity[j+1];
     dup_scale_halfvelocity[count]=scale_halfvelocity[j+1];
     count++;
     dup_scale_position[count]=scale_position[j+2];
     dup_scale_force[count]=scale_force[j+2];
     dup_scale_velocity[count]=scale_velocity[j+2];
     dup_scale_halfvelocity[count]=scale_halfvelocity[j+2];
     count++;
     dup_scale_position[count]=scale_position[j+3];
     dup_scale_force[count]=scale_force[j+3];
     dup_scale_velocity[count]=scale_velocity[j+3];
     dup_scale_halfvelocity[count]=scale_halfvelocity[j+3];
     count++;
     dup_scale_position[count]=scale_position[j+4];
     dup_scale_force[count]=scale_force[j+4];
     dup_scale_velocity[count]=scale_velocity[j+4];
     dup_scale_halfvelocity[count]=scale_halfvelocity[j+4];
     count++;
    }
   }
   arr_lastrow[rearr]=count_core;
  }
  
  arr_startrow[0]=0;
  for(rearr=1;rearr<N;rearr++)
  {
   arr_startrow[rearr]=arr_lastrow[rearr-1];
  }



  fprintf(fp,"%d\n",npart);
  for(j=0;j<5*npart;j+=5)
  {
    fprintf(fp,"Atom%d ",(int)dup_scale_position[j+1]);
    fprintf(fp,"%f %f %f\n",dup_scale_position[j+2],dup_scale_position[j+3],dup_scale_position[j+4]);
   
  }

  printf("\n%d Force after exchange -------------------\n",t);
  for(i=0;i<5*sumpart;i+=5)
  printf("%d %d %f %f %f\n",(int)dup_scale_force[i],(int)dup_scale_force[i+1],dup_scale_force[i+2],dup_scale_force[i+3],dup_scale_force[i+4]);

  printf("\n%d Position after exchange -------------------\n",t);
  for(i=0;i<5*sumpart;i+=5)
  printf("%d %d %f %f %f\n",(int)dup_scale_position[i],(int)dup_scale_position[i+1],dup_scale_position[i+2],dup_scale_position[i+3],dup_scale_position[i+4]);

  printf("\n%d Arr_startrow after exchange--------------------\n",t);
  for(i=0;i<N;i++)
  printf("%d ",arr_startrow[i]);

  printf("\n%d Arr_lastrow after exchange--------------------\n",t);
  for(i=0;i<N;i++)
  printf("%d ",arr_lastrow[i]);



  cudaMemcpy(dev_position,dup_scale_position,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(dev_arr_startrow,arr_startrow,N*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(dev_arr_lastrow,arr_lastrow,N*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(dev_halfvelocity,dup_scale_halfvelocity,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(dev_velocity,dup_scale_velocity,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(dev_force,dup_scale_force,5*sumpart*sizeof(double),cudaMemcpyHostToDevice);

 }

 fclose(fp);
 free(particle);free(mass);free(radius);
 free(scale_force);free(scale_halfvelocity);
 free(scale_velocity);free(scale_position);
 cudaFree(dev_halfvelocity);cudaFree(dev_velocity);cudaFree(dev_force);cudaFree(dev_position);cudaFree(dev_mass);cudaFree(dev_radius);
 cudaFree(dev_arr_startrow);cudaFree(dev_arr_lastrow);cudaFree(dev_arr_startx);cudaFree(dev_arr_starty);cudaFree(dev_arr_startz);
 cudaFree(dev_arr_width);cudaFree(dev_arr_height);cudaFree(dev_arr_depth);
 cudaDeviceReset();
 return 0;

}

 
 
