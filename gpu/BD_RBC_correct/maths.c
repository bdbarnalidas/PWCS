__host__ __device__ double DISTV(double s,double *u,double *v)
{

 double tmp;
 int i;                                                                
 tmp = 0.0;                                                            
 for (i = 0; i < NDIM; i++)                                          
 tmp += (u[i]-v[i]) *(u[i]-v[i]);                         
 
  s = sqrt(tmp); 

 return s;
}
 
__host__ __device__ double* SUBV(double *v, double *u,double *w)
{
 int i;                                                             
 for (i = 0; i < NDIM; i++)                                       
 v[i] = u[i] - w[i];
 return v;
}

__host__ __device__ double* MULVS(double *v,double *u,double s)
{
 int i;                                                             
 for (i = 0; i < NDIM; i++)                                       
 v[i] = u[i] * s;         
 return v;                               
}
