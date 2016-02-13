#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

//gcc -O5 -DPOISSON berk.c -lm -o berk.exe;gcc -O5 -DAMPERE -DJ0 berk.c -lm -o berka0.exe;
//gcc -O5 -DAMPERE -DJ1 berk.c -lm -o berka1.exe;gcc -O5 -DAMPERE -DJ2 berk.c -lm -o berka2.exe

//time ./berk.exe 128 128 8000 0.1 1 0.03 0.3 9;mv diag.plot diagp.plot
//time ./berka0.exe 128 128 8000 0.1 1 0.03 0.3 9;mv diag.plot diagpa0.plot
//time ./berka1.exe 128 128 8000 0.1 1 0.03 0.3 9;mv diag.plot diagpa1.plot
//time ./berka2.exe 128 128 8000 0.1 1 0.03 0.3 9;mv diag.plot diagpa2.plot

//time ./berka0.exe 64 256 10000 0.1 1 0.01 0.3 8 1 0.03162 -1;mv diag.plot diaga0.plot
//time ./berka1.exe 64 256 10000 0.1 1 0.01 0.3 8 1 0.03162 -1;mv diag.plot diaga1.plot
//time ./berka2.exe 64 256 10000 0.1 1 0.01 0.3 8 1 0.03162 -1;mv diag.plot diaga2.plot

#define Nbdr 20
#define Nval 2

#define domf(){\
domlandau();\
}

#define f(){\
fberk();\
}

#define domlandau(){\
xmin=0.;xmax=2*M_PI/landau_k;vmin=-vmax;\
}

#define fberk(){\
result= 1./sqrt(2.*M_PI)*(0.9*exp(-.5*v*v)+0.2*exp(-0.5*(v-4.5)*(v-4.5)/(0.5*0.5)))*(1.+landau_alpha*cos(landau_k*x));\
}

#define fberk0(){\
result= 1./sqrt(2.*M_PI)*(0.9*exp(-.5*v*v)+0.2*exp(-0.5*(v-4.5)*(v-4.5)/(0.5*0.5)));\
}

int interpol_per0(double *dtab,double *ltab,double *mtab,int N){
  //init periodic splines
  int i;
  dtab[0] = 4.;
  mtab[0] = 0.25;
  for(i=0;i<N-1;i++){
    ltab[i]=1./dtab[i];
    dtab[i+1]=4.-ltab[i];
    mtab[i+1]=-mtab[i]/dtab[i+1];
  }
  dtab[N-1]-=ltab[N-2] + 2.0*mtab[N-2];
  for(i=0;i<N;i++)dtab[i]=1./dtab[i];
  return 0;  
}

int interpol_nat0(double *dnat,double *lnat,int N){
  //init natural splines
  int i;
  dnat[0]=4.;
  for(i=0;i<N-2;i++){
    lnat[i]=1./dnat[i];
    dnat[i+1]=4.-lnat[i];
  }
  return 0;
}

int interpol_per(double* bufin, double* bufout, int N, double alpha,double* dtab, double* ltab, double* mtab)
{
	double * ftab;
	double x, w0, w1, w2, w3;
	double uns6=1.0/6.0;
	int i, ix, ix0; 
	
	ftab=&bufin[Nbdr];
  /*resolution of L phi = ftab*/
  for(i=1;i<N;i++) 
		ftab[i]-=ftab[i-1]*ltab[i-1];
  
  for(i=0;i<N-1;i++)
		ftab[N-1]-=mtab[i]*ftab[i];
  /*resolution of U eta =phi*/
  for(i=0;i<N;i++)
		ftab[i]*=6.;
  
  ftab[N-1]*=dtab[N-1];
  ftab[N-2]=dtab[N-2]*(ftab[N-2]-(1.-mtab[N-3])*ftab[N-1]);
  for(i=N-3;i>=1;i--) ftab[i]=dtab[i]*(ftab[i]-ftab[i+1]+mtab[i-1]*ftab[N-1]);
  ftab[0]=dtab[0]*(ftab[0]-ftab[1]-ftab[N-1]);


  for(i=1;i<=Nbdr;i++)bufin[Nbdr-i]=bufin[N-i+Nbdr];
  for(i=0;i<Nbdr;i++)bufin[N+i+Nbdr]=bufin[i+Nbdr];

  x=-alpha;
  if(x<0) x+=1;
  if(x>=1)x-=1;
  ix = (int) (x*N);
  assert(ix>=0 &&ix<N);
  x = x*N-ix;
  assert(x>=0 &&x<1);
  ix0=ix;
  //if(alpha==0)if(ix!=i){fprintf(stderr,"Pb\n");return 1;}
  //if(alpha==0)if(xloc!=0){fprintf(stderr,"Pb\n");return 1;}
  
  w0=uns6*(1.0-x)*(1.0-x)*(1.0-x);
  w1=uns6+0.5*(1.0-x)*(-(1.0-x)*(1.0-x)+(1.0-x)+1.0);
  w2=uns6+0.5*x*(-x*x+x+1.0);
  w3=uns6*x*x*x;

  ix=ix0-1;
  for(i=0;i<N;i++){
    ix++;if(ix>=N)ix-=N;
    /*interpolation*/
    bufout[i]=w0*ftab[ix-1]+w1*ftab[ix]+w2*ftab[ix+1]+w3*ftab[ix+2];
  }
  return 0;
}


int interpol_nat(double* bufin, double* p, int N, double alpha,double* dnat,double* lnat)
{
  //warning: f_0=0...f_{N}=0 N intervals and N+1 points
  int j,ix;
  double *coef;
  double x,dx=1./N;
  double w[4];
  coef=&bufin[Nbdr];
  for(j=0;j<N-1;j++)p[j]=6.*coef[j+1];p[0]/=dnat[0];for(j=1;j<N-1;j++)p[j]=(p[j]-p[j-1])/dnat[j];
  for(j=N-3;j>=0;j--) p[j]-=p[j+1]*lnat[j];
   //from 0 eta_1 eta_2 ... eta_{N -1} add eta_{-1} = -eta_1 eta_N=0 and eta_{N+1}=-eta_{N-1}
  coef[0]=-p[0];coef[1]=0.;for(j=0;j<N-1;j++)coef[j+2]=p[j];coef[N+1]=0.;coef[N+2]=-coef[N];
  p[0]=0.;
  for(j=1;j<N;j++){
    x=j*dx-alpha;
    if(!(x>=0.&&x<1.))p[j]=0;else{
      ix=(int)(x*N);
      assert(x>=0.&&x<1.);assert(ix>=0 &&ix<N);x=x*N-ix;assert(x>=0 &&x<1.);
      w[0]=(1./6.)*(1-x)*(1-x)*(1-x);w[1]=1./6.+0.5*(1-x)*(-(1-x)*(1-x)+(1-x)+1);w[2]=1./6.+0.5*x*(-x*x+x+1);w[3]=(1./6.)*x*x*x;      
      p[j]=w[0]*coef[ix]+w[1]*coef[ix+1]+w[2]*coef[ix+2]+w[3]*coef[ix+3];
    }
  }
  return 0;
}


 
int main(int argc,char **argv){
  
  //parameters
  int nbarg=1;
  int Nx=128,Nv=128,interp=1,nbiter=100,ndiag=-1;
  double dt=0.1,landau_alpha=0.01,landau_k=0.5,vmax=9.,gamma=0.,nu=0.;
  
  int i,j,s,N,date;
  
  double xmin,xmax,vmin,x,v,dx,dv,result,dtmp;
  double alpha,ee,emax;
  
  double *dper,*lper,*mper,*dnat,*lnat,*dperv,*lperv,*mperv;
  double *bufin,*bufout;
  double *ftab,*E,*rho,*jrho;
  double *diag;
  double *rho1;
  
  FILE *file;
  char str[256];
  
  if(argc>nbarg)Nx=atoi(argv[nbarg++]);
  if(argc>nbarg)Nv=atoi(argv[nbarg++]);
  if(argc>nbarg)nbiter=atoi(argv[nbarg++]);
  if(argc>nbarg)dt=atof(argv[nbarg++]);
  if(argc>nbarg)interp=atoi(argv[nbarg++]);
  if(argc>nbarg)landau_alpha=atof(argv[nbarg++]);  
  if(argc>nbarg){landau_k=atof(argv[nbarg++]);}
  if(argc>nbarg){vmax=atof(argv[nbarg++]);}
  if(argc>nbarg){gamma=atof(argv[nbarg++]);}
  if(argc>nbarg){nu=atof(argv[nbarg++]);}
  if(argc>nbarg)ndiag=atoi(argv[nbarg++]);
  
  fprintf(stderr,"Nx=%d Nv=%d nbiter=%d dt=%lg interp=%d\n",Nx,Nv,nbiter,dt,interp);
  fprintf(stderr,"landau_alpha=%lg landau_k=%lg vmax=%lg gamma=%lg nu=%lg ndiag=%d\n",landau_alpha,landau_k,vmax,gamma,nu,ndiag);
  
  N=Nx;if(N<Nv)N=Nv;
  
  domf();

  dx = (xmax-xmin)/(double) Nx;
  dv = (vmax-vmin)/(double) Nv;

  ftab=(double*)malloc(sizeof(double)*Nx*Nv);

  E=(double*)malloc(sizeof(double)*Nx);
  rho=(double*)malloc(sizeof(double)*Nx);
  jrho=(double*)malloc(sizeof(double)*Nx);
  rho1=(double*)malloc(sizeof(double)*Nx);
  
  dper=(double*)malloc(sizeof(double)*Nx);
  lper=(double*)malloc(sizeof(double)*(Nx-1));
  mper=(double*)malloc(sizeof(double)*Nx);
  #ifdef PER
  dperv=(double*)malloc(sizeof(double)*Nv);
  lperv=(double*)malloc(sizeof(double)*(Nv-1));
  mperv=(double*)malloc(sizeof(double)*Nv);
  #else
  dnat=(double*)malloc(sizeof(double)*(Nv-1));
  lnat=(double*)malloc(sizeof(double)*(Nv-2));
  #endif
  bufin=(double*)malloc(sizeof(double)*(N+2*Nbdr+1));
  bufout=(double*)malloc(sizeof(double)*N);
  
  diag=(double*)malloc(sizeof(double)*Nval*(nbiter+1));

  dtmp=0.;for(j=0;j<Nv;j++){v=vmin+(double)j*dv;fberk0();dtmp+=result*dv;}  
  fprintf(stderr,"rho0=%1.20lg %1.20lg\n",dtmp,1-dtmp);
  //exit(1);
  //init of distribution function
  for(j=0;j<Nv;j++)for(i=0;i<Nx;i++){
    x=xmin+(double)i*dx;v=vmin+(double)j*dv;
    f();ftab[i+Nx*j]=result;
  }
  
  interpol_per0(dper,lper,mper,Nx);  
  #ifdef PER
  interpol_per0(dperv,lperv,mperv,Nv);
  #else
  interpol_nat0(dnat,lnat,Nv);
  #endif
  for(i=0;i<Nx;i++)rho[i]=0.;
  for(j=0;j<Nv;j++)for(i=0;i<Nx;i++)rho[i]+=ftab[i+Nx*j]*dv;   

  E[0]=0.;for(i=1;i<Nx;i++) E[i]=E[i-1]+0.5*dx*(rho[i-1]+rho[i]-2.);
  dtmp=0.;for(i=1;i<Nx;i++)dtmp+=E[i];dtmp=-dtmp/(double)Nx;
  for(i=0;i<Nx;i++)E[i]+=dtmp;
  
  FILE * pFile;
  pFile = fopen ("df.dat","w");
  if (pFile!=NULL)
  {
    for(j=0;j<Nv;j++) {
      v=vmin+j*dv;
      for(i=0;i<Nx;i++) {
        x=xmin+i*dx;
        f();
        fprintf(pFile, "%f %f %f \n", x, v, result);
      }
      fprintf(pFile, "\n");
    }
    fclose (pFile);
  }

  pFile = fopen ("F0.dat","w");
  if (pFile!=NULL)
  {
    for(j=0;j<Nv;j++) {
      v=vmin+j*dv;
      fberk0();
      fprintf(pFile, "%f %f \n", v, result);
    }
    fclose (pFile);
  }

  pFile = fopen ("E0.dat","w");
  if (pFile!=NULL)
  {
    for(i=0;i<Nx;i++) {
      x=xmin+i*dx;
      fprintf(pFile, "%f %f %f \n", x, E[i], rho[i]);
    }
    fclose (pFile);
  }


  ee=0;for(i=0;i<Nx;i++)ee+=E[i]*E[i];ee*=1./(double)Nx;//ee*=dx;
  printf("\n ee = %15.12f \n ", ee);


  for(date=1;date<=nbiter;date++){  
    
    #ifdef J0
    /*compute J directly*/
    for(i=0;i<Nx;i++)jrho[i]=0.;
    for(j=0;j<Nv;j++){v=vmin+(double)j*dv;for(i=0;i<Nx;i++){jrho[i]+=ftab[i+Nx*j]*v*dv;}}
    dtmp=0.;for(i=0;i<Nx;i++)dtmp+=jrho[i];dtmp*=(1./(double)Nx);for(i=0;i<Nx;i++)jrho[i]-=dtmp;    
    #endif
    
    /*advection of a 1/2 time step in x*/
    for(j=0;j<Nv;j++){
      alpha=-0.5*dt*(vmin+(double)j*dv)/(xmax-xmin);
      for(i=0;i<Nx;i++)bufin[i+Nbdr]=ftab[i+Nx*j];
      interpol_per(bufin,ftab+Nx*j,Nx,alpha,dper,lper,mper);
    }
    
    /*collision operator 1/2 time step*/
    //\partial_tf=-nu (f-F0)
    dtmp=exp(-nu*0.5*dt);
    for(j=0;j<Nv;j++){
      v=vmin+j*dv;fberk0();result*=(1-dtmp);
      for(i=0;i<Nx;i++)ftab[i+Nx*j]=ftab[i+Nx*j]*dtmp+result;
    }
    
    
    #ifdef J1
    /*compute J from rho*/
    for(i=0;i<Nx;i++)rho1[i]=0.;
    for(j=0;j<Nv;j++)for(i=0;i<Nx;i++)rho1[i]+=ftab[i+Nx*j]*dv;       
    for(i=0;i<Nx;i++){dtmp=rho1[i];rho1[i]=(rho1[i]-rho[i])/dt+0.5*nu*(rho1[i]+rho[i]-2.);rho[i]=dtmp;}
    jrho[0]=0.;for(i=1;i<Nx;i++) jrho[i]=jrho[i-1]+0.5*dx*(rho1[i-1]+rho1[i]);
    dtmp=0.;for(i=1;i<Nx;i++)dtmp+=jrho[i];dtmp=-dtmp/(double)Nx;
    for(i=0;i<Nx;i++)jrho[i]+=dtmp;
    #endif
    
    #ifdef J2
    /*compute J like in Vann_PhdThesis03  WARNING: rho1 should be called jrho1*/
    for(i=0;i<Nx;i++)jrho[i]=0.;
    for(j=0;j<Nv;j++){v=vmin+(double)j*dv;for(i=0;i<Nx;i++){jrho[i]+=ftab[i+Nx*j]*v*dv;}}
    dtmp=0.;for(i=0;i<Nx;i++)dtmp+=jrho[i];dtmp=-dtmp/(double)Nx;
    for(i=0;i<Nx;i++)jrho[i]+=dtmp;
    for(i=0;i<Nx;i++)jrho[i]=0.5*(rho1[i]+jrho[i]);    
    #endif
    
    #ifdef POISSON
    for(i=0;i<Nx;i++)rho[i]=0.;
    for(j=0;j<Nv;j++)for(i=0;i<Nx;i++)rho[i]+=ftab[i+Nx*j]*dv;       
    E[0]=0.;for(i=1;i<Nx;i++) E[i]=E[i-1]+0.5*dx*(rho[i-1]+rho[i]-2.);
    dtmp=0.;for(i=1;i<Nx;i++)dtmp+=E[i];dtmp=-dtmp/(double)Nx;
    for(i=0;i<Nx;i++)E[i]+=dtmp;
    #endif
    
    #ifdef AMPERE
    if(date<2){
      for(i=0;i<Nx;i++)rho[i]=0.;
      for(j=0;j<Nv;j++)for(i=0;i<Nx;i++)rho[i]+=ftab[i+Nx*j]*dv;       
      E[0]=0.;for(i=1;i<Nx;i++) E[i]=E[i-1]+0.5*dx*(rho[i-1]+rho[i]-2.);
      dtmp=0.;for(i=1;i<Nx;i++)dtmp+=E[i];dtmp=-dtmp/(double)Nx;
      for(i=0;i<Nx;i++)E[i]+=dtmp;
    }
    else for(i=0;i<Nx;i++){E[i]=((1.-0.5*gamma*dt)/(1.+0.5*gamma*dt))*E[i]+(dt/(1.+0.5*gamma*dt))*jrho[i];}
    #endif
    
    emax=0.;for(i=0;i<Nx;i++)if(fabs(E[i])>emax)emax=fabs(E[i]);
    ee=0;for(i=0;i<Nx;i++)ee+=(E[i]*E[i]);ee*=1./(double)Nx;//ee*=dx;
    diag[Nval*date]=ee;diag[Nval*date+1]=emax;

    /* advection of a 1 time step in v*/
    for(i=0;i<Nx;i++){
      alpha=-E[i]*dt/(vmax-vmin);
      for(j=0;j<Nv;j++)bufin[j+Nbdr]=ftab[i+Nx*j];
      #ifdef PER
      interpol_per(bufin,bufout,Nv,alpha,dperv,lperv,mperv);
      #else
      interpol_nat(bufin,bufout,Nv,alpha,dnat,lnat);      
      #endif
      for(j=0;j<Nv;j++)ftab[i+Nx*j]=bufout[j];    
    }

    #ifdef J2
    /*compute J like in Vann_PhdThesis03 WARNING: rho1 should be jrho1*/
    for(i=0;i<Nx;i++)rho1[i]=0.;
    for(j=0;j<Nv;j++){v=vmin+(double)j*dv;for(i=0;i<Nx;i++){rho1[i]+=ftab[i+Nx*j]*v*dv;}}
    dtmp=0.;for(i=0;i<Nx;i++)dtmp+=rho1[i];dtmp=-dtmp/(double)Nx;
    for(i=0;i<Nx;i++)rho1[i]+=dtmp;
    #endif

    /*collision operator 1/2 time step*/
    //\partial_tf=-nu (f-F0)
    dtmp=exp(-nu*0.5*dt);
    for(j=0;j<Nv;j++){
      v=vmin+j*dv;fberk0();result*=(1-dtmp);
      for(i=0;i<Nx;i++)ftab[i+Nx*j]=ftab[i+Nx*j]*dtmp+result;
    }

    /* advection of a 1/2 time step in x*/
    for(j=0;j<Nv;j++){
      alpha=-0.5*dt*(vmin+(double)j*dv)/(xmax-xmin);
      for(i=0;i<Nx;i++)bufin[i+Nbdr]=ftab[i+Nx*j];
      interpol_per(bufin,ftab+Nx*j,Nx,alpha,dper,lper,mper);
    }
  }
  
  file=fopen("diag.plot","w");
  for(date=1;date<=nbiter;date++){
    fprintf(file,"%lg",date*dt);
    for(i=0;i<Nval;i++)fprintf(file," %1.20lg",diag[Nval*date+i]);
    fprintf(file,"\n");
  }
  fclose(file);

  
  return 0;
}
