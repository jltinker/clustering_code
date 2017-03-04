#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include "nrutil.h"

#define c_on_H0 3000
#define PI 3.14159
#define OMEGA_M 0.29
#define RA_MIN 100
#define RA_MAX 260
#define Q0 2.0
#define Q1 -1.0
#define QZ0 0.1
#define COLLISION_ANGLE (62./3600.*3.14159/180.)

#define ind(i,j) (i*njack+j)
#define mabs(A) ((A) < 0.0 ? -(A) : (A))

float SPEED_OF_LIGHT  = 1; 
float COLLISION_WEIGHT = 2.64;
float COLLISION_WEIGHT2 = 1;

//internal functions
double redshift_distance(double z);
double func_dr1(double z);
float random_radial_position(float zlo, float zhi);
float angular_separation(float a1, float d1, float a2, float d2);
float random_redshift_from_data(int n, float *z);
float collision_weight(float theta, int iflag);

// external functions
double qromo(double (*func)(double), double a, double b,
	     double (*choose)(double(*)(double), double, double, int));
double midpnt(double (*func)(double), double a, double b, int n);
void sort2(int n, float arr[], int id[]);

void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart, int i0);

double second(void);
double timediff(double t0,double t1);

int main(int argc, char **argv)
{
  FILE *fp, *fpmag, *fpkcorr;
  char aa[1000];

  float rmin, rmax, dlogr, r, dx, dy, dz, xt, rcube = 125, BUF = 20, 
    x1, r1, y1, z1, zz1, ra1, x11, kcorr, deltaz;
  int nr = 14;
  int ngal,ngal2,i,j,k,nmock,ibin,**nrand,nrandoms,**npairs_dr2, ngal_tot;
  float *x,*y,*z;
  float *x2,*y2,*z2;
  float *rx1,*rx2,*ry1,*ry2,*rz1,*rz2;
  float *xg, *yg;
  int nrand1, nrand2;
  double npairs_tot;
  double *rbar,nratio,DR1,DR2,RR,DD,xi,**pibar;
  double **npairs, **npairs_dr1;

  float maglo=-10, maghi=-30, zlo=0, zhi=1;
  float zmin, zmax, theta;
  int nz, jbin, ntmp;
  int CALC_RR_COUNTS=0;
  int input_code = 2;

  float lx,ly,lz,pi,phi, PI_MAX, weight, *gweight;

  //variables for jackknife calculations
  int ***nrand_jack, **nx, ***nxj, njack1, nby5,nsby5,ii, *subindx, ijack, n;
  double ***npairs_dr1_jack, ***npairs_jack;
  int ijack1, ijack2, jjack, njack_ra, njack_dec, njack_tot, njack=10;
  float dx_jack=15, dy_jack=0.1, ra_min=100, dec_min=-0.07;
  float *gal_ra, *gal_dec, *gal_dec2;
  float *ran_ra, *ran_dec;
  int id, idp, id2;
  int *rancnt_jack;
  float *galcnt_jack;
  float xi_err, xi_jack, xi_bar;
  double ngal_tmp, nrand_tmp;
  float xi_data[30], xij[200*200][30], covar[30][30];
  FILE *fpcovar, *fpjack;
  int *jackvect,*rjackvect, *isurvey;
  double ngalx,xx;

  double **nd, ***ndj;

  //variables for dn4k split
  float dn4k_split, *gal_dn4k, *gal_z;
  int ABOVE_SPLIT, SPLIT=0;
  FILE *fpdn4k;


  // variables for the meshlink
  float *rsqr, box_min, box_max, rmax_small, rmax_small_sqr, xmin, xmax, rsearch;
  int *meshparts, ***meshstart,nmesh,meshfac,nbrmax,*indx;
  int *meshpartsr, ***meshstartr,nmeshr,meshfacr,j1,i1;
  int *meshpartsr2, ***meshstartr2,nmeshr2,meshfacr2, ismall, nrand1_large, nrandx;

  // omp stuff
  int irank, nrank, flag=0;
  float t0,t1,ctime1,ctime2,*rx,**pix;

  int p, ix,iy,iz,iix,iiy,iiz,ir,nbr,iiix,iiiy,iiiz;
  float rmax2, side, side2, sinv;

  // using actual redshift, not cz
  SPEED_OF_LIGHT = 1;

  if(argc<4)
    {
      fprintf(stderr,"wp_LS_weight galdat1 rand1.dat RR_file [covarfile] [collision_weight1] [collision_weight2] [rmin] [rmax] [nrbin] [njack_per_side] [pi_max]> wp.dat\n");
      exit(0);
    }

  COLLISION_WEIGHT = 2.64;
  if(argc>5)
    COLLISION_WEIGHT = atof(argv[5]);
  if(argc>6)
    COLLISION_WEIGHT2 = atof(argv[6]);
  fprintf(stderr,"collision_weight= %f\n",COLLISION_WEIGHT);
  fprintf(stderr,"collision_weight2= %f\n",COLLISION_WEIGHT2);

  rmin = 0.2;
  rmax = 60.0;
  nr = 10;
  dlogr = log(rmax/rmin)/(nr);

  zmin = 0;
  PI_MAX = zmax = 60;
  if(argc>11)
    PI_MAX=zmax=atof(argv[11]);
  nz = (int)zmax;
  deltaz = (zmax-zmin)/nz;

  // let's see if we want to input the binning
  if(argc>7)
    rmin = atof(argv[7]);
  if(argc>8)
    rmax = atof(argv[8]);
  if(argc>9)
    nr = atoi(argv[9]);
  dlogr = log(rmax/rmin)/(nr);
  fprintf(stderr,"rmin= %f, rmax= %f, nrbin= %d\n",rmin, rmax, nr);

  njack1 = 5;
  if(argc>10)
    njack1 = atoi(argv[10]);
  njack_tot=njack1*njack1;
  fprintf(stderr,"Njack1= %d (total jacks= %d)\n",njack1, njack_tot);


  rsearch = sqrt(rmax*rmax + PI_MAX*PI_MAX);


  //stats arrays
  npairs = dmatrix(1,nr,1,nz);
  npairs_dr1 = dmatrix(1,nr,1,nz);
  npairs_dr2 = imatrix(1,nr,1,nz);
  rbar = dvector(1,nr);
  pibar = dmatrix(1,nr,1,nz);

  for(i=1;i<=nr;++i)
    for(j=1;j<=nz;++j)
	  npairs[i][j] = rbar[i] = npairs_dr1[i][j] = npairs_dr2[i][j] = pibar[i][j] = 0;

  //jacks arrays
  npairs_jack = d3tensor(0,njack_tot-1,1,nr,1,nz);
  npairs_dr1_jack = d3tensor(0,njack_tot-1,1,nr,1,nz);

  for(i=0;i<njack_tot;++i)
    for(j=1;j<=nr;++j)
      for(k=1;k<=nz;++k)
	npairs_jack[i][j][k] = npairs_dr1_jack[i][j][k] = 0;

  // read in all galaxies
  fp = openfile(argv[1]);
  ngal_tot = filesize(fp);

  // mock input
  fp = openfile(argv[1]);
  ngal = filesize(fp);
  
  x = vector(1,ngal);
  y = vector(1,ngal);
  z = vector(1,ngal);
  jackvect = ivector(1,ngal);
  isurvey = ivector(1,ngal);
  gweight = vector(1,ngal);
  gal_ra = vector(1,ngal);
  gal_dec = vector(1,ngal);
  gal_dec2 = vector(1,ngal);
  gal_z = vector(1,ngal);
  
  // jackknife array
  jackvect = ivector(1,ngal);

  // temp array
  xg = vector(1,ngal);
  yg = vector(1,ngal);
  subindx = ivector(1,ngal);
  indx = ivector(1,ngal);

  for(j=1;j<=ngal;++j)
    {
      fscanf(fp,"%f %f %f %f %d",&x1,&y1,&z1,&gweight[j],&isurvey[j]);
      xg[j] = x1;
      yg[j] = y1;

      jackvect[j] = 1;
      z1 /= SPEED_OF_LIGHT;
      fgets(aa,1000,fp);

      gal_dec2[j] = y1*PI/180.;      
      phi = x1;
      theta = PI/2.0 - y1;

      // assuming RA/DEC input in degrees
      phi = x1*PI/180.0;
      theta = PI/2.0 - y1*PI/180.0;
      
      zz1 = z1;
      z1 = redshift_distance(z1);
      
      x[j] = z1*sin(theta)*cos(phi);
      y[j] = z1*sin(theta)*sin(phi);
      z[j] = z1*cos(theta);
      
      gal_ra[j] = phi;
      gal_dec[j] = cos(theta);
      gal_z[j] = zz1;
      indx[j] = j;
    }
  fprintf(stderr,"Read [%d/%d] galaxies from file [%s]\n",j,ngal,argv[1]);
  fflush(stdout);
  

  // sort everything in dec (xg is not actually used)
  sort2(ngal, yg, indx);

  // now go in quintiles of dec and sort by ra.
  nby5 = ngal/njack1;
  for(i=0;i<njack1;++i)
    {
      for(n=0,j=i*nby5+1;j<=(i+1)*nby5;++j)
	{
	  n++;
	  id = indx[j];
	  yg[n] = gal_ra[id];
	  subindx[n] = id;
	}
      // sort by ra
      sort2(n,yg,subindx);
      // now divide these into groups of five
      nsby5 = n/njack1;
      for(ii=0;ii<njack1;++ii)
	for(k=0,j=ii*nsby5+1;j<=(ii+1)*nsby5;++j)
	  {
	    k++;
	    id = subindx[j];
	    jackvect[id] = i + ii*njack1;
	    //printf("GAL %f %f %d %d %d %d\n",gal_ra[id],gal_dec2[id],jackvect[id],i,ii,id);
	  }      
    }
  // temp array
  free_vector(xg,1,ngal);
  free_vector(yg,1,ngal);
  free_ivector(subindx,1,ngal);
  free_ivector(indx,1,ngal);
  for(i=1;i<=-ngal;++i)
    printf("GAL %f %f %d\n",gal_ra[i],gal_dec2[i],jackvect[i]);


  // NBNBNB replace temporarily
  for(i=1;i<=ngal;++i)
    {
      //jackvect[i] = isurvey[i];
      isurvey[i] = 1;
    }

  // read in first list of randoms
  fp = openfile(argv[2]);
  nrand1 = filesize(fp);

  rx1 = vector(1,nrand1);
  ry1 = vector(1,nrand1);
  rz1 = vector(1,nrand1);
  rjackvect = ivector(1,nrand1);

  ran_ra = vector(1,nrand1);
  ran_dec = vector(1,nrand1);

  //temp
  xg = vector(1,nrand1);
  yg = vector(1,nrand1);
  subindx = ivector(1,nrand1);
  indx = ivector(1,nrand1);

  for(i=1;i<=nrand1;++i)
    {
      fscanf(fp,"%f %f",&rx1[i],&ry1[i]);
      //fscanf(fp,"%d",&rjackvect[i]);
      yg[i] = ry1[i];
      indx[i] = i;
      //rjackvect[i] = 1; //NB: don't currently have a jackknife samples made up yet.
      ran_ra[i] = rx1[i]*PI/180.;
      ran_dec[i] = ry1[i]*PI/180.;
      phi = rx1[i]*(PI/180.0);
      theta = PI/2.0 - ry1[i]*(PI/180.0);
      // give this a random redshift from the galaxy sample
      zz1 = z1 = random_redshift_from_data(ngal,gal_z);
      
      z1 = redshift_distance(z1);
      fgets(aa,1000,fp);

      rx1[i] = z1*sin(theta)*cos(phi);
      ry1[i] = z1*sin(theta)*sin(phi);
      rz1[i] = z1*cos(theta);
    }
  fclose(fp);
  fprintf(stderr,"Read [%d] randoms from file [%s]\n",nrand1,argv[2]);

  // NBNBNB-- skip this
  //goto SKIP_JACKS;


  // sort everything in dec (xg is not actually used)
  sort2(nrand1, yg, indx);

  // now go in quintiles of dec and sort by ra.
  nby5 = nrand1/njack1;
  for(i=0;i<njack1;++i)
    {
      for(n=0,j=i*nby5+1;j<=(i+1)*nby5;++j)
	{
	  n++;
	  id = indx[j];
	  yg[n] = ran_ra[id];
	  subindx[n] = id;
	}
      // sort by ra
      sort2(n,yg,subindx);
      // now divide these into groups of five
      nsby5 = n/njack1;
      for(ii=0;ii<njack1;++ii)
	for(k=0,j=ii*nsby5+1;j<=(ii+1)*nsby5;++j)
	  {
	    k++;
	    id = subindx[j];
	    rjackvect[id] = i + ii*njack1;
	  }
    }
  for(i=1;i<=-nrand1;++i)
    printf("RAND %f %f %d\n",ran_ra[i],ran_dec[i],rjackvect[i]);

  // temp array
  free_vector(xg,1,nrand1);
  free_vector(yg,1,nrand1);
  free_ivector(subindx,1,nrand1);
  free_ivector(indx,1,nrand1);

 SKIP_JACKS:

  // find the range of the data
  xmax = 0;
  xmin = 100000;
  for(i=1;i<=ngal;++i)
    {
      if(x[i]<xmin)xmin=x[i];
      if(y[i]<xmin)xmin=y[i];
      if(z[i]<xmin)xmin=z[i];
      if(x[i]>xmax)xmax=x[i];
      if(y[i]>xmax)xmax=y[i];
      if(z[i]>xmax)xmax=z[i];
    }
  for(i=1;i<=nrand1;++i)
    {
      if(rx1[i]<xmin)xmin=rx1[i];
      if(ry1[i]<xmin)xmin=ry1[i];
      if(rz1[i]<xmin)xmin=rz1[i];
      if(rx1[i]>xmax)xmax=rx1[i];
      if(ry1[i]>xmax)xmax=ry1[i];
      if(rz1[i]>xmax)xmax=rz1[i];
    }
  box_min = xmin-rsearch;
  box_max = xmax+rsearch;
  fprintf(stderr,"Dimensions: %f %f\n",box_min, box_max);


  // make the meshlink for the data
  nmesh=0;
  meshlink2(ngal,&nmesh,box_min,box_max,rsearch,&x[1],&y[1],&z[1],&meshparts,&meshstart,meshfac);

  // make the meshlink for the randoms--> outer bins
  nmeshr=0;
  meshlink2(nrand1,&nmeshr,box_min,box_max,rsearch,&rx1[1],&ry1[1],&rz1[1],
	    &meshpartsr,&meshstartr,meshfacr);


  
  // if there are no ran-ran pairs, make them
  CALC_RR_COUNTS=1;
  if(CALC_RR_COUNTS) {
    fprintf(stderr,"Calculating RR pairs...\n");

    nrand = imatrix(1,nr,1,nz);
    nrand_jack = i3tensor(0,njack_tot-1,1,nr,1,nz);
    
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	nrand[i][j] = 0;

    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  nrand_jack[k][i][j] = 0;

    ctime1 = omp_get_wtime();

    /*
#pragma omp parallel
    {

#pragma omp for schedule(static)
    for(i=1;i<=16;++i)
      {
	printf("Thread %d is doing iteration %d.\n", omp_get_thread_num( ), i);
      }
    }
    exit(0);
    */

#pragma omp parallel  shared(rx1,ry1,rz1,meshpartsr, meshstartr,nrand1,flag) \
  private(nx, nxj, i, j, k, j1, i1, irank, nrank,indx,rsqr,nbrmax,ibin,jbin,dx,dy,dz,lx,ly,lz,pi,r,ctime2,rx2,ry2,rz2,meshpartsr2,meshstartr2, p, ix,iy,iz,iix,iiy,iiz,ir,rmax2,nbr,side,side2,sinv,iiix,iiiy,iiiz,id,id2)
{

    irank=omp_get_thread_num();
    nrank=omp_get_num_threads();
    if(!irank)fprintf(stderr,"Num threads: %d %d\n",nrank,nmeshr);

    // create local pair count arrays
    nx = imatrix(1,nr,1,nz);
    nxj = i3tensor(0,njack_tot-1,1,nr,1,nz);
    
    for(i1=1;i1<=nr;++i1)
      for(j=1;j<=nz;++j)
	nx[i1][j] = 0;

    for(i1=1;i1<=nr;++i1)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  nxj[k][i1][j] = 0;

    // these should be good for the randoms as well
    // do I need to declare these in the parallel region?
    indx=malloc(nrand1*sizeof(int));
    rsqr=malloc(nrand1*sizeof(float));


#pragma omp barrier

    for(i=1+irank;i<=nrand1;i+=nrank) {

      if(i%10000==1 && !irank){ fprintf(stderr,"%d\n",i); }
      nbrmax=nrand1;
      nbrsfind2(box_min,box_max,rsearch,nmeshr,rx1[i],ry1[i],rz1[i],&nbrmax,indx,rsqr,
		&rx1[1],&ry1[1],&rz1[1],meshpartsr,meshstartr,-1);

      for(j1=0;j1<nbrmax;++j1)
	{
	  j = indx[j1]+1;
	  if(j<=i)continue;
	  dx = (rx1[i] - rx1[j]);
	  //if(fabs(dx)>rmax*2)continue; // CHECK: why comment these out?
	  dy = (ry1[i] - ry1[j]);
	  //if(fabs(dy)>rmax*2)continue;
	  dz = (rz1[i] - rz1[j]);
	  //if(fabs(dz)>rmax*2)continue;
	  
	  lx = 0.5*(rx1[i] + rx1[j]);
	  ly = 0.5*(ry1[i] + ry1[j]);
	  lz = 0.5*(rz1[i] + rz1[j]);

	  pi = mabs((dx*lx + dy*ly + dz*lz)/sqrt(lx*lx + ly*ly + lz*lz));
	  if(pi>=PI_MAX)continue;

	  r = sqrt(dx*dx + dy*dy + dz*dz - pi*pi);
	  if(r<rmin)continue;

	  ibin = (int)(log(r/rmin)/dlogr)+1;
	  jbin = (int)(pi/deltaz) + 1;
	  
	  if(ibin<=0)continue;
	  if(ibin>nr)continue;
	  if(jbin<=0)continue;
	  if(jbin>nz)continue;

	  nx[ibin][jbin]+=1;

	  // only count this pair if both points within the same subsample
	  id = rjackvect[i];
	  id2 = rjackvect[j];
	  //if(id==id2) nxj[id2][ibin][jbin]+=1;
	  //continue;
	  nxj[id][ibin][jbin]+=1;
	  if(id!=id2) nxj[id2][ibin][jbin]+=1;
	}      
    }

    free(indx);
    free(rsqr);

    //system("date");
    if(!flag) {
      flag = 1;
      ctime2 = omp_get_wtime();
      //printf("TT %.2f\n",ctime2-ctime1);
    }

    // now combine all the counts
    #pragma omp critical
    {
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	{
	  nrand[i][j] += nx[i][j];
	}
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  nrand_jack[k][i][j] += nxj[k][i][j];
    } 

  } // end pragma section


    ctime2 = omp_get_wtime();
    fprintf(stderr,"TIME %.2f\n",ctime2-ctime1);

    // output the data to a file
    fp = fopen(argv[3],"w");
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	fprintf(fp,"%d\n",nrand[i][j]);
    for(ijack=0;ijack<njack_tot;++ijack)
      for(i=1;i<=nr;++i)
	for(j=1;j<=nz;++j)
	  fprintf(fp,"%d\n",nrand_jack[ijack][i][j]);
    fclose(fp);
    fprintf(stderr,"Done with the random-random pairs.\n");
    fprintf(stderr,"Outputted random-random pairs to file [%s].\n",argv[3]);
  }

  // now find out how many galaxies/randoms are in each subsample
  galcnt_jack = vector(0,njack_tot-1);
  for(i=0;i<njack_tot;++i)
    galcnt_jack[i] = 0;

  //this is the WEIGHTED number of galaxies
  ngalx = 0;
  for(i=1;i<=ngal;++i)
    {
      ngalx += gweight[i];
      id = jackvect[i];
      if(id<0)continue;
      galcnt_jack[id]+=gweight[i];
    }
  j=0;
  for(i=0;i<njack_tot;++i)
    {
      //printf("%d %d\n",i,galcnt_jack[i]);
      j+=galcnt_jack[i];
    }
  fprintf(stderr,"done with galcnt %d %.0f %d\n",ngal,ngalx,j);

  rancnt_jack = ivector(0,njack_tot-1);
  for(i=0;i<njack_tot;++i)
    rancnt_jack[i] = 0;

  for(i=1;i<=nrand1;++i)
    {
      id = rjackvect[i];
      if(id<0)continue;
      rancnt_jack[id]++;
    }
  j=0;
  for(i=0;i<njack_tot;++i)
    j+=rancnt_jack[i];
  fprintf(stderr,"done with rancnt %d %d\n",nrand1,j);


  

#pragma omp parallel  shared(x,y,z,meshparts, meshstart,ngal,flag, gal_ra, gal_dec2) \
  private(nd, ndj,rx,pix, i, j, k, j1, i1, irank, nrank,indx,rsqr,nbrmax,ibin,jbin,dx,dy,dz,lx,ly,lz,pi,r, \
	  ctime2,rx2,ry2,rz2, p, ix,iy,iz,iix,iiy,iiz,ir,rmax2,nbr,side,side2, \
	  sinv,iiix,iiiy,iiiz,id,id2,theta, weight)
  {

    irank=omp_get_thread_num();
    nrank=omp_get_num_threads();


    // create local pair count arrays
    nd = dmatrix(1,nr,1,nz);
    ndj = d3tensor(0,njack_tot-1,1,nr,1,nz);
    rx = vector(1,nr);
    pix = matrix(1,nr,1,nz);

    for(i=1;i<=nr;++i)
      {
	rx[i] = 0;
	for(j=1;j<=nz;++j)
	  pix[i][j] = 0;
      }

    for(i1=1;i1<=nr;++i1)
      for(j=1;j<=nz;++j)
	nd[i1][j] = 0;

    for(i1=1;i1<=nr;++i1)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  ndj[k][i1][j] = 0;

    // these should be good for the randoms as well
    // do I need to declare these in the parallel region?
    indx=malloc(ngal*sizeof(int));
    rsqr=malloc(ngal*sizeof(float));

  // gal-gal pairs
    for(i=irank+1;i<=ngal;i+=nrank) {
    if(i%10000==0)fprintf(stderr,"%d\n",i);
    nbrmax=ngal;
    nbrsfind2(box_min,box_max,rsearch,nmesh,x[i],y[i],z[i],&nbrmax,indx,rsqr,&x[1],&y[1],&z[1],
	      meshparts,meshstart,i-1);
    for(j1=0;j1<nbrmax;++j1)
      {
	j = indx[j1]+1;
	dx = (x[i] - x[j]);
	//if(fabs(dx)>rmax*2)continue;
	dy = (y[i] - y[j]);
	//if(fabs(dy)>rmax*2)continue;
	dz = (z[i] - z[j]);
	//if(fabs(dz)>rmax*2)continue;
	
	lx = 0.5*(x[i] + x[j]);
	ly = 0.5*(y[i] + y[j]);
	lz = 0.5*(z[i] + z[j]);
	
	pi = fabs((dx*lx + dy*ly + dz*lz)/sqrt(lx*lx + ly*ly + lz*lz));
	if(pi>=PI_MAX)continue;

	r = sqrt(dx*dx + dy*dy + dz*dz - pi*pi);
	if(r<rmin)continue;

	ibin = (int)(log(r/rmin)/dlogr)+1;
	jbin = (int)(pi/deltaz) + 1;
	//jbin = (int)(log(pi/zmin)/dlogz) + 1;


	if(jbin<=0)continue;
	if(jbin>nz)continue;
	if(ibin<=0)continue;
	if(ibin>nr)continue;

	if(isurvey[i]== isurvey[j]) 
	  {
	    theta = angular_separation(gal_ra[i],gal_dec2[i],gal_ra[j],gal_dec2[j]);
	    weight = collision_weight(theta, isurvey[i])*gweight[i]*gweight[j];
	  }
	else
	  {
	    weight = gweight[i]*gweight[j];
	  }

	nd[ibin][jbin]+=weight;
	rx[ibin]+=r*weight;
	pix[ibin][jbin]+=pi*weight;
	      
	id = jackvect[i];
	id2 = jackvect[j];
	//if(id==id2)ndj[id2][ibin][jbin]+=weight;
	//continue;
	ndj[id][ibin][jbin]+=weight;
	if(id!=id2)
	  ndj[id2][ibin][jbin]+=weight;
      }
  }
  free(indx);
  free(rsqr);

  // now combine all the counts
#pragma omp critical
  {
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	npairs[i][j] += nd[i][j];
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	pibar[i][j] += pix[i][j];
    for(i=1;i<=nr;++i)
      rbar[i] += rx[i];

    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  npairs_jack[k][i][j] += ndj[k][i][j];
  } 

  free_dmatrix(nd,1,nr,1,nz);
  free_d3tensor(ndj, 0,njack_tot-1,1,nr,1,nz);
#pragma omp barrier
  }

  fprintf(stderr,"done with gal-gal pairs\n");

#pragma omp parallel  \
  private(nd, ndj,rx,pix, i, j, k, j1, i1, irank, nrank,indx,rsqr,nbrmax,ibin,jbin,dx,dy,dz,lx,ly,lz,pi,r, \
	  ctime2,rx2,ry2,rz2, p, ix,iy,iz,iix,iiy,iiz,ir,rmax2,nbr,side,side2, \
	  sinv,iiix,iiiy,iiiz,id,id2,theta, weight)
  {

  indx=malloc(nrand1*sizeof(int));
  rsqr=malloc(nrand1*sizeof(float));


    irank=omp_get_thread_num();
    nrank=omp_get_num_threads();

    // create local pair count arrays
    nd = dmatrix(1,nr,1,nz);
    ndj = d3tensor(0,njack_tot-1,1,nr,1,nz);

    for(i1=1;i1<=nr;++i1)
      for(j=1;j<=nz;++j)
	nd[i1][j] = 0;

    for(i1=1;i1<=nr;++i1)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  ndj[k][i1][j] = 0;

  //ga1-rand1 pairs
  for(i=1+irank;i<=ngal;i+=nrank) {
    if(i%10000==1)fprintf(stderr,"%d\n",i);
    nbrmax=nrand1;
    nbrsfind2(box_min,box_max,rsearch,nmeshr,x[i],y[i],z[i],&nbrmax,indx,rsqr,&rx1[1],&ry1[1],&rz1[1],
		meshpartsr,meshstartr,-1);
    for(j1=0;j1<nbrmax;++j1)
      {
	j = indx[j1]+1;
	dx = (x[i] - rx1[j]);
	//if(fabs(dx)>rmax*2)continue;
	dy = (y[i] - ry1[j]);
	//if(fabs(dy)>rmax*2)continue;
	dz = (z[i] - rz1[j]);
	//if(fabs(dz)>rmax*2)continue;
	
	lx = 0.5*(x[i] + rx1[j]);
	ly = 0.5*(y[i] + ry1[j]);
	lz = 0.5*(z[i] + rz1[j]);
	
	pi = fabs((dx*lx + dy*ly + dz*lz)/sqrt(lx*lx + ly*ly + lz*lz));
	if(pi>=PI_MAX)continue;

	r = sqrt(dx*dx + dy*dy + dz*dz - pi*pi);
	if(r<rmin)continue;

	ibin = (int)(log(r/rmin)/dlogr)+1;
	jbin = (int)(pi/deltaz) + 1;
	//jbin = (int)(log(pi/zmin)/dlogz) + 1;

	if(jbin<=0)continue;
	if(jbin>nz)continue;

	if(ibin<=0)continue;
	if(ibin>nr)continue;

	weight = gweight[i];
	nd[ibin][jbin]+=weight;

	id = jackvect[i];
	id2 = rjackvect[j];
	//if(id==id2)ndj[id2][ibin][jbin]+=weight;
	//continue;
	ndj[id][ibin][jbin]+=weight;
	if(id!=id2)ndj[id2][ibin][jbin]+=weight;
      }	
  }

  // now combine all the counts
#pragma omp critical
  {
    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	npairs_dr1[i][j] += nd[i][j];

    for(i=1;i<=nr;++i)
      for(j=1;j<=nz;++j)
	for(k=0;k<njack_tot;++k)
	  npairs_dr1_jack[k][i][j] += ndj[k][i][j];
  }

  }
  fprintf(stderr,"done with gal-rand pairs\n");




  fp = fopen("xi2d.dat","w");
  
  if(argc>4)
    fpcovar = fopen(argv[4],"w");
  else
    fpcovar = fopen("covar.dat","w");
  
  j = 1;
  // cut out the first r-bin?

  for(i=1;i<=nr;++i)
    {
      npairs_tot = 0;
      for(j=1;j<=nz;++j)
	npairs_tot+= npairs[i][j];
      r = rbar[i]/npairs_tot;

      if(!npairs_tot)r = exp(dlogr*(i-0.5))*rmin;

      xi_data[i] = 0;
      for(j=1;j<=nz;++j)
	{
	  DD = npairs[i][j]/(1.*ngalx*ngalx/2);
	  RR = nrand[i][j]/(1.*nrand1*nrand1/2);
	  DR1 = npairs_dr1[i][j]/(1.*nrand1*ngalx);
	  xi = (DD - 2*DR1 + RR)/RR;
	  if(xi<-1)xi==-1;
	  xi_data[i] += xi*(zmax/nz);
	  fprintf(fp,"%e %e %e %e %e %e %.1f %d\n",r,pibar[i][j]/npairs[i][j],
		  (DD - 2*DR1 + RR)/RR,DD,RR,DR1,npairs[i][j],nrand[i][j]);
	}
      xi_data[i] *= 2;
      xi = xi_data[i];
      //printf("%11.5f %10.1f %.4e\n",r,npairs_tot,xi);
      fflush(stdout);

      // do the jack knife errors
      xi_err = 0;
      xi_bar = 0;
      for(k=0;k<njack_tot;++k)
	{
	  ngal_tmp = ngalx - galcnt_jack[k];
	  nrand_tmp = nrand1 - rancnt_jack[k];
	  xi_jack = 0;
	  for(j=1;j<=nz;++j)
	    {
	      DD = (npairs[i][j]-npairs_jack[k][i][j])/(1.*ngal_tmp*ngal_tmp/2.);
	      RR = (nrand[i][j]-nrand_jack[k][i][j])/(1.*nrand_tmp*nrand_tmp/2.);
	      DR1 = (npairs_dr1[i][j]-npairs_dr1_jack[k][i][j])/(1.*nrand_tmp*ngal_tmp);
	      xx = (DD - 2*DR1 + RR)/RR;
	      if(isnan(xx) || isinf(xx) || xx<-1)xx=-1;
	      xi_jack += xx;
	    }
	  xi_jack *= 2;
	  if(isnan(xi_jack))xi_jack=xi;
	  xi_bar += xi_jack/njack_tot;
	  if(i==-nr)
	    printf("%d %d %e %e %e %.0f %d\n",i,k,xi_jack,xi,xi_bar,galcnt_jack[k],rancnt_jack[k]);
	}
      xi_data[i] = xi_bar;

      xi_err = 0;
      for(k=0;k<njack_tot;++k)
	{
	  ngal_tmp = ngalx - galcnt_jack[k];
	  nrand_tmp = nrand1 - rancnt_jack[k];
	  xi_jack = 0;
	  for(j=1;j<=nz;++j)
	    {
	      DD = (npairs[i][j]-npairs_jack[k][i][j])/(1.*ngal_tmp*ngal_tmp/2.);
	      RR = (nrand[i][j]-nrand_jack[k][i][j])/(1.*nrand_tmp*nrand_tmp/2.);
	      DR1 = (npairs_dr1[i][j]-npairs_dr1_jack[k][i][j])/(1.*nrand_tmp*ngal_tmp);
	      xx = (DD - 2*DR1 + RR)/RR;
	      if(isnan(xx) || isinf(xx) || xx<-1)xx=-1;
	      xi_jack += xx;
	    }
	  xi_jack *= 2;
	  if(isnan(xi_jack))xi_jack=xi;
	  xi_err += (xi_jack-xi_bar)*(xi_jack-xi_bar);
	  xij[k][i]  = xi_jack;
	}
      xi_err = sqrt((njack_tot-1.)/njack_tot*xi_err);

      printf("%11.5f %.4e %.4e %10.1f %.4e\n",r,xi,xi_err,npairs_tot,xi_bar);
    }
  fclose(fp);

    for(i=1;i<=nr;++i)
      for(j=1;j<=nr;++j)
	{
	  covar[i][j] = 0;
	  for(ijack=0;ijack<njack_tot;++ijack)
	    covar[i][j] += (xi_data[i] - xij[ijack][i])*(xi_data[j] - xij[ijack][j]);
	  covar[i][j] *= (njack_tot-1.0)/njack_tot;
	  fprintf(fpcovar,"%d %d %e\n",i,j,covar[i][j]);
	}
    fclose(fpcovar);
}



double redshift_distance(double z)
{
  if(z<=0)return 0;
  return c_on_H0*qromo(func_dr1,0.0,z,midpnt);
}
double func_dr1(double z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M),-0.5);
}

float random_radial_position(float zlo, float zhi)
{
  static int flag=1;
  static float zloprev=-1, zhiprev=-1, rlo, rhi;
  float r;

  if(flag)
    {
      flag = 0;
      srand48(2344234);
    }

  if(zloprev!=zlo || zhiprev!=zhi)
    {
      rlo = redshift_distance(zlo);
      rhi = redshift_distance(zhi);
      zloprev = zlo;
      zhiprev = zhi;
    }

  while(1)
    {
      r = drand48()*(rhi-rlo)+rlo;
      if(drand48()<pow(r/rhi,2.0))return r;
    }
  return 0;

}

float angular_separation(float a1, float d1, float a2, float d2)
{
  float cd1,cd2,sd1,sd2,ca1a2,sa1a2;

  return atan((sqrt(cos(d2)*cos(d2)*sin(a2-a1)*sin(a2-a1) + 
		    pow(cos(d1)*sin(d2) - sin(d1)*cos(d2)*cos(a2-a1),2.0)))/
	      (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a2-a1)));
}

float random_redshift_from_data(int n, float *z)
{
  static int flag=1;
  if(flag)
    {
      srand48(234099);
      flag=0;
    }
  return z[(int)(drand48()*n)+1];
}

float collision_weight(float theta, int iflag)
{
  if(theta<COLLISION_ANGLE)
    {
      if(iflag==1)
	return COLLISION_WEIGHT;
      if(iflag==2)
	return COLLISION_WEIGHT2;
    }
  return 1;
}

