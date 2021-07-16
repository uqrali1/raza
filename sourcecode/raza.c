
/* rubbiya
 * 3d_log - Three Dimentional Laplacian of Gaussian, z-crossings Filter, 3D
 *
 * Author: Rubbiya Ali
 *
 * Copyright Institute for Molecular Bioscience, The University of Queensland, Australia
 * Incorporated into IMOD with permission
 * Changes for IMOD are confined to main() amd associated routines at end
 *
 *  $Id: 3d_log_zerocrossing.c,v 1419291aacca??? 2017/08/31 10:28:11 rubbiya 
 */

/* IMOD modifications above here */
/*--------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include <stdbool.h>
#define FREE_ARG char*

#include "imodel.h"
#include "b3dutil.h"

#include <sys/stat.h>
//s#include <string.h> // for strcat() call

#define PI 3.14159

typedef struct {
    float ObjectLabel;
    long index;

    long minx;
    long maxx;
    long minz;
    long miny;
    long maxy;
    long maxz;

    long majorx1;
    long majorx2;

    long majory1;
    long majory2;

    long majorz1;
    long majorz2;


    long minorx1;
    long minorx2;

    long minory1;
    long minory2;

    long minorz1;
    long minorz2;

    long middlex1;
    long middlex2;

    long middley1;
    long middley2;

    long middlez1;
    long middlez2;

    long surf;
    long vol;
    int cont;

    float POI_maj_min_x;
    float POI_maj_min_y;
    float POI_maj_min_z;



} LabelInfo;

int sort(const void *x, const void *y) {
    return (*(int*)x - *(int*)y);
}

float euclidean(long a, long b, long c, long d, long e, long f) {
    return sqrt((a-d)*(a-d)+(b-e)*(b-e)+(c-f)*(c-f));
}

long * findcenter(long* arrayx1, long* arrayy1,long* arrayz1,long* arrayx2,long* arrayy2,long* arrayz2, long newcounter1,long newcounter2){
    long x1, x2, y1, y2, z1, z2, i,j;
    x1=0;
    x2=0;
    y1=0;
    y2=0;
    z1=0;
    z2=0;

    float min=0;
    for (i=0; i<newcounter1; i++) {
        for (j=0; j<newcounter2; j++) {
            float val=euclidean(arrayx1[i],arrayy1[i],arrayz1[i],arrayx2[j],arrayy2[j],arrayz2[j]);
            if (val>min){
                min=val;
                x1=arrayx1[i];
                y1=arrayy1[i];
                z1=arrayz1[i];

                x2=arrayx2[j];
                y2=arrayy2[j];
                z2=arrayz2[j];
            }

        }
    }
    long *res[] ={x1, y1, z1, x2, y2, z2};
    int slice=(z1+z2)/2;
    //printf(" \n%ld %ld %ld %ld %ld %ld %f Slice Number %d", x1,y1,z1,x2,y2,z2,min, slice);
    return res;
}


/*********
 ** createKernelLoG
 *********/
float ***createKernelLoG_NonSeparable3(long     nx,          /* image dimension in x direction */
                                       long     ny,          /* image dimension in y direction */
                                       long     nz,          /* image dimension in z direction */
                                       float    sigma,			/* Sigma value for LoG kernel */
                                       float    ***kernel) {

    float cst = -(1.0 / (M_PI * pow(sigma,4)));
    float dem = 2.0 * pow(sigma,2);
    int size = (sigma*6.0)+1;
    float x, y, z;
    int  k, l, m;
    int halfsize = (size) / 2;

    printf(" /***********************************/\n");
    printf(" /****** KERNEL INITIALIZATION ******/\n");
    printf(" /***********************************/\n");

    printf(" Kernel Width:    %d\n Kernel Size:     %d \n", size, size * size * size );

    for(k=0; k<size; k++)
        for(l=0; l<size; l++)
            for(m=0; m<size; m++) {
                x = (k-halfsize)*(k-halfsize);
                y = (l-halfsize)*(l-halfsize);
                z = (m-halfsize)*(m-halfsize);
                kernel[k][l][m] = cst * (1.0-(x+y+z)/dem) * exp(-(x+y+z)/dem);
                //                	printf("%f\n",kernel[k][l][m]);
            }

    return kernel;
}

float ***getNeighborhoodXYZ(long     x,          /* xCoordinate */
                            long     y,          /* yCoordinate */
                            long     z,          /* zCoordinate */
                            long     nx,         /* image dimension in x direction */
                            long     ny,         /* image dimension in y direction */
                            long     nz,         /* image dimension in z direction */
                            int      size,		 /* Sigma value for LoG kernel */
                            float    *data,   /* input: original image;   */
                            float    ***neigh)   /* 3D neighbors */
{

    //    printf("/********************IN NEIGHBORHOOD******************/");
    int counter=0;
    long     m, n, o;         /* loop variables */
    int     pixelIndex =   0;

    int radius= size/2;
    //    printf(" radius= %d \n", radius );

    for (m=-radius; m<=radius; m++)
        for (n=-radius; n<=radius; n++)
            for (o=-radius; o<=radius; o++)
            {
                pixelIndex = ((z+m) * nx * ny) + ((y+n) * nx) + (x+o);
                neigh[m+radius][n+radius][o+radius]=data[pixelIndex];
                //                printf("%d.    %d   %d  %d  %d       %f       %f \n", counter++, x, y, z, pixelIndex,neigh[m+radius][n+radius][o+radius], data[pixelIndex] );

            }
    counter=0;
    return neigh;
}

/**
 * Apply a Laplacian of Gaussian 3D.
 * Non Separable implementation.
 */
float threed_log_NonSeparable(long     nx,          /* image dimension in x direction */
                              long     ny,          /* image dimension in y direction */
                              long     nz,          /* image dimension in z direction */
                              float    sigmaX,			/* Sigma value for kernel */
                              float    ***u)        /* input: original image;  output: smoothed */
{

    long    i, j, k, x, y, z;         /* loop variables */

    int l, m;

    float   ***input;                              /* work copy of u */
    float   ***kernel;                             /* work copy of kernel */
    float   ***neigh;
    float   ***output; //short or float?        /* work copy of kernel */

 //   printf("log 1   \n");


    int size = (sigmaX*6.0)+1;
    int pixelIndex =   0;
   // printf("log 2   \n");

    input = f3tensor (0, nx+1,0,ny+1,0,nz+1 );
  //  printf("log 3   \n");

    //kernel = f3tensor (0, size+1, 0, size+1, 0, size+1);
    kernel = f3tensor (0, size, 0, size, 0, size);
  //  printf("log 4   \n");

    output = f3tensor (0, nx+1,0 ,ny+1,0 ,nz+1 );
  //  printf("log 5   \n");


    float   *data       = (float*)malloc((nx) * (ny) * (nz) * sizeof(float));
  //  printf("log 6   \n");


    if (data == 0)
    {
        printf("ERROR: Out of memory\n");
        return 1;
    }


    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                input[i][j][k] = u[i][j][k];

            }

    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++){
                //    printf("%f    \n", input[i][j][k] );

            }

    createKernelLoG_NonSeparable3(nx, ny, nz, sigmaX, kernel);

    //            for (k=0; k<size; k++)
    //        		for (j=0; j<size; j++)
    //        			for (i=0; i<size; i++)
    //        			{
    //                        printf(" %f \n", kernel[i][j][k] );
    //        			}
    int dimx = size;
    // neigh = f3tensor (0,dimx+1,0,dimx+1,0,dimx+1); // it can not be free_f3tensor because it is referenced from another function 'getNeighborhoodXYZ'
    neigh = f3tensor (0,size,0,size,0,size); // it can not be free_f3tensor because it is referenced from another function 'getNeighborhoodXYZ'
    printf(" *\n");

    for (k=0; k<nz; k++)
        for (j=0; j<ny; j++)
            for (i=0; i<nx; i++)
            {
                pixelIndex = k * (nx * ny) + j * nx + i;
                data[pixelIndex]= 0;
                // printf(" %d.  %f    \n", pixelIndex, data[pixelIndex]);
            }
    printf(" **\n");

    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                output[i][j][k]= 0;
                //                printf(" %d.  %f    \n", pixelIndex, data[pixelIndex]);
            }
    /* Comment: in 1D array the loops should start from zero */
    printf(" ***\n");

    for (k=size; k<=nz-size; k++)
        for (j=size; j<=ny-size; j++)
            for (i= size; i<=nx-size; i++)
            {
                pixelIndex = k * (nx * ny) + j * nx + i;
                data[pixelIndex]= input[i][j][k];
                //                 printf("%d.   %f   \n", pixelIndex, data[pixelIndex] );
            }
    printf(" ****\n");

    /*************** working perfectly till here ****************/

    printf("\n Step I:   Kernel initialization done...\n");

    float pix;


    printf(" *\n");


    for (z=size; z<=nz-size; z++)
        for (y=size; y<=ny-size; y++)
            for (x=size; x<=nx-size; x++)
            {
                //                if ((x==7)&&(y==7)&&(z==7))
                //                    {
                pix = 0.0;
                getNeighborhoodXYZ(x, y, z, nx, ny, nz, size , data, neigh );

                for (k=0; k<dimx; k++)
                    for (l=0; l<dimx; l++)
                        for (m=0; m<dimx; m++)
                            //        printf("%d.     %d      %d      %d      %f      %f\n", counter++,x,y,z,neigh[k][l][m], kernel[k][l][m]);
                            pix += neigh[k][l][m]*kernel[k][l][m];
                output[x][y][z]= pix;
                //                    output[x][y][z]= 255; //only    [7][7][7] is 255

            }
    printf(" **\n");

    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                u[i][j][k] =output[i][j][k];
                //      printf(" %f    \n", u[i][j][k]);
            }
    printf(" ***\n");

    free_f3tensor (input, 0,nx+1,0,ny+1,0,nz+1);
    free_f3tensor (kernel, 0, size , 0, size, 0, size);
    free_f3tensor (neigh, 0, size , 0, size, 0, size);
    free_f3tensor (output, 0,nx+1,0,ny+1,0,nz+1);
    free(data);
    printf(" Step II:  LoG Done ...\n");

    return;
}

/* ---------------------------------------------------------------------- */

float threed_zcrossings

(float    ht,          /* time step size, 0 < ht <= 0.25 */
 long     nx,          /* image dimension in x direction */
 long     ny,          /* image dimension in y direction */
 long     nz,          /* image dimension in z direction */
 float    hx,          /* pixel size in x direction */
 float    hy,          /* pixel size in y direction */
 float    hz,          /* pixel size in z direction */
 // float    sigma,			/* Sigma value for LoG kernel */
 float    zcvalue,			/* Z crossing value */

 float    ***u)        /* input: original image;  output: binary */

/* apply 3 dimensional z crossings */
{
    /* arbitrary zero crossing on laplacian of gaussian */

    printf(" \t\tApplying arbitrary 3D Zcrossing ...\n");
    float   ***f;                                          /* work copy of u */
    float   ***logVolume;
    float   ***transformedlogVolume;

    float   ***zclogVolume;
    long    i, j, k, x, y, z;         /* loop variables */
    float   zcThreshold             =   zcvalue;
    bool ts=false;

    float   max, min;               /* largest, smallest grey value */
    float   mean;                   /* average grey value */
    float   vari;
    int size =0;


    /* ---- allocate storage ---- */

    f		= f3tensor (0,nx+1,0,ny+1,0,nz+1);
    logVolume   = f3tensor (0,nx+1,0,ny+1,0,nz+1);

    transformedlogVolume   = f3tensor (0,nx+1,0,ny+1,0,nz+1);
    zclogVolume = f3tensor (0,nx+1,0,ny+1,0,nz+1);

    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                zclogVolume[i][j][k]= 0;
                // printf(" %d. %f \n", pixelIndex, data[pixelIndex]);
            }
    const int NUM_NEIGHBOURS = 26;

    //float   *neighbour  = (float*)malloc((26 + 1) * sizeof(float)); // valgrind error ' Invalid write of size 4 '
    float   *neighbour  = malloc((NUM_NEIGHBOURS+1) * sizeof(*neighbour)); // valgrind solution

    printf(" before break ..... \n");

    for (size = 0; size <= NUM_NEIGHBOURS; size++)       // valgrind
        neighbour[size] = 0;

    printf(" after break ..... \n");

    analyse (u, nx, ny, nz, &min, &max, &mean, &vari);
    printf(" \t\tMean:          %1.10f \n", mean);

    /* ---- copy u into f ---- */

    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                logVolume[i][j][k] = u[i][j][k];

            }

    for (z=1; z<nz; z++) {
        for (y=1; y<ny; y++) {
            for (x=1; x<nx; x++){
                //printf("-  %f  \n", f[x][y][z]);

                neighbour[0]=  logVolume[x-1][y-1][z-1];
                neighbour[1]=  logVolume[x][y-1][z-1];
                neighbour[2]=  logVolume[x+1][y-1][z-1];
                neighbour[3]=  logVolume[x-1][y][z-1];
                neighbour[4]=  logVolume[x][y][z-1];
                neighbour[5]=  logVolume[x+1][y][z-1];
                neighbour[6]=  logVolume[x-1][y+1][z-1];
                neighbour[7]=  logVolume[x][y+1][z-1];
                neighbour[8]=  logVolume[x+1][y+1][z-1];

                neighbour[9] =  logVolume[x-1][y-1][z];
                neighbour[10]=  logVolume[x][y-1][z];
                neighbour[11]=  logVolume[x+1][y-1][z];
                neighbour[12]=  logVolume[x-1][y][z];
                neighbour[13]=  logVolume[x][y][z];
                neighbour[14]=  logVolume[x+1][y][z];
                neighbour[15]=  logVolume[x-1][y+1][z];
                neighbour[16]=  logVolume[x][y+1][z];
                neighbour[17]=  logVolume[x+1][y+1][z];

                neighbour[18]=  logVolume[x-1][y-1][z+1];
                neighbour[19]=  logVolume[x][y-1][z+1];
                neighbour[20]=  logVolume[x+1][y-1][z+1];
                neighbour[21]=  logVolume[x-1][y][z+1];
                neighbour[22]=  logVolume[x][y][z+1];
                neighbour[23]=  logVolume[x+1][y][z+1];
                neighbour[24]=  logVolume[x-1][y+1][z+1];
                neighbour[25]=  logVolume[x][y+1][z+1];
                neighbour[26]=  logVolume[x+1][y+1][z+1];
                //                printf("-  %f  \n", zcThreshold);

                //                printf("-  %f  \n", logVolume[x][y][z]);

                // zero-crossing
                if (zcThreshold==0) {
                    ts=false;
                    if (logVolume[x][y][z]<=0){
                        for (i=0;i<=26;i++){
                            if ( neighbour[i]>0)
                                ts=true;
                        }
                    }
                    if (ts)
                        zclogVolume[x][y][z]=255;
                    else
                        zclogVolume[x][y][z]=0;

                }
                // z-crossing negative
                else if (zcThreshold<0){
                    ts=false;
                    if (logVolume[x][y][z]<=zcThreshold){
                        for (i=0;i<=26;i++){
                            //	printf("-  %f  \n", neighbour[i]);
                            if ( neighbour[i]>zcThreshold)
                                ts=true;
                        }
                    }

                    if (ts)
                        zclogVolume[x][y][z]=255;
                    else
                        zclogVolume[x][y][z]=0;
                }
                //z-crossing positive
                else if (zcThreshold>0){
                    ts=false;
                    if (logVolume[x][y][z]>zcThreshold){
                        for (i=0;i<=26;i++){
                            if ( neighbour[i]<=zcThreshold)
                                ts=true;
                        }
                    }
                    if (ts)
                        zclogVolume[x][y][z]=255;
                    else
                        zclogVolume[x][y][z]=0;
                }
                //				printf("-  %f  \n", zclogVolume[x][y][z]);

            }
        }
    }

    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                if ((k==2)|(k==1)){
                    u[i-1][j-1][k]=zclogVolume[i][j][3];
                    u[i][j][k]=zclogVolume[i][j][3];

                    //                    u[i-1][j-1][k]=0;
                    //                    u[i][j][k]=0;

                }

                else if ((k==nz)|(k==(nz-1))){
                    u[i-1][j-1][k]=zclogVolume[i][j][(nz-2)];
                    u[i][j][k]=zclogVolume[i][j][(nz-2)];
                    //                    u[i-1][j-1][k]=0;
                    //                    u[i][j][k]=0;

                }
                else{
                    u[i][j][k] =zclogVolume[i][j][k];
                }
            }

    for (k=1; k<=nz; k++)
    for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++)
    {
        if((k<=5)||(k>=(nz-5))||(j<=5)||(j>=(ny-5))||(i<=5)||(i>=(nx-5)))
        u[i][j][k] =0;

    }




    free_f3tensor (zclogVolume, 0,nx+1,0,ny+1,0,nz+1);
    free_f3tensor (f, 0,nx+1,0,ny+1,0,nz+1);
    free_f3tensor (logVolume, 0,nx+1,0,ny+1,0,nz+1);
    free_f3tensor (transformedlogVolume, 0,nx+1,0,ny+1,0,nz+1);
    free(neighbour);

    printf(" Step III: Z Crossings Done ...\n");
    return;


}




void analyse

(float   ***u,         /* image, unchanged */
 long    nx,          /* pixel number in x direction */
 long    ny,          /* pixel number in x direction */
 long    nz,          /* pixel number in x direction */
 float   *min,        /* minimum, output */
 float   *max,        /* maximum, output */
 float   *mean,       /* mean, output */
 float   *vari)       /* variance, output */

/*
 calculates minimum, maximum, mean and variance of an image u
 */

{
    long    i, j, k;       /* loop variables */
    float   help;       /* auxiliary variable */
    double  help2;      /* auxiliary variable */

    *min  = u[1][1][1];
    *max  = u[1][1][1];
    help2 = 0.0;

    for (i=1; i<=nx; i++)
        for (j=1; j<=ny; j++)
            for (k=1; k<=nz; k++)
            {
                if (u[i][j][k] < *min) *min = u[i][j][k];
                if (u[i][j][k] > *max) *max = u[i][j][k];
                help2 = help2 + (double)u[i][j][k];
            }
    *mean = (float)help2 / (nx * ny * nz);

    *vari = 0.0;
    for (i=1; i<=nx; i++)
        for (j=1; j<=ny; j++)
            for (k=1; k<=nz; k++)
            {
                help  = u[i][j][k] - *mean;
                *vari = *vari + help * help;
            }
    *vari = *vari / (nx * ny * nz);

    return;

} /* analyse */


/*--------------------------------------------------------------------------*/
// IMOD modifications all below here

#include "b3dutil.h"
#include "parse_params.h"
#include "mrcfiles.h"
#include "mrcslice.h"

void usage(char *progname, float ht, float sigma, float zcvalue)
{
    printf("%s by Rubbiya Akram Ali (adapted for IMOD)\n"
           "Usage: %s [options] <input file> <output file>\nOptions:\n"
           "\t-z #\tZ Crossing value, threshold for gradients (default %.2f)\n"
           "\t-s #\tSigma of LoG filter (default %.2f) \n"
           "\t-o #\tOutput only the given Z slice (numbered from 1)\n"
           "\t-r #\tTolerance rate (default 10)\n"
           "\t-c #\tCondition for testing (1-7)\n"
           "\t\t  \t1. All following condition together\n"
           "\t\t  \t2. 3D Axes together\n"
           "\t\t  \t3. 3D volume\n"
           "\t\t  \t4. 3D surface area\n"
           "\t\t  \t5. Major axis only\n"
           "\t\t  \t6. Middle axis only\n"
           "\t\t  \t7. Minor axis only\n",
           progname, progname, sigma, zcvalue, ht);
    exit(1);
}

void testNumericEntry(char *endptr, char *argv, char *option)
{
    if (endptr == argv)
        exitError("Option %s must be followed by a number, not by %s", option, argv);
}


void threed_point_of_intersection
(
 int  x1,                   /* point1.x */
 int  y1,                   /* point1.Y */
 int  z1,                   /* point1.Z */
 int  x2,                   /* point2.x */
 int  y2,                   /* point2.Y */
 int  z2,                   /* point2.z */
 int  x3,                   /* point3.x */
 int  y3,                   /* point3.y */
 int  z3,                   /* point3.y */
 int  x4,                   /* point4.x */
 int  y4,                   /* point4.y */
 int  z4,                   /* point4.z */
 float *pointOfIntersection /* pointOfIntersection */
)

/* calculates the point of intersection between two perpendicular lines. */

{
    //    x1=1;
    //	y1=0;
    //	z1=0;
    //
    //	x3=0;
    //	y3=5;
    //	z3=5;


    /* Step 1:  Calcualate vector of each line from two given points A and B.
     V = B-A

     Follow example for explanation
     http://mathforum.org/library/drmath/view/63719.html
     */

    float vectorLineX1, vectorLineY1, vectorLineZ1, vectorLineX2, vectorLineY2, vectorLineZ2;


    float t1=0.0;
    float t2=0.0;


    vectorLineX1  = x2 - x1;
    vectorLineY1  = y2 - y1;
    vectorLineZ1  = z2 - z1;

    vectorLineX2  = x4 - x3;
    vectorLineY2  = y4 - y3;
    vectorLineZ2  = z4 - z3;

    t1 = (( (y3 - y1) * vectorLineX2) + ( (x1 - x3) * vectorLineY2) )/(vectorLineY1 * vectorLineX2 - vectorLineX1 * vectorLineY2);

    t2 = (x1 + vectorLineX1 * t1 - x3)/(vectorLineX2);

    /* Step 3: Point of intersection

     Plugging those into all six equations for x, y, and z, we get

     x = 1 + 2*2 = 5     x = 0 + 5*1 = 5
     y = 0 + 3*2 = 6     y = 5 + 1*1 = 6
     z = 0 + 1*2 = 2     z = 5 - 3*1 = 2

     So this is indeed the intersection of the lines.

     */

    //    printf("t1=%f t2=%f\n",t1,t2);


    float pointOfIntersectionX1, pointOfIntersectionY1, pointOfIntersectionZ1, pointOfIntersectionX2, pointOfIntersectionY2, pointOfIntersectionZ2;

    pointOfIntersectionX1 = x1 + vectorLineX1 * t1;
    pointOfIntersectionY1 = y1 + vectorLineY1 * t1;
    pointOfIntersectionZ1 = z1 + vectorLineZ1 * t1;

    pointOfIntersectionX2 = x3 + vectorLineX2 * t2;
    pointOfIntersectionY2 = y3 + vectorLineY2 * t2;
    pointOfIntersectionZ2 = z3 + vectorLineZ2 * t2;

    if ((pointOfIntersectionX1==pointOfIntersectionX2)&&(pointOfIntersectionY1==pointOfIntersectionY2)&&(pointOfIntersectionZ1==pointOfIntersectionZ2)) {

        pointOfIntersection[0] = pointOfIntersectionX1;
        pointOfIntersection[1] = pointOfIntersectionY1;
        pointOfIntersection[2] = pointOfIntersectionZ1;

        //		printf("POIx=%f POIy=%f POIz=%f\n", pointOfIntersection[0],  pointOfIntersection[1],  pointOfIntersection[2]);
    }
    return ;
}
/* ---------------------------------------------------------------------- */
/*
 Connected component labeling and cener finding.
 Rubbiya Ali*/
float threed_ObjectCenter
(
 float ***u,        /* image matrix */
 long  nx,          /* size in x direction */
 long  ny,          /* size in y direction */
 long  nz,          /* size in z direction */
 char *outputFilesPath,  /* Input directory path */
 float condition,
 float toleranceRate
 )

/* creates dummy boundaries by periodical continuation */

{
    printf("\t\t Finding objects...\n");

    float min;
    float max;
    float mean;
    float vari;
    float tol_min=0.0;
    float tol_max=0.0;

    float tolRate_min_major=0.0;
    float tolRate_max_major=0.0;
    float tolRate_min_middle=0.0;
    float tolRate_max_middle=0.0;
    float tolRate_min_minor=0.0;
    float tolRate_max_minor=0.0;
    float tolRate_min_surfaceArea=0.0;
    float tolRate_max_surfaceArea=0.0;
    float tolRate_min_volume=0.0;
    float tolRate_max_volume=0.0;

    int totcontour=0;
    int size=0;

    int cont=2;
    float myk=0;
    int obj_c=1;
    long max_x1=0,min_x1=0,max_y1=0,min_y1=0,max_z1=0,min_z1=0;

    analyse (u, nx, ny, nz, &min, &max, &mean, &vari);


    long    i, j, k, x, y, z, i1, index;                                       /* loop variables */
    float   ***f;                                          /* work copy of u */

    int neighbors=26;
    float  *array = (float*)malloc((neighbors+1) * sizeof(array));

    for (size = 0; size <= 27; size++)       // valgrind
        array[size] = 0;


    /* ---- allocate storage ---- */

    f		    = f3tensor (0,nx+1,0,ny+1,0,nz+1);
    float mymin;
    /* ---- copy u into f ---- */

    for (i=1; i<=nx; i++)
        for (j=1; j<=ny; j++)
            for (k=1; k<=nz; k++)
                f[i][j][k] = u[i][j][k];
    float label=-32000;



    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {

                if (i==1||i==2||i==3||i==nx-2||i==nx-1||i==nx||j==1||j==2||j==3||j==ny-2||j==ny-1||j==ny) {


                    // if (k==1||k==2||k==nz-1||k==nz){

                    f[i][j][k]=0;

                    // }
                }

            }

    /************** *********************** *****************/
    /*             CONNECTED COMPONENT LABELING             */
    /************** First Pass - Horizontal *****************/
    /*  First pass will assign the labels to all the voxels  */

    for (z=1; z<nz-1; z++) {
        for (y=1; y<ny-1; y++) {
            for (x=1; x<nx-1; x++){
                //printf("-  %f  \n", f[x][y][z]);
                long thisPixel=f[x][y][z];
                if (thisPixel!=0)

                {
                    mymin=-5;
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];



                    qsort(array, 27, sizeof(float), sort);
                    int ts=1;
                    for (i1=0; i1<=26; i1++) {
                        if ((array[i1]!=0)&(array[i1]!=255)&(array[i1]!=-5)){
                            mymin=array[i1];
                            ts=2;
                            //				printf("%d ",mymin);
                        }
                        if (ts==2)
                            break;
                    }
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];
                    if ((mymin!=0)&(mymin!=255)&(mymin!=-5)) {// we do not want background
                        //We do not want to touch background
                        if (f[x][y][z]==mymin){

                            for (i1=0; i1<=26; i1++) {
                                if ((array[i1]!=0))
                                    array[i1]=mymin;

                            }
                        }
                        else{
                            array[13]=mymin;
                            for (i1=0; i1<=26; i1++) {
                                if ((mymin!=0)&(mymin!=255)&(mymin!=-5))
                                    array[i1]=mymin;

                            }
                        }

                    }
                    else if ((f[x][y][z]==255)){
                        if ((label==-1)|(label==254))
                            label+=2;
                        label=label+0.1;
                        array[13]=label;      // incrementing value of label and assigning it to current voxel

                        if (f[x-1][y-1][z-1]==255){
                            f[x-1][y-1][z-1]=array[13];
                        }
                        if (f[x][y-1][z-1]==255){
                            f[x][y-1][z-1]=array[13];
                        }
                        if (f[x+1][y-1][z-1]==255){
                            f[x+1][y-1][z-1]=array[13];
                        }
                        if (f[x-1][y][z-1]==255){
                            f[x-1][y][z-1]=array[13];
                        }
                        if (f[x][y][z-1]==255){
                            f[x][y][z-1]=array[13];
                        }
                        if (f[x+1][y][z-1]==255){
                            f[x+1][y][z-1]=array[13];
                        }
                        if (f[x-1][y+1][z-1]==255){
                            f[x-1][y+1][z-1]=array[13];
                        }
                        if (f[x][y+1][z-1]==255){
                            f[x][y+1][z-1]=array[13];
                        }
                        if (f[x+1][y+1][z-1]==255){
                            f[x+1][y+1][z-1]=array[13];
                        }
                        if (f[x-1][y-1][z]==255){
                            f[x-1][y-1][z]=array[13];
                        }
                        if (f[x+1][y-1][z]==255){
                            f[x+1][y-1][z]=array[13];
                        }
                        if (f[x-1][y][z]==255){
                            f[x-1][y][z]=array[13];
                        }
                        if (f[x+1][y][z]==255){
                            f[x+1][y][z]=array[13];
                        }
                        if (f[x-1][y+1][z]==255){
                            f[x-1][y+1][z]=array[13];
                        }
                        if (f[x][y+1][z]==255){
                            f[x][y+1][z]=array[13];
                        }
                        if (f[x+1][y+1][z]==255){
                            f[x+1][y+1][z]=array[13];
                        }
                        if (f[x-1][y-1][z+1]==255){
                            f[x-1][y-1][z+1]=array[13];
                        }
                        if (f[x][y-1][z+1]==255){
                            f[x][y-1][z+1]=array[13];
                        }
                        if (f[x+1][y-1][z+1]==255){
                            f[x+1][y-1][z+1]=array[13];
                        }
                        if (f[x-1][y][z+1]==255){
                            f[x-1][y][z+1]=array[13];
                        }
                        if (f[x][y][z+1]==255){
                            f[x][y][z+1]=array[13];
                        }
                        if (f[x+1][y][z+1]==255){
                            f[x+1][y][z+1]=array[13];
                        }
                        if (f[x-1][y+1][z+1]==255){
                            f[x-1][y+1][z+1]=array[13];
                        }
                        if (f[x][y+1][z+1]==255){
                            f[x][y+1][z+1]=array[13];
                        }
                        if (f[x+1][y+1][z+1]==255){
                            f[x+1][y+1][z+1]=array[13];
                        }



                    }

                }

            }
        }

    }



    /************** Second Pass - Vertical *****************/
    /*  Second pass will assign the minimum label to all the connected voxel  */

    // check why the loops are working backword and their effect on output

    for (z=nz-1; z>1; z--) {
        for (x=nx-1; x>1; x--){
            for (y=ny-1; y>1; y--) {
                //printf("-  %f  \n", f[x][y][z]);
                long thisPixel=f[x][y][z];
                if (thisPixel!=0)

                {
                    mymin=-5;
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];



                    qsort(array, 27, sizeof(float), sort);
                    int ts=1;
                    for (i1=0; i1<=26; i1++) {
                        if ((array[i1]!=0)&(array[i1]!=255)&(array[i1]!=-5)){
                            mymin=array[i1];
                            ts=2;
                            //				printf("%d ",mymin);
                        }
                        if (ts==2)
                            break;
                    }
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];
                    if ((mymin!=0)&(mymin!=255)&(array[i1]!=-5)) {
                        // we do not want background
                        //We do not want to touch background
                        if (f[x][y][z]==mymin){

                            for (i1=0; i1<=26; i1++) {
                                if ((array[i1]!=0))
                                    array[i1]=mymin;

                            }
                        }
                        else{
                            array[13]=mymin;
                            for (i1=0; i1<=26; i1++) {
                                if ((mymin!=0)&(mymin!=255)&(mymin!=-5))
                                    array[i1]=mymin;

                            }
                        }

                    }

                    f[x-1][y-1][z-1]=array[0];
                    f[x][y-1][z-1]=array[1];
                    f[x+1][y-1][z-1]=array[2];
                    f[x-1][y][z-1]=array[3];
                    f[x][y][z-1]=array[4];
                    f[x+1][y][z-1]=array[5];
                    f[x-1][y+1][z-1]=array[6];
                    f[x][y+1][z-1]=array[7];
                    f[x+1][y+1][z-1]=array[8];

                    f[x-1][y-1][z]=array[9];
                    f[x][y-1][z]=array[10];
                    f[x+1][y-1][z]=array[11];
                    f[x-1][y][z]=array[12];
                    f[x][y][z]=array[13];
                    f[x+1][y][z]=array[14];
                    f[x-1][y+1][z]=array[15];
                    f[x][y+1][z]=array[16];
                    f[x+1][y+1][z]=array[17];

                    f[x-1][y-1][z+1]=array[18];
                    f[x][y-1][z+1]=array[19];
                    f[x+1][y-1][z+1]=array[20];
                    f[x-1][y][z+1]=array[21];
                    f[x][y][z+1]=array[22];
                    f[x+1][y][z+1]=array[23];
                    f[x-1][y+1][z+1]=array[24];
                    f[x][y+1][z+1]=array[25];
                    f[x+1][y+1][z+1]=array[26];

                }

            }
        }
    }
    int labelcounter=0;


    /************** Third Pass - Horizontal again *****************/
    /*  Third pass will rescan all the voxels and assign minimum labels to all connected labels */

    // check why the loops are working backword and their effect on output

    for (z=nz-1; z>1; z--) {
        for (y=ny-1; y>1; y--) {
            for (x=nx-1; x>1; x--){

                //printf("-  %f  \n", f[x][y][z]);
                long thisPixel=f[x][y][z];
                if (thisPixel!=0)
                {
                    mymin=-5;
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];



                    qsort(array, 27, sizeof(float), sort);
                    int ts=1;
                    for (i1=0; i1<=26; i1++) {
                        if ((array[i1]!=0)&(array[i1]!=255)&(array[i1]!=-5)){
                            mymin=array[i1];
                            ts=2;
                            //				printf("%d ",mymin);
                        }
                        if (ts==2)
                            break;
                    }
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];
                    if ((mymin!=0)&(mymin!=255)&(mymin!=-5)) {// we do not want background
                        //We do not want to touch background
                        if (f[x][y][z]==mymin){

                            for (i1=0; i1<=26; i1++) {
                                if ((array[i1]!=0))
                                    array[i1]=mymin;

                            }
                        }
                        else{
                            array[13]=mymin;
                            for (i1=0; i1<=26; i1++) {
                                if ((mymin!=0)&(mymin!=255)&(mymin!=-5))
                                    array[i1]=mymin;

                            }
                        }

                    }

                    f[x-1][y-1][z-1]=array[0];
                    f[x][y-1][z-1]=array[1];
                    f[x+1][y-1][z-1]=array[2];
                    f[x-1][y][z-1]=array[3];
                    f[x][y][z-1]=array[4];
                    f[x+1][y][z-1]=array[5];
                    f[x-1][y+1][z-1]=array[6];
                    f[x][y+1][z-1]=array[7];
                    f[x+1][y+1][z-1]=array[8];

                    f[x-1][y-1][z]=array[9];
                    f[x][y-1][z]=array[10];
                    f[x+1][y-1][z]=array[11];
                    f[x-1][y][z]=array[12];
                    f[x][y][z]=array[13];
                    f[x+1][y][z]=array[14];
                    f[x-1][y+1][z]=array[15];
                    f[x][y+1][z]=array[16];
                    f[x+1][y+1][z]=array[17];

                    f[x-1][y-1][z+1]=array[18];
                    f[x][y-1][z+1]=array[19];
                    f[x+1][y-1][z+1]=array[20];
                    f[x-1][y][z+1]=array[21];
                    f[x][y][z+1]=array[22];
                    f[x+1][y][z+1]=array[23];
                    f[x-1][y+1][z+1]=array[24];
                    f[x][y+1][z+1]=array[25];
                    f[x+1][y+1][z+1]=array[26];

                }

            }
        }
    }



    /************** Third Pass - Horizontal again *****************/
    /*  Third pass will rescan all the voxels and assign minimum labels to all connected labels */

    // check why the loops are working backword and their effect on output

    for (z=nz-1; z>1; z--) {
        for (x=nx-1; x>1; x--){
            for (y=ny-1; y>1; y--) {


                //printf("-  %f  \n", f[x][y][z]);
                long thisPixel=f[x][y][z];
                if (thisPixel!=0)
                {
                    mymin=-5;
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];



                    qsort(array, 27, sizeof(float), sort);
                    int ts=1;
                    for (i1=0; i1<=26; i1++) {
                        if ((array[i1]!=0)&(array[i1]!=255)&(array[i1]!=-5)){
                            mymin=array[i1];
                            ts=2;
                            //				printf("%d ",mymin);
                        }
                        if (ts==2)
                            break;
                    }
                    array[0]= f[x-1][y-1][z-1];
                    array[1]= f[x][y-1][z-1];
                    array[2]= f[x+1][y-1][z-1];
                    array[3]= f[x-1][y][z-1];
                    array[4]= f[x][y][z-1];
                    array[5]= f[x+1][y][z-1];
                    array[6]= f[x-1][y+1][z-1];
                    array[7]= f[x][y+1][z-1];
                    array[8]= f[x+1][y+1][z-1];

                    array[9]= f[x-1][y-1][z];
                    array[10]= f[x][y-1][z];
                    array[11]= f[x+1][y-1][z];
                    array[12]= f[x-1][y][z];
                    array[13]= f[x][y][z];
                    array[14]= f[x+1][y][z];
                    array[15]= f[x-1][y+1][z];
                    array[16]= f[x][y+1][z];
                    array[17]= f[x+1][y+1][z];

                    array[18]= f[x-1][y-1][z+1];
                    array[19]= f[x][y-1][z+1];
                    array[20]= f[x+1][y-1][z+1];
                    array[21]= f[x-1][y][z+1];
                    array[22]= f[x][y][z+1];
                    array[23]= f[x+1][y][z+1];
                    array[24]= f[x-1][y+1][z+1];
                    array[25]= f[x][y+1][z+1];
                    array[26]= f[x+1][y+1][z+1];
                    if ((mymin!=0)&(mymin!=255)&(mymin!=-5)) {// we do not want background
                        //We do not want to touch background
                        if (f[x][y][z]==mymin){

                            for (i1=0; i1<=26; i1++) {
                                if ((array[i1]!=0))
                                    array[i1]=mymin;

                            }
                        }
                        else{
                            array[13]=mymin;
                            for (i1=0; i1<=26; i1++) {
                                if ((mymin!=0)&(mymin!=255)&(mymin!=-5))
                                    array[i1]=mymin;
                                
                            }
                        }
                        
                    }
                    
                    f[x-1][y-1][z-1]=array[0];
                    f[x][y-1][z-1]=array[1];
                    f[x+1][y-1][z-1]=array[2];
                    f[x-1][y][z-1]=array[3];
                    f[x][y][z-1]=array[4];
                    f[x+1][y][z-1]=array[5];
                    f[x-1][y+1][z-1]=array[6];
                    f[x][y+1][z-1]=array[7];
                    f[x+1][y+1][z-1]=array[8];
                    
                    f[x-1][y-1][z]=array[9];
                    f[x][y-1][z]=array[10];
                    f[x+1][y-1][z]=array[11];
                    f[x-1][y][z]=array[12];
                    f[x][y][z]=array[13];
                    f[x+1][y][z]=array[14];
                    f[x-1][y+1][z]=array[15];
                    f[x][y+1][z]=array[16];
                    f[x+1][y+1][z]=array[17];
                    
                    f[x-1][y-1][z+1]=array[18];
                    f[x][y-1][z+1]=array[19];
                    f[x+1][y-1][z+1]=array[20];
                    f[x-1][y][z+1]=array[21];
                    f[x][y][z+1]=array[22];
                    f[x+1][y][z+1]=array[23];
                    f[x-1][y+1][z+1]=array[24];
                    f[x][y+1][z+1]=array[25];
                    f[x+1][y+1][z+1]=array[26];
                    
                }
                
            }
        }
    }
    
    for (z=nz; z>(nz-1); z--) {
        for (y=ny-1; y>1; y--) {
            for (x=nx-1; x>1; x--){
                f[x][y][z]=f[x][y][z-1];
            }
        }
    }
    
    /*  After Third pass f[i][j][k] contains all the labels or zeroos  */


    float  *labelarray = (float*)malloc(((nx * ny * nz) +1) * sizeof(labelarray));


    for (size = 0; size <= ((nx * ny * nz) +1); size++)       // valgrind
        labelarray[size] = 0;


    int tss=0;
    long jj=0;

    for (i=1; i<=nx; i++)
        for (j=1; j<=ny; j++)
            for (k=1; k<=nz; k++){
                if ((f[i][j][k]!=0)&(f[i][j][k]!=255)) {
                    tss=1;

                    for (index=0; index<labelcounter; index++) {
                        if (f[i][j][k]==labelarray[index]) {
                            tss=2;
                        }
                        if (tss==2)
                            break;
                    }
                    if (tss==1){
                        labelarray[labelcounter++] = f[i][j][k];  //doing labelcounter ++ to count the total number of labels assigned
                        /*  Here all the labels from f[i][j][k] will be copied in labelarray  */
                        //						printf("objects %f %f %d\n", f[i][j][k], labelarray[labelcounter-1],k);
                    }
                }
            }

    //  LabelInfo *Object =malloc(labelcounter*sizeof(LabelInfo)+1000) ;//=(long*)malloc(sizeof(long)); // declared array(Equal to total labels) to store each label detail in

    char line[1000]="";
    int lineLength=0;
    long Surface[1000];
    long Volume[1000];

    int objectCount=0;
    LabelInfo Object_datau[1000] ;//=(long*)malloc(sizeof(long)); // declared array(Equal to total labels) to store each label detail in

    char *outputDirectory= (char*)malloc((strlen(outputFilesPath)+20) * sizeof(outputDirectory));

    char *outputDirectory_obj= (char*)malloc((strlen(outputFilesPath)+20) * sizeof(outputDirectory_obj));


    for (size = 0; size < (strlen(outputFilesPath)+20); size++)
        outputDirectory[size] = NULL;

    for (size = 0; size < (strlen(outputFilesPath)+20); size++)
        outputDirectory_obj[size] = NULL;

    FILE *writeOutFile1=NULL,*writeOutFile2=NULL,*writeOutFile3=NULL, *infile=NULL;

    char fileNameWithPath1[300];
    FILE *objectFile1=NULL;
    char fileName1[500]="";

    char fileNameWithPath[300]="";

    LabelInfo mmydata1[5500];// =malloc(labelcounter*sizeof(LabelInfo)+100) ;//=(long*)malloc(sizeof(long));  //
    int raza=0;

    /* ---- allocate storage ---- */

    int ***myObject;
    myObject		    = f3tensor (0,nx*ny*nz+1,0,3,0,2);


    strcpy(outputDirectory,outputFilesPath);
    strcat(outputDirectory, "");

    strcpy(fileNameWithPath1, outputDirectory);
    strcpy(fileName1, "UserInputModel_Model2Point.txt");//fileName1="UserInputModel_Model2Point.txt";

    strcat(fileNameWithPath1, fileName1);
    infile      = fopen (&fileNameWithPath1, "rt");  /* open the file for reading text*/
    int val=0;;

    printf("\t\tInput 1: %s - ", fileName1);

    if(infile)
    {

        printf("Started reading... \n");

        while((fgets(line, sizeof(line), infile)) != NULL){

            //     printf("Entered while loop below above \n");

            lineLength= strlen(line);

            char tempobject [100]="";
            char tempobject1 [100]="";

            char tempobjecta [100]="";
            char tempobject1a [100]="";

            char tempobjectx [100]="";
            char tempobject1x [100]="";

            char tempobjecty [100]="";
            char tempobject1y [100]="";

            char tempobjectz [100]="";
            char tempobject1z [100]="";

            raza=0;
            while ((line[raza]==' ')){
                raza++;
            }

            while ((line[raza]!=' ')){
                sprintf(tempobject, "%c", line[raza]);
                strcat(tempobject1,tempobject);
                raza++;
            }

            while ((line[raza]==' ')){

                raza++;
            }

            while ((line[raza]!=' ')){
                sprintf(tempobjecta, "%c", line[raza]);
                strcat(tempobject1a,tempobjecta);
                raza++;
            }

            while ((line[raza]==' ')){

                raza++;
            }
            while ((line[raza]!=' ')){
                sprintf(tempobjectx, "%c", line[raza]);
                strcat(tempobject1x,tempobjectx);
                raza++;
            }

            while ((line[raza]==' ')){

                raza++;
            }

            while ((line[raza]!=' ')){
                sprintf(tempobjecty, "%c", line[raza]);
                strcat(tempobject1y,tempobjecty);
                raza++;
            }

            while ((line[raza]==' ')){

                raza++;
            }

            while (raza!=(lineLength-1)){
                sprintf(tempobjectz, "%c", line[raza]);
                strcat(tempobject1z,tempobjectz);
                raza++;
            }
            //tempobject1
            i=atoi(tempobject1x);   // converting x corrdinate (which was read as string) to integer
            j=atoi(tempobject1y);   // converting y corrdinate (which was read as string) to integer
            k=atoi(tempobject1z);   // converting z corrdinate (which was read as string) to integer
            val=atoi(tempobject1);  // Checkpoint what is this for?
            //	printf("%d %d %d \n",i,j,k);
            bool ts=false;

            for (raza=0;raza<objectCount;raza++){
                if (f[i][j][k]==Object_datau[raza]. ObjectLabel)
                    ts=true;
            }
            if (!ts){
                if (f[i][j][k]!=0){
                    Object_datau[objectCount]. ObjectLabel=f[i][j][k];
                    objectCount=objectCount+1;
                }
            }


        }
        fclose(infile);  /* close the file prior to exiting the routine */
        //        fclose(writeOutFile1);
    }


    /** making a new folder to put all the text files in **/

    strcpy(outputDirectory,outputFilesPath);
    strcat(outputDirectory, "OutputDirectory/");
   

    strcat(outputDirectory_obj,outputFilesPath);
    strcat(outputDirectory_obj, "OutputDirectory/obj/");

    // printf(" final directory: %s\n ", outputDirectory);

    struct stat st = {0};

    if (stat(outputDirectory, &st) == -1) {
        mkdir(outputDirectory, 0700);
        printf("%\n Output Directory 'OutputDirectory' generated....\n");
    }

    if (stat(outputDirectory_obj, &st) == -1) {
        mkdir(outputDirectory_obj, 0700);
        printf("% Output Directory 'OutputDirectory/obj' generated....\n");
    }

    free_f3tensor (myObject, 0,nx*ny*nz+1,0,3,0,2);

    long mmax=0.0;

    for (i=0; i<5000; i++) {
        mmydata1[i].surf=0;
        mmydata1[i].vol=0;
        mmydata1[i].cont=0;

        mmydata1[i].POI_maj_min_x=0.0;
        mmydata1[i].POI_maj_min_y=0.0;
        mmydata1[i].POI_maj_min_z=0.0;

        mmydata1[i].minx=50000;
        mmydata1[i].maxx=0;
        mmydata1[i].miny=50000;
        mmydata1[i].maxy=0;
        mmydata1[i].minz=50000;
        mmydata1[i].maxz=0;
        mmydata1[i].majorx1=0;
        mmydata1[i].majorx2=0;
        mmydata1[i].minorx1=0;
        mmydata1[i].minorx2=0;
        mmydata1[i].middlex1=0;
        mmydata1[i].middlex2=0;

        mmydata1[i].majory1=0;
        mmydata1[i].majory2=0;
        mmydata1[i].minory1=0;
        mmydata1[i].minory2=0;
        mmydata1[i].middley1=0;
        mmydata1[i].middley2=0;

        mmydata1[i].majorz1=0;
        mmydata1[i].majorz2=0;
        mmydata1[i].minorz1=0;
        mmydata1[i].minorz2=0;
        mmydata1[i].middlez1=0;
        mmydata1[i].middlez2=0;
    }
    mmax=1;
    for (k=1; k<=nz; k++)
        for (j=1; j<=ny; j++)
            for (i=1; i<=nx; i++)
            {
                //              printf("%f\n",f[i][j][k]);
                u[i][j][k] = 0.0;//f[i][j][k];  // because we want to write the picked particles on output file. fxyz already has the LoG output
            }

    // Below we are writing each object to a file. This will save the processing time.

    int counter_x=0;
    long limit_counter=0;


    for (jj=0; jj<labelcounter; jj++) {
        long overall_vol=0.0;

        // printf("%f object name: in file %d \n",mmydata1[jj]. ObjectLabel, jj );
        sprintf(fileName1,"%ld",(jj+1));
        strcat(fileName1, ".txt");
        strcpy(fileNameWithPath1,fileName1);

        strcpy(fileNameWithPath1, outputDirectory_obj);
        strcat(fileNameWithPath1, fileName1);
        objectFile1 = fopen(&fileNameWithPath1, "w");   //change the path to global path
        if (objectFile1 == NULL)
        {
            printf("Error opening file! %s \n", fileNameWithPath);
            exit(1);
        }

        for (k=1; k<=nz; k++){
            for (i=1; i<=nx; i++){
                bool tm =false;
                long y_start=0;
                long y_end=0;
                for (j=1; j<=ny; j++){
                    if (f[i][j][k]!=0.0){

                        if (f[i][j][k]==labelarray[jj]){


                            if (!tm){
                                y_start=j;

                                tm=true;
                            }
                            y_end=j;

                            fprintf(objectFile1, "%ld %ld %ld\n",i,j,k);
                            limit_counter++;
                        }

                    }
                }
                if ((y_start!=0)&(y_end!=0))
                    overall_vol=overall_vol+(y_end-y_start)+1;
            }
        }
        fclose(objectFile1);

     //   if (limit_counter<1000000){
            bool   ts=true;
            mmydata1[counter_x].vol=overall_vol; // while writing file, we can find the volume.
            mmydata1[counter_x].index=jj;
            mmydata1[counter_x++].ObjectLabel=  labelarray[jj];
            ts=false;
       // }
    }
    labelcounter=counter_x;

    // Below we are opening the written object files again. Computational complextity will be linear as compared to prevoius version of code where we were not writing the files and doing the below mentioned steps using fxyz

    // We are calculating minor, middle and major axes in parallel as well.

    FILE *fileUserInput=NULL;
 
    int tempints=3000000;

    long *myx = (long *)malloc((tempints+1) * sizeof(*myx));      // allocate 50 ints
    long *myy = (long *)malloc((tempints+1) * sizeof(*myy));      // allocate 50 ints
    long *myz = (long *)malloc((tempints+1) * sizeof(*myz));      // allocate 50 ints

    for (size = 0; size <= 3000000; size++)       // valgrind
        myx[size] = 0;
    for (size = 0; size <= 3000000; size++)       // valgrind
        myy[size] = 0;
    for (size = 0; size <= 3000000; size++)       // valgrind
        myz[size] = 0;

     printf("Entered reading.......?\n");

    float  *pointOfIntersection     =   (float*)malloc(4 * sizeof(*pointOfIntersection));	//1D pointer


    for (size = 0; size <= 3; size++)       // valgrind
        pointOfIntersection[size] = 0;

    for (jj=0;jj<labelcounter;jj++){
        mmydata1[jj].POI_maj_min_x=0;
        mmydata1[jj].POI_maj_min_y=0;
        mmydata1[jj].POI_maj_min_z=0;
    }
    for (jj=0; jj<labelcounter; jj++) { // labelcounter=total number of objects

        // reading objects from text files start from here jth object

        if (mmydata1[jj].vol!=0){
            long myVolumeCounter=0;
            sprintf(fileName1,"%ld",(mmydata1[jj].index+1));
            //            printf("\t\tReading file  %ld (%ld) out of %d\n",fileCounter++,(mmydata1[jj].index+1), labelcounter);
            printf("\t\tReading object # %ld out of %d\n",(mmydata1[jj].index+1), labelcounter);

            strcat(fileName1, ".txt");
            strcpy(fileNameWithPath1,fileName1);
            strcpy(fileNameWithPath1, outputDirectory_obj);
            strcat(fileNameWithPath1, fileName1);

            fileUserInput      = fopen(&fileNameWithPath1, "rt");  /* open the file for reading text*/

            //   printf("Entered reading.......1\n");

            if(fileUserInput)
            {
                while((fgets(line, sizeof(line), fileUserInput)) != NULL)
                {
                    lineLength= strlen(line);

                    char tempobjectx [100]="";
                    char tempobject1x [100]="";

                    char tempobjecty [100]="";
                    char tempobject1y [100]="";

                    char tempobjectz [100]="";
                    char tempobject1z [100]="";

                    raza=0;
                    while ((line[raza]!=' ')){
                        sprintf(tempobjectx, "%c", line[raza]);
                        strcat(tempobject1x,tempobjectx);
                        raza++;
                    }

                    while ((line[raza]==' ')){

                        raza++;
                    }

                    while ((line[raza]!=' ')){
                        sprintf(tempobjecty, "%c", line[raza]);
                        strcat(tempobject1y,tempobjecty);
                        raza++;
                    }

                    while ((line[raza]==' ')){

                        raza++;
                    }
                    while (raza!=(lineLength-1)){
                        sprintf(tempobjectz, "%c", line[raza]);
                        strcat(tempobject1z,tempobjectz);
                        raza++;
                    }

                    i=atoi(tempobject1x);
                    j=atoi(tempobject1y);
                    k=atoi(tempobject1z);

                    if (i<mmydata1[jj].minx) {
                        mmydata1[jj].minx=i;
                    }
                    if (j<mmydata1[jj].miny) {
                        mmydata1[jj].miny=j;
                    }
                    if (i>mmydata1[jj].maxx) {
                        mmydata1[jj].maxx=i;
                    }
                    if (j>mmydata1[jj].maxy) {
                        mmydata1[jj].maxy=j;
                    }

                    if (k<mmydata1[jj].minz) {
                        mmydata1[jj].minz=k;
                    }
                    if (k>mmydata1[jj].maxz) {
                        mmydata1[jj].maxz=k;
                    }
                    //    printf("%d %d %d\n",i,j,k);
                    myx[myVolumeCounter]=i;//x
                    myy[myVolumeCounter]=j;//y
                    myz[myVolumeCounter++]=k;//z

                }
            }

            // jth object completely read and saved in myx,myy,myz
            //total number of voxels in myVolumeCounter

            double major    =   -100;
            double middle   =   -100;
            double minor    =   -100;

            //      printf("%d \n",Object_data[mRaza]. ObjectLabel);
            double dis=0.0;
            //we used two loops to compare first one with all next possible voxels ????
            for (i =0;i<myVolumeCounter-1;i++){

                for (j =i+1;j<myVolumeCounter;j++){

                    dis=sqrt(((myx[i]-myx[j])*(myx[i]-myx[j]))+((myy[i]-myy[j])*(myy[i]-myy[j]))+((myz[i]-myz[j])*(myz[i]-myz[j])));
                    if (dis>major){
                        major=dis;

                        mmydata1[jj].majorx1=myx[i];
                        mmydata1[jj].majorx2=myx[j];

                        mmydata1[jj].majory1=myy[i];
                        mmydata1[jj].majory2=myy[j];

                        mmydata1[jj].majorz1=myz[i];
                        mmydata1[jj].majorz2=myz[j];

                    }
                }
            }

            // Middle axis calculations
            double slp1=0.0;
            double vLineX1 = 0.0,vLineX2 = 0.0,vLineY1 = 0.0,vLineY2 = 0.0,vLineZ1 = 0.0,vLineZ2 = 0.0;
            double amag=1.0;
            double bmag=1.0;
            double angle1 =0.0, angle2=0.0;
            double val = 180.0 / PI;

            for (i =0;i<myVolumeCounter-1;i++){
                for (j =i+1;j<myVolumeCounter;j++){
                    // So, vector between two points is x2-x1,y2-y1,z2-z1
                    // now check whether vector generated by Major axis is perpendiculer to each other or not by their dot product. v.w=0
                    vLineX1= mmydata1[jj].majorx1 - mmydata1[jj].majorx2;
                    vLineY1= mmydata1[jj].majory1 - mmydata1[jj].majory2;
                    vLineZ1= mmydata1[jj].majorz1 - mmydata1[jj].majorz2;

                    vLineX2= myx[i]- myx[j];
                    vLineY2= myy[i]- myy[j];
                    vLineZ2= myz[i]- myz[j];

                    slp1= vLineX1 * vLineX2 + vLineY1 * vLineY2 + vLineZ1 * vLineZ2;
                    amag=sqrt(vLineX1*vLineX1+vLineY1*vLineY1+vLineZ1*vLineZ1);
                    bmag=sqrt(vLineX2*vLineX2+vLineY2*vLineY2+vLineZ2*vLineZ2);
                    angle1=acos(slp1/(amag*bmag))*val;
                    //                slp1=((mmydata1[jj].majorx2 - mmydata1[jj].majorx1)*(myx[j]- myx[i]))+ ((mmydata1[jj].majory2-mmydata1[jj].majory1)*(myy[j]-myy[i]))+ ((mmydata1[jj].majorz2-mmydata1[jj].majorz1)*(myz[j]-myz[i]));


                    pointOfIntersection[0]=-100;
                    pointOfIntersection[1]=-100;
                    pointOfIntersection[2]=-100;

                    //        threed_point_of_intersection(mmydata1[jj].majorx1, mmydata1[jj].majory1, mmydata1[jj].majorz1, mmydata1[jj].majorx2, mmydata1[jj].majory2, mmydata1[jj].majorz2, myx[j],myy[j],myz[j], myx[i],myy[i],myz[i],pointOfIntersection);

                    //		if (((angle1>89.5)&(angle1<90.5)&(pointOfIntersection[0]!=-100))){
                    if (((angle1>89.5)&(angle1<90.5))){

                        //      printf("xa=%f y=%f z=%f\n", pointOfIntersection[0],  pointOfIntersection[1],  pointOfIntersection[2]);

                        dis=sqrt(((myx[i]-myx[j])*(myx[i]-myx[j]))+((myy[i]-myy[j])*(myy[i]-myy[j]))+((myz[i]-myz[j])*(myz[i]-myz[j])));
                        //	if ((pointOfIntersection[0]!=-100)& (pointOfIntersection[1]!=-100)& (pointOfIntersection[2]!=-100))

                        if (dis > middle)
                        {
                            middle=dis;
                            mmydata1[jj].POI_maj_min_x=pointOfIntersection[0];
                            mmydata1[jj].POI_maj_min_y=pointOfIntersection[1];
                            mmydata1[jj].POI_maj_min_z=pointOfIntersection[2];

                            mmydata1[jj].middlex1=myx[i];
                            mmydata1[jj].middlex2=myx[j];


                            mmydata1[jj].middley1=myy[i];
                            mmydata1[jj].middley2=myy[j];


                            mmydata1[jj].middlez1=myz[i];
                            mmydata1[jj].middlez2=myz[j];

                            //printf("slope=%f midx1=%ld midx2=%ld midy1=%ld midy2=%ld midz1=%ld midz2=%ld \n", angle1, mmydata1[jj].middlex1,  mmydata1[jj].middlex2,  mmydata1[jj].middley1, mmydata1[jj].middley2, mmydata1[jj].middlez1, mmydata1[jj].middlez2);

                        }
                    }
                }
            }

            double slp2=0.0;
            amag=1.0;
            bmag=1.0;

            float midMinorx = 0.0, midMinory = 0.0, midMinorz = 0.0;
            float vectorLineMinorX = 0.0, vectorLineMinorY = 0.0, vectorLineMinorZ = 0.0;

            // Minor axis calculations

            for (i =0;i<myVolumeCounter-1;i++){
                for (j =i+1;j<myVolumeCounter;j++){

                    vLineX1= mmydata1[jj].majorx1 - mmydata1[jj].majorx2;
                    vLineY1= mmydata1[jj].majory1 - mmydata1[jj].majory2;
                    vLineZ1= mmydata1[jj].majorz1 - mmydata1[jj].majorz2;

                    vLineX2= myx[i]- myx[j];
                    vLineY2= myy[i]- myy[j];
                    vLineZ2= myz[i]- myz[j];


                    slp1= vLineX1 * vLineX2 + vLineY1 * vLineY2 + vLineZ1 * vLineZ2;

                    amag=sqrt(vLineX1*vLineX1+vLineY1*vLineY1+vLineZ1*vLineZ1);
                    bmag=sqrt(vLineX2*vLineX2+vLineY2*vLineY2+vLineZ2*vLineZ2);
                    angle1=acos(slp1/(amag*bmag))*val;

                    vLineX1= mmydata1[jj].middlex1 - mmydata1[jj].middlex2;
                    vLineY1= mmydata1[jj].middley1 - mmydata1[jj].middley2;
                    vLineZ1= mmydata1[jj].middlez1 - mmydata1[jj].middlez2;


                    slp2= vLineX1 * vLineX2 + vLineY1 * vLineY2 + vLineZ1 * vLineZ2;

                    amag=(sqrt((vLineX1*vLineX1)+(vLineY1*vLineY1)+(vLineZ1*vLineZ1)));
                    bmag=(sqrt((vLineX2*vLineX2)+(vLineY2*vLineY2)+(vLineZ2*vLineZ2)));
                    angle2=(acos(slp2/(amag*bmag))*val);

                    if ((((angle1>89.5)&(angle1<90.5)))
                        & (((angle2>89.5)&(angle2<90.5)))){

                        pointOfIntersection[0]=-100;
                        pointOfIntersection[1]=-100;
                        pointOfIntersection[2]=-100;

                        //	threed_point_of_intersection(mmydata1[jj].majorx1, mmydata1[jj].majory1, mmydata1[jj].majorz1, mmydata1[jj].majorx2, mmydata1[jj].majory2, mmydata1[jj].majorz2, myx[j],myy[j],myz[j], myx[i],myy[i],myz[i],pointOfIntersection);
                        //         if(pointOfIntersection[0]!=-100){
                        //  if ((pointOfIntersection[0]!=-100)& (pointOfIntersection[1]!=-100)& (pointOfIntersection[2]!=-100)){
                        pointOfIntersection[0]=-100;
                        pointOfIntersection[1]=-100;
                        pointOfIntersection[2]=-100;

                        //		threed_point_of_intersection(mmydata1[jj].middlex1, mmydata1[jj].middley1, mmydata1[jj].middlez1, mmydata1[jj].middlex2, mmydata1[jj].middley2, mmydata1[jj].middlez2, myx[j],myy[j],myz[j], myx[i],myy[i],myz[i],pointOfIntersection);

                        //  if ((pointOfIntersection[0]!=-100)& (pointOfIntersection[1]!=-100)& (pointOfIntersection[2]!=-100)){


                        dis=sqrt(((myx[i]-myx[j])*(myx[i]-myx[j]))+((myy[i]-myy[j])*(myy[i]-myy[j]))+((myz[i]-myz[j])*(myz[i]-myz[j])));


                        // printf("%ld %ld %ld ",coordinatesOfThirdAxis[0],coordinatesOfThirdAxis[1],coordinatesOfThirdAxis[1]);
                        if (dis>minor){
                            minor=dis;

                            mmydata1[jj].minorx1=myx[i];
                            mmydata1[jj].minorx2=myx[j];

                            mmydata1[jj].minory1=myy[i];
                            mmydata1[jj].minory2=myy[j];

                            mmydata1[jj].minorz1=myz[i];
                            mmydata1[jj].minorz2=myz[j];

                            // calculating midpoint of minor axis

                            midMinorx = ( mmydata1[jj].minorx1 + mmydata1[jj].minorx2 )/2;
                            midMinory = ( mmydata1[jj].minory1 + mmydata1[jj].minory2 )/2;
                            midMinorz = ( mmydata1[jj].minorz1 + mmydata1[jj].minorz2 )/2;

                            // calculting vector between midpoint of minor axis and point of intersection between major and middle axes.

                            vectorLineMinorX= mmydata1[jj].POI_maj_min_x - midMinorx;
                            vectorLineMinorY= mmydata1[jj].POI_maj_min_y - midMinory;
                            vectorLineMinorZ= mmydata1[jj].POI_maj_min_z - midMinorz;

                        }
                    }
                }
            }

            mmydata1[jj].surf=((myVolumeCounter-1)); //
            //     printf("%f Object Number ",mmydata1[jj].ObjectLabel);
            //     printf("    %ld Surface Area",mmydata1[jj].surf);
            //     printf("    %ld Volume \n",mmydata1[jj].vol);

            fclose(fileUserInput);

        }
        else {
            mmydata1[jj].minorx1=0.0;
            mmydata1[jj].minorx2=0.0;


            mmydata1[jj].minory1=0.0;
            mmydata1[jj].minory2=0.0;


            mmydata1[jj].minorz1=0.0;
            mmydata1[jj].minorz2=0.0;

            mmydata1[jj].majorx1=0.0;
            mmydata1[jj].majorx2=0.0;

            mmydata1[jj].majory1=0.0;
            mmydata1[jj].majory2=0.0;

            mmydata1[jj].majorz1=0.0;
            mmydata1[jj].majorz2=0.0;

            mmydata1[jj].middlex1=0.0;
            mmydata1[jj].middlex2=0.0;


            mmydata1[jj].middley1=0.0;
            mmydata1[jj].middley2=0.0;


            mmydata1[jj].middlez1=0.0;
            mmydata1[jj].middlez2=0.0;
        }
    }
    //	printf("pased 4\n");

    int mymiddlex1u[objectCount];
    int mymiddley1u[objectCount];
    int mymiddlez1u[objectCount];

    int mymiddlex2u[objectCount];
    int mymiddley2u[objectCount];
    int mymiddlez2u[objectCount];

    int mymajorx1u[objectCount];
    int mymajory1u[objectCount];
    int mymajorz1u[objectCount];

    int mymajorx2u[objectCount];
    int mymajory2u[objectCount];
    int mymajorz2u[objectCount];

    int myminorminor_x1u[objectCount];
    int myminorminor_x2u[objectCount];
    int myminorminor_y1u[objectCount];
    int myminorminor_y2u[objectCount];
    int myminorminor_z1u[objectCount];
    int myminorminor_z2u[objectCount];

    for (k=0; k<objectCount; k++) {
        for (jj=0; jj<labelcounter; jj++) {
            if (mmydata1[jj].ObjectLabel==Object_datau[k].ObjectLabel){

                mymiddlex1u[k]=mmydata1[jj].middlex1;
                mymiddlex2u[k]=mmydata1[jj].middlex2;

                mymiddley1u[k]=mmydata1[jj].middley1;
                mymiddley2u[k]=mmydata1[jj].middley2;

                mymiddlez1u[k]=mmydata1[jj].middlez1;
                mymiddlez2u[k]=mmydata1[jj].middlez2;

                mymajorx1u[k]=mmydata1[jj].majorx1;
                mymajorx2u[k]=mmydata1[jj].majorx2;

                mymajory1u[k]=mmydata1[jj].majory1;
                mymajory2u[k]=mmydata1[jj].majory2;

                mymajorz1u[k]=mmydata1[jj].majorz1;
                mymajorz2u[k]=mmydata1[jj].majorz2;

                myminorminor_x1u[k]=mmydata1[jj].minorx1;
                myminorminor_x2u[k]=mmydata1[jj].minorx2;

                myminorminor_y1u[k]=mmydata1[jj].minory1;
                myminorminor_y2u[k]=mmydata1[jj].minory2;

                myminorminor_z1u[k]=mmydata1[jj].minorz1;
                myminorminor_z2u[k]=mmydata1[jj].minorz2;

                Volume[k]=mmydata1[jj].vol;
                Surface[k]=mmydata1[jj].surf;

            }
        }
    }


    float sumInputMajorLength  = 0.0;
    float sumInputMiddleLength = 0.0;
    float sumInputMinorLength  = 0.0;
    float sumInputVolume       = 0.0;
    float sumInputSurface      = 0.0;

    float averageInputMajorLength  = 0.0;
    float averageInputMiddleLength = 0.0;
    float averageInputMinorLength  = 0.0;
    float averageInputVolume       = 0.0;
    float averageInputSurface      = 0.0;


    printf("\n User particles: major middle minor Volume Surface \n");

	long disregardedUser=0;
    int nondisregardedUser=0;


    for (raza=0;raza<objectCount;raza++){

        float mj1_u= sqrt((mymajorx1u[raza]-mymajorx2u[raza])*(mymajorx1u[raza]-mymajorx2u[raza])+(mymajory1u[raza]-mymajory2u[raza])*(mymajory1u[raza]-mymajory2u[raza])+(mymajorz1u[raza]-mymajorz2u[raza])*(mymajorz1u[raza]-mymajorz2u[raza]));

        float mdl1_u= sqrt((mymiddlex1u[raza]-mymiddlex2u[raza])*(mymiddlex1u[raza]-mymiddlex2u[raza])+(mymiddley1u[raza]-mymiddley2u[raza])*(mymiddley1u[raza]-mymiddley2u[raza])+(mymiddlez1u[raza]-mymiddlez2u[raza])*(mymiddlez1u[raza]-mymiddlez2u[raza]));

        float mnr1_u= sqrt((myminorminor_x1u[raza]-myminorminor_x2u[raza])*(myminorminor_x1u[raza]-myminorminor_x2u[raza])+(myminorminor_y1u[raza]-myminorminor_y2u[raza])*(myminorminor_y1u[raza]-myminorminor_y2u[raza])+(myminorminor_z1u[raza]-myminorminor_z2u[raza])*(myminorminor_z1u[raza]-myminorminor_z2u[raza]));

	if ((mj1_u<=0)||(mj1_u>32000)||(mdl1_u<=0)||(mdl1_u>32000)||(mnr1_u<=0)||(mnr1_u>32000)||(Volume[raza]<=0)|| (Surface[raza]<=0)){
	
	disregardedUser=disregardedUser+1;
        printf(" Disregarded User particles: %f %f %f %ld %ld \n", mj1_u, mdl1_u, mnr1_u, Volume[raza], Surface[raza]);

	} else {
        printf(" User particles: %f %f %f %ld %ld \n", mj1_u, mdl1_u, mnr1_u, Volume[raza], Surface[raza]);

        sumInputMajorLength     += mj1_u;
        sumInputMiddleLength    += mdl1_u;
        sumInputMinorLength     += mnr1_u;
        sumInputVolume          += Volume[raza];
        sumInputSurface         += Surface[raza];

        nondisregardedUser=nondisregardedUser+1;
        }

    }


    printf(" No of user particles selected by user :%ld \n", (objectCount+1));
    printf(" No of user particles disregarded due to low quality :%ld \n", disregardedUser);

    averageInputMajorLength    = sumInputMajorLength/(nondisregardedUser);
     averageInputMiddleLength   = sumInputMiddleLength/(nondisregardedUser);
     averageInputMinorLength    = sumInputMinorLength/(nondisregardedUser);
     averageInputVolume         = sumInputVolume/(nondisregardedUser);
     averageInputSurface        = sumInputSurface/(nondisregardedUser);

    printf("\n average particles: %f %f %f %f %f \n", averageInputMajorLength, averageInputMiddleLength, averageInputMinorLength, averageInputVolume, averageInputSurface);

//    printf("\n sigma particles: %f %f %f %f %f \n", sigmaInputMajorLength, sigmaInputMiddleLength, sigmaInputMinorLength, sigmaInputVolume, sigmaInputSurface);

    printf(" Step IV: generating output ...\n");

    float tolerance_rate = toleranceRate;   // tolerance loop incremented by 5

    for(myk=0;myk<=toleranceRate;myk+=2) {// tolerance rate loop
        totcontour=0;
        tolerance_rate = (float)myk;
        tol_min = 1-(tolerance_rate /100);  // (0.7-1.30 means 30% above and below),(0.9-1.1 means 10% above and below)
        tol_max= 1 + (tolerance_rate /100);


        tolRate_min_major       = averageInputMajorLength * tol_min;
        tolRate_max_major       = averageInputMajorLength * tol_max;

        tolRate_min_middle      = averageInputMiddleLength * tol_min;
        tolRate_max_middle      = averageInputMiddleLength * tol_max;

        tolRate_min_minor       = averageInputMinorLength * tol_min;
        tolRate_max_minor       = averageInputMinorLength * tol_max;

        tolRate_min_surfaceArea = averageInputSurface * tol_min;
        tolRate_max_surfaceArea = averageInputSurface * tol_max;

        tolRate_min_volume      = averageInputVolume * tol_min;
        tolRate_max_volume      = averageInputVolume * tol_max;

        printf(" min: %f %f %f %f %f %f \n", tolerance_rate, tolRate_min_major, tolRate_min_middle, tolRate_min_minor, tolRate_min_surfaceArea, tolRate_min_volume);
        printf(" max: %f %f %f %f %f %f \n", tolerance_rate, tolRate_max_major, tolRate_max_middle, tolRate_max_minor, tolRate_max_surfaceArea, tolRate_max_volume);


        for (jj=0; jj<labelcounter; jj++) {
            mmydata1[jj].cont=1;
        }

        for (jj=0; jj<labelcounter; jj++) {

            float mj1_o= sqrt((mmydata1[jj].majorx1-mmydata1[jj].majorx2)*(mmydata1[jj].majorx1-mmydata1[jj].majorx2)+(mmydata1[jj].majory1-mmydata1[jj].majory2)*(mmydata1[jj].majory1-mmydata1[jj].majory2)+(mmydata1[jj].majorz1-mmydata1[jj].majorz2)*(mmydata1[jj].majorz1-mmydata1[jj].majorz2)); // length of major axis of the objects detected by raza

            float mdl1_o= sqrt((mmydata1[jj].middlex1-mmydata1[jj].middlex2)*(mmydata1[jj].middlex1-mmydata1[jj].middlex2)+(mmydata1[jj].middley1-mmydata1[jj].middley2)*(mmydata1[jj].middley1-mmydata1[jj].middley2)+(mmydata1[jj].middlez1-mmydata1[jj].middlez2)*(mmydata1[jj].middlez1-mmydata1[jj].middlez2)); // length of middle axis of the objects detected by raza

            float mnr1_o= sqrt((mmydata1[jj].minorx1-mmydata1[jj].minorx2)*(mmydata1[jj].minorx1-mmydata1[jj].minorx2)+(mmydata1[jj].minory1-mmydata1[jj].minory2)*(mmydata1[jj].minory1-mmydata1[jj].minory2)+(mmydata1[jj].minorz1-mmydata1[jj].minorz2)*(mmydata1[jj].minorz1-mmydata1[jj].minorz2)); // length of minor axis of the objects detected by raza

                /*************** CHANGE TOLERANCE RATE BELOW ********************/

                // printf(" tolerance rate %f, %f\n ", tol_min, tol_max);
                if ((mj1_o>0) && (mdl1_o>0) && (mnr1_o>0) && (mmydata1[jj].surf >0) && (mmydata1[jj].vol>0) ){

                    printf(" particles: %f %f %f %ld %ld \n", mj1_o, mdl1_o, mnr1_o, mmydata1[jj].surf, mmydata1[jj].vol);

                    if (condition==1.0){    // all combined

                        if ((mj1_o >= tolRate_min_major)&(mj1_o<=tolRate_max_major)&(mdl1_o>=tolRate_min_middle)&(mdl1_o<=tolRate_max_middle)&(mnr1_o>=tolRate_min_minor)&(mnr1_o<=tolRate_max_minor)&(mmydata1[jj].vol>tolRate_min_volume)&(mmydata1[jj].vol<tolRate_max_volume)&(mmydata1[jj].surf>tolRate_min_surfaceArea)&(mmydata1[jj].surf<tolRate_max_surfaceArea)){

                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            //                        printf(" files are %ld\n ", ((mmydata1[jj].index+1)));
                            totcontour++;
                        }
                    }
                    else if (condition==2.0){   // 3d Axes together

                        if ((mj1_o>=tolRate_min_major)&(mj1_o<=tolRate_max_major)&(mdl1_o>=tolRate_min_middle)&(mdl1_o<=tolRate_max_middle)&(mnr1_o>=tolRate_min_minor)&(mnr1_o<=tolRate_max_minor)){

                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            totcontour++;
                        }
                    }
                    else if (condition==3.0){   //3D volume

                        if ((mmydata1[jj].vol>=tolRate_min_volume)&(mmydata1[jj].vol<=tolRate_max_volume)){
                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            totcontour++;
                        }
                    }
                    else if (condition==4.0){   //3D Surface Area

                        if ((mmydata1[jj].surf>=tolRate_min_surfaceArea)&(mmydata1[jj].surf<=tolRate_max_surfaceArea)){
                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            totcontour++;
                        }
                    }
                    else if (condition==5.0){   //major axis

                        if ((mj1_o>=tolRate_min_major)&(mj1_o<=tolRate_max_major)){

                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            totcontour++;
                        }
                    }
                    else if (condition==6.0){   //middle axis

                        if ((mdl1_o>=tolRate_min_middle)&(mdl1_o<=tolRate_max_middle)){

                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            totcontour++;
                        }
                    }
                    else if (condition==7.0){   // minor axis
                        if ((mnr1_o>=tolRate_min_minor)&(mnr1_o<=tolRate_max_minor)){
                            mmydata1[jj].cont=1214; //Aimah=12, Masoomeen=14;
                            totcontour++;
                        }
                    }
                }


//            }
        }

        /************************************************************************/
        /***** all information related to Maximum Tolerance rate will be saved in following output files  *******/
        /************************************************************************/


        printf(" /******************************/\n");
        printf(" /** Tolerance threshold = %d **/\n", (int)(myk));
        printf(" /******************************/\n");

        cont=2;
        obj_c=1;
        max_x1=0,min_x1=0,max_y1=0,min_y1=0,max_z1=0,min_z1=0;

        /************************************************************************/

        fileName1[500]="";
        sprintf(fileName1,"%d",(int)(myk));
        strcat(fileName1, "2_Final_subvolume_coordinates.txt");

        strcpy(outputDirectory,outputFilesPath);
        strcat(outputDirectory, "OutputDirectory/");
        strcpy(fileNameWithPath1, outputDirectory);
        strcat(fileNameWithPath1, fileName1);

        printf("\tFile 2: 'Final_subvolume_coordinates.txt'...\n");

        //        printf("Final_subvolume_coordinates path:  %s\n",fileNameWithPath1);
        //        printf("Final_subvolume_coordinates:  %s\n",fileName1);
        //

        writeOutFile2 = fopen(&fileNameWithPath1, "w");   //change the path to global path

        if (writeOutFile2 == NULL)
        {
            //        printf("Error opening 2_Final_subvolume_coordinates.txt file!\n");
            //        exit(1);
        }

        long averagex = 0.0,averagey = 0.0,averagez = 0.0;
        averagex=(max_x1-min_x1)/2;
        averagey=(max_y1-min_y1)/2;
        averagez=(max_z1-min_z1)/2;



        for (jj=0; jj<labelcounter; jj++) {
            if(mmydata1[jj].cont==1214){


                min_x1=    mmydata1[jj].minx;
                max_x1=    mmydata1[jj].maxx;
                min_y1=    mmydata1[jj].miny;
                max_y1=    mmydata1[jj].maxy;
                min_z1=    mmydata1[jj].minz;
                max_z1=    mmydata1[jj].maxz;
                //
                //            printf("%ld %ld \n",max_x1,min_x1);
                //            printf("%ld %ld \n",max_y1,min_y1);
                //            printf("%ld %ld \n",max_z1,min_z1);
                //

                if ((max_x1<=nx)&(min_x1>=0)&(max_y1<=ny)&(min_y1>=0)&(max_z1<=nz)&(min_z1>=0)){
                    fprintf(writeOutFile2, "%ld %ld %ld %ld %ld %ld\n",min_x1,max_x1,min_y1,max_y1,min_z1,max_z1);
                    printf("%ld %ld %ld %ld %ld %ld\n",min_x1,max_x1,min_y1,max_y1,min_z1,max_z1);
                    //
                }
            }
        }

        fclose(writeOutFile2);
        /************************************************************************/
        printf("\tFile 3: 'StructuralDetails.txt'...\n");

        long lengthMajorAxis = 0.0, lengthMiddleAxis = 0.0, lengthMinorAxis = 0.0;


        fileName1[500]="";
        sprintf(fileName1,"%d",(int)(myk));
        strcat(fileName1, "3_StructuralDetails.txt");

        strcpy(outputDirectory,outputFilesPath);
        strcat(outputDirectory, "OutputDirectory/");
        strcpy(fileNameWithPath1, outputDirectory);
        strcat(fileNameWithPath1, fileName1);



        writeOutFile3 = fopen(&fileNameWithPath1, "w");   //change the path to global path
        //            printf("Structural details %s",fileNameWithPath1);
        //            printf("Structural details %s\n",fileName1);
        //


        if (writeOutFile3 == NULL)
        {
            printf(" final directory: %s\n ", outputDirectory);
            printf("Error opening 3_StructuralDetails.txt file!\n");
            //		exit(1);
        }

        //	fprintf(writeOutFile3, "Object_No.	Length_MajorAxis	Length_MiddleAxis	Length_MinorAxis	Surface_Area	Volume\n");
        printf("Object_No.	Length_MajorAxis	Length_MiddleAxis	Length_MinorAxis	Surface_Area	Volume ObjectLabel\n");

        for (jj=0; jj<labelcounter; jj++) {
            if(mmydata1[jj].cont==1214){

                lengthMajorAxis=0;
                lengthMiddleAxis=0;
                lengthMinorAxis=0;

                lengthMajorAxis= sqrt(((mmydata1[jj].majorx2 - mmydata1[jj].majorx1)*(mmydata1[jj].majorx2 - mmydata1[jj].majorx1))+ ((mmydata1[jj].majory2 - mmydata1[jj].majory1)*(mmydata1[jj].majory2 - mmydata1[jj].majory1))+ ((mmydata1[jj].majorz2 - mmydata1[jj].majorz1)*(mmydata1[jj].majorz2 - mmydata1[jj].majorz1)));
                lengthMiddleAxis= sqrt(((mmydata1[jj].middlex2 - mmydata1[jj].middlex1)*(mmydata1[jj].middlex2 - mmydata1[jj].middlex1))+ ((mmydata1[jj].middley2 - mmydata1[jj].middley1)*(mmydata1[jj].middley2 - mmydata1[jj].middley1))+ ((mmydata1[jj].middlez2 - mmydata1[jj].middlez1)*(mmydata1[jj].middlez2 - mmydata1[jj].middlez1)));
                lengthMinorAxis= sqrt(((mmydata1[jj].minorx2 - mmydata1[jj].minorx1)*(mmydata1[jj].minorx2 - mmydata1[jj].minorx1))+ ((mmydata1[jj].minory2 - mmydata1[jj].minory1)*(mmydata1[jj].minory2 - mmydata1[jj].minory1))+ ((mmydata1[jj].minorz2 - mmydata1[jj].minorz1)*(mmydata1[jj].minorz2 - mmydata1[jj].minorz1)));


                fprintf(writeOutFile3, "%d\t%ld\t%ld\t%ld\t%ld\t%ld\n",obj_c, lengthMajorAxis, lengthMiddleAxis, lengthMinorAxis,mmydata1[jj].surf, mmydata1[jj].vol);

                printf("%d\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\n",obj_c, lengthMajorAxis, lengthMiddleAxis, lengthMinorAxis,mmydata1[jj].surf, mmydata1[jj].vol,mmydata1[jj].ObjectLabel);

                obj_c++;

            }

        }

        fclose(writeOutFile3);

        tolerance_rate = (1 - tol_min) * 100;
        //
        //        fprintf(writeOutFile3, "Total no objects detected %d out of %d \n", totcontour, labelcounter);
        //        fprintf(writeOutFile3, "Tolerance Rate %0.2f percent \n", tolerance_rate);
        //        fprintf(writeOutFile3, "Total Processing Time:-----------------------\n");

        printf("Total no objects detected %d out of %d \n", totcontour, labelcounter);
        printf("Tolerance Rate %0.2f percent \n", tolerance_rate);
        printf(" ---------------------------------------------------------------------\n");

    }

    /************************************************************************/

    //      printf("\n*** OUTPUT FILE ***\nFile '1_model_file_major_minor.txt' starts here...\n");

    printf("\tFile 1: 'model_file_major_minor.txt'...\n");

    fileName1[500]="";
    sprintf(fileName1,"%d",(int)(myk)); //myk is tolerance rate
    strcat(fileName1, "1_model_file_major_minor.txt");

    strcpy(outputDirectory,outputFilesPath);
    strcat(outputDirectory, "OutputDirectory/");
    strcpy(fileNameWithPath1, outputDirectory);
    strcat(fileNameWithPath1, fileName1);

    //        printf("The minor major file write path %s \n",fileNameWithPath1);
    writeOutFile1 = fopen(&fileNameWithPath1, "w");   //change the path to global path

    if (writeOutFile1 == NULL)
    {
        printf(" final directory: %s\n ", outputDirectory);
        printf("Error opening 1_model_file_major_minor.txt file!\n");
        exit(1);
    }

    for (jj=(labelcounter-1); jj>=0; jj--) {
        if(mmydata1[jj].cont==1214){

            if ((mmydata1[jj].maxx-mmydata1[jj].minx) > (max_x1-min_x1)){
                min_x1=mmydata1[jj].minx;
                max_x1=mmydata1[jj].maxx;
            }
            if ((mmydata1[jj].maxy-mmydata1[jj].miny) > (max_y1-min_y1)){
                min_y1=mmydata1[jj].miny;
                max_y1=mmydata1[jj].maxy;
            }
            if ((mmydata1[jj].maxz-mmydata1[jj].minz) > (max_z1-min_z1)){
                min_z1=mmydata1[jj].minz;
                max_z1=mmydata1[jj].maxz;
            }

            /*
             Following lines are printing in following order on text file '1_model_file_major_minor.txt'

             line 1: first point on major axis of current object
             line 2: second point on major axis of current object
             line 3: first point on middle axis of current object
             line 4: second point on middle axis of current object
             line 5: first point on minor axis of current object
             line 6: second point on minor axis of current object
             */

            //           for (razu=mmydata1[jj].majorz1; razu< mmydata1[jj].majorz2; razu++) {

            fprintf(writeOutFile1, "%d %d %ld %ld %ld\n",obj_c,cont,mmydata1[jj].majorx1,mmydata1[jj].majory1,mmydata1[jj].majorz1);
            fprintf(writeOutFile1, "%d %d %ld %ld %ld\n",obj_c,cont,mmydata1[jj].majorx2,mmydata1[jj].majory2,mmydata1[jj].majorz2);

            //          }
            obj_c++;

            //          for (razu=mmydata1[jj].middlez1; razu< mmydata1[jj].middlez2; razu++) {

            fprintf(writeOutFile1, "%d %d %ld %ld %ld\n",obj_c,cont,mmydata1[jj].middlex1,mmydata1[jj].middley1,mmydata1[jj].middlez1);
            fprintf(writeOutFile1, "%d %d %ld %ld %ld\n",obj_c,cont,mmydata1[jj].middlex2,mmydata1[jj].middley2,mmydata1[jj].middlez2);

            //            }
            obj_c++;

            //            for (razu=mmydata1[jj].minorz1; razu< mmydata1[jj].minorz2; razu++) {

            fprintf(writeOutFile1, "%d %d %ld %ld %ld\n",obj_c,cont,mmydata1[jj].minorx1,mmydata1[jj].minory1,mmydata1[jj].minorz1);
            fprintf(writeOutFile1, "%d %d %ld %ld %ld\n",obj_c,cont,mmydata1[jj].minorx2,mmydata1[jj].minory2,mmydata1[jj].minorz2);
            //          }
            //            printf("Code reached Here \n");

            obj_c++;

            totcontour++;
            //		printf("final:%ld\n",mmydata1[jj].surf);
            //		printf("final:%ld\n",mmydata1[jj].vol);

            sprintf(fileName1,"%ld",(mmydata1[jj].index+1));
            strcat(fileName1, ".txt");
            strcpy(fileNameWithPath1,fileName1);
            strcpy(fileNameWithPath1, outputDirectory_obj);
            strcat(fileNameWithPath1, fileName1);

            fileUserInput      = fopen (&fileNameWithPath1, "rt");  /* open the file for reading text*/

            //          printf("Entered reading.......1\n");

            if(fileUserInput)
            {
                //				printf("inside..");
                while((fgets(line, sizeof(line), fileUserInput)) != NULL)
                {
                    lineLength= strlen(line);

                    char tempobjectx [100]="";
                    char tempobject1x [100]="";

                    char tempobjecty [100]="";
                    char tempobject1y [100]="";

                    char tempobjectz [100]="";
                    char tempobject1z [100]="";

                    raza=0;
                    while ((line[raza]!=' ')){
                        sprintf(tempobjectx, "%c", line[raza]);
                        strcat(tempobject1x,tempobjectx);
                        raza++;
                    }

                    while ((line[raza]==' ')){

                        raza++;
                    }

                    while ((line[raza]!=' ')){
                        sprintf(tempobjecty, "%c", line[raza]);
                        strcat(tempobject1y,tempobjecty);
                        raza++;
                    }

                    while ((line[raza]==' ')){

                        raza++;
                    }
                    while (raza!=(lineLength-1)){
                        sprintf(tempobjectz, "%c", line[raza]);
                        strcat(tempobject1z,tempobjectz);
                        raza++;
                    }

                    i=atoi(tempobject1x);
                    j=atoi(tempobject1y);
                    k=atoi(tempobject1z);

                    u[i][j][k]=255;
                }
            }
            fclose(fileUserInput);

        }

    }

    fclose(writeOutFile1);



    for (z=1; z<nz-1; z++) {
        u[1][1][z]=255;
        u[1][2][z]=255;
    }


    //printf("\n Total number of contours using z-crosing: %d \n", labelcounter);
    //    printf("\n Total number of contours using Particle picking: %d \n", totcontour);
    //

    free_f3tensor (f, 0,nx+1,0,ny+1,0,nz+1);
    free(labelarray);
    free(array);
    free(outputDirectory );
    free(outputDirectory_obj );
    free(myx );
    free(myy );
    free(myz );
    free(pointOfIntersection );
    return;
}

#define STRING_MAX 65536

int main (int argc, char **argv)
{

    clock_t begin, end;
    double time_spent;

   // printf("\t I just started in main....\n");

    begin=clock();

    MrcHeader header;
    float  ***u;                   /* image */
    int   i = 0, j = 0, k = 0, p = 0;             /* loop variables */
    int   nx = 0, ny = 0, nz = 0;             /* image size in x, y direction */
    FILE   *fp_infile, *fp_outfile = NULL;    /* input file, output file */
//    FILE   *fp_logoutfile = NULL;    /* input file, output file */
//    FILE   *fp_zcrossoutfile = NULL;    /* input file, output file */

    float  ht = 0.1f;              /* time step size */
    //float  zcvalue=0;
    //float  sigma = 1.0;             /* kernel size */
    float  max = 0.0, min = 0.0;               /* largest, smallest grey value */
    float  mean = 0.0;                   /* average grey value */
    float  vari = 0.0;
    int    pmax = 1;              /* largest iteration number */

    float sigma=1.0;//1.65;
    // it's better if user put value between 1 to 6
    // 1.2 gives all values positivve and 3.2 gives all values negative
    //2.2 all negativ
    float zcvalue= 0;
    float condition=1; //   1-> All (2,3 4) , 2-> 3D Axes, 3-> 3D volume, 4-> Surface Area
    float toleranceRate=10; //  10 %

    /* variance */
    Islice *sl;
    int sliceMode;
    char *progname = imodProgName(argv[0]);
    char *endptr;
    int nWrite = 0;
    int doWrite,doWritelog, doWritezcross, iarg, nzout, kst, knd, writeArg;
    int* writeList;
    char outFile[STRING_MAX]="";
    int oneSlice = 0;
    int outMode = -1;
    float minout = 0.0, maxout = 0.0, sumout = 0.0;

    struct tm *startlocalTime;
    struct tm *finishlocalTime;
    //	struct tm *totalProcessingTime;
    int hours = 0, minutes = 0, seconds = 0;
    int tmp = 0;


    time_t t;

    if (argc < 2)
        usage(progname, ht, sigma, zcvalue); //if total number of arguments are less than 2, then it will show decription of filter

    setExitPrefix("ERROR: 3d_log_zeroCrossings -");
    t = time(NULL);
    startlocalTime = localtime(&t);
    printf ("\nProgram %s\n", progname);
    printf (" started at: %s\n", asctime(startlocalTime));

    if (argc < 3) {
        printf("ERROR: %s - incorrect number of input arguments\n", progname);
        usage(progname,  ht , sigma, zcvalue);
        exit(3);
    }



    for (iarg = 1; iarg < argc; iarg++){
        if (argv[iarg][0] == '-'){
            switch (argv[iarg][1]){

                case 'z':
                    zcvalue = strtod(argv[++iarg], &endptr);
                    printf(" Zcross value:     %f\n", zcvalue);

                    testNumericEntry(endptr, argv[iarg], "-z");

                    break;
                case 's':
                    sigma = strtod(argv[++iarg], &endptr);
                    printf(" Sigma:            %f\n", sigma);

                    testNumericEntry(endptr, argv[iarg], "-s");
                    break;

                case 'c':
                    condition = strtod(argv[++iarg], &endptr);
                    if (condition==1) {

                        tmp= (int) condition;

                        printf(" Search condition: %d - All major, middle and minor axes, volume, surface area combined \n", tmp);
                    } else if(condition==2) {

                        tmp= (int) condition;
                        printf(" Search condition: %d - Major, middle and minor combined \n", tmp);
                    } else if(condition==3) {

                        tmp= (int) condition;
                        printf(" Search condition: %d - 3D Volume \n", tmp);

                    } else if(condition==4) {

                        tmp= (int) condition;
                        printf(" Search condition: %d - 3D Surface area \n", tmp);

                    } else if(condition==5) {

                        tmp= (int) condition;
                        printf(" Search condition: %d - Major axis \n", tmp);

                    } else if(condition==6) {

                        tmp= (int) condition;
                        printf(" Search condition: %d - Middle axis \n", tmp);

                    } else if(condition==7) {

                        tmp= (int) condition;
                        printf(" Search condition: %d - Minor axis \n", tmp);
                    }
                    testNumericEntry(endptr, argv[iarg], "-c");

                    break;

                case 'r':

                    toleranceRate = strtod(argv[++iarg], &endptr);
                    printf(" Tolerance rate(main):          %f\n", toleranceRate);

                    testNumericEntry(endptr, argv[iarg], "-r");

                    break;

                case 'o':
                    oneSlice = strtol(argv[++iarg], &endptr, 10);
                    testNumericEntry(endptr, argv[iarg], "-o");
                    break;
                case 'm':
                    outMode = strtol(argv[++iarg], &endptr, 10);
                    testNumericEntry(endptr, argv[iarg], "-m");
                    if (sliceModeIfReal(outMode) < 0)
                        exitError("Output mode %d not allowed\n", outMode);
                    break;
                case 't':
                    ht = strtod(argv[++iarg], &endptr);
                    testNumericEntry(endptr, argv[iarg], "-t");
                    break;
                case 'i':
                    writeList = parselist(argv[++iarg], &nWrite);
                    writeArg = iarg;
                    break;
                case 'P':
                    pidToStderr();
                    break;
                default:
                    exitError("Invalid option %s", argv[iarg]);
                    break;
            }

        }else{
            break;
        }
    }

    // Override pmax with maximum to write
    if (nWrite > 0) {
        pmax = 0;
        for (i = 0; i < nWrite; i++)
            pmax = B3DMAX(pmax, writeList[i]);
    }

    if (iarg != argc - 2)
        exitError("Command line should end with input and output files");

    //	printf("input file:          %s\n", argv[iarg]);
    //	printf("output file:         %s\n", argv[iarg + 1]);

    /* ---- read input directory path ---- */

    const char *srcFullPath=argv[iarg];
    char *srcInputFile =argv[iarg];
    int src_path_len=strlen(srcFullPath);

    char *ssc = NULL;
    int l = 0;
    int size =0;
    ssc = strstr(srcInputFile, "/");	//char *strstr(const char *haystack, const char *needle) function finds the first occurrence of the substring needle in the string haystack. The terminating '\0' characters are not compared.
//    printf ("ssc: %s\n", ssc);
    do{
        l = strlen(ssc) + 1;
//        printf ("l: %d\n", l);
        srcInputFile = &srcInputFile[strlen(srcInputFile)-l+2];
//        printf ("strlen(srcInputFile)-l+2: %d\n", strlen(srcInputFile)-l+2);
//        printf ("srcInputFile: %s\n", srcInputFile);
        ssc = strstr(srcInputFile, "/");
//        printf ("ssc: %s\n", ssc);
    }while(ssc);

    int remove_string_len=strlen(srcInputFile);
//    printf ("remove_string_len: %d\n", remove_string_len);
    int result_string_len= src_path_len-remove_string_len;
//    printf ("result_string_len: %d\n", result_string_len);
    char *outputFilesPath = ( char * )malloc( (result_string_len +1) * sizeof( char ) );    //added 1 because of //valgrind initialization

    for (size = 0; size < result_string_len+1; size++)      //valgrind initialization
        outputFilesPath[size] = NULL;

    strncpy(outputFilesPath, srcFullPath, result_string_len);


    /* ---- read input image ---- */
    if ((fp_infile = fopen (argv[iarg], "rb")) == 0)
        exitError("Could not open input file %s", argv[iarg]);

    /* read header */
    if (mrc_head_read(fp_infile, &header))
        exitError("Reading header of input file %s", argv[iarg]);

    // Check if it is the correct data type and set slice type
    sliceMode = sliceModeIfReal(header.mode);
    if (sliceMode < 0)
        exitError("File mode is %d; only byte, short, integer allowed",
                  header.mode);

    nx = header.nx;
    ny = header.ny;
    nz = header.nz;
    if (oneSlice < 1 || oneSlice > nz)
        oneSlice = 0;

    /* allocate storage */
    u=f3tensor( 0,nx+1,0,ny+1,0,nz+1);

    /* read image data */
    for (k=1; k<=nz; k++) {

        // Create a slice and read into it
        sl = sliceCreate(nx, ny, sliceMode);
        if (!sl)
            exitError("Creating slice for input");
        if (mrc_read_slice(sl->data.b, fp_infile, &header, k - 1, 'Z'))
            exitError("Reading slice %d", k);

        // Convert slice to floats
        if (sliceMode != SLICE_MODE_FLOAT)
            if (sliceNewMode(sl, SLICE_MODE_FLOAT) < 0)
                exitError("Converting slice to float");

        // Copy data into array
        for (j=0; j<ny; j++)
            for (i=0; i<nx; i++)
                u[i+1][j+1][k] = sl->data.f[i + j * nx];

        sliceFree(sl);
    }
    fclose(fp_infile);

    /* ---- Image ---- */
    analyse (u, nx, ny, nz, &min, &max, &mean, &vari);

    printf(" Input file name:  %s\n", srcInputFile);
    printf(" Input file path:  %s\n", outputFilesPath);   //here outputFilesPath contains the path of the input directory
    printf(" Dimensions:       %d x %d x %d\n", nx, ny, nz);
    printf(" Minimum:          %1.10f \n", min);
    printf(" Maximum:          %1.10f \n", max);
    printf(" Mean:             %1.10f \n", mean);
    printf(" Variance:         %1.10f \n\n", vari);

    // Take care of header
    if (oneSlice && nWrite)
        sprintf(outFile, "%s: Z %d, iter %s", progname, oneSlice,
                argv[writeArg]);
    else
        sprintf(outFile, "%s: Edge-enhancing laplacian of Gaussian z Crossings", progname);
    mrc_head_label(&header, outFile);

//    printf(" 1 \n");

    // Adjust output mode if valid entry made (it was tested on arg processing)
    if (outMode >= 0) {
        sliceMode = sliceModeIfReal(outMode);
        header.mode = outMode;
    }

//    printf(" 2 \n");

    // Fix things in header for an output file.  You have to set header size
    // not just set next to 0
    mrcInitOutputHeader(&header);

//    printf(" 3 \n");

    /* ---- process image ---- */

    nzout = 0;
    minout = 1.e30;
    maxout = -minout;
    sumout = 0;

//    printf(" 4 \n");

    for (p=1; p<=pmax; p++) {
        /* perform one iteration */

        threed_log_NonSeparable(nx, ny, nz, sigma, u);

        threed_zcrossings(ht, nx, ny, nz, 1.0, 1.0, 1.0, zcvalue, u);

	threed_ObjectCenter (u, nx, ny, nz, outputFilesPath, condition, toleranceRate );

        /* check minimum, maximum, mean, variance */
        analyse (u, nx, ny, nz, &min, &max, &mean, &vari);
        printf("minimum:       %1.10f \n", min);
        printf("maximum:       %1.10f \n", max);
        printf("mean:          %1.10f \n", mean);
        printf("variance:      %1.10f \n\n", vari);

        // Write data if it is an iteration on list or the last iteration
        doWrite = 0;
        for (i = 0; i < nWrite; i++)
            if (p == writeList[i])
                doWrite = 1;
        if (doWrite|| p == pmax) {

            header.amean = mean;
            header.amin = min;
            header.amax = max;

            strncpy(outFile, argv[iarg + 1], STRING_MAX - 10);
            outFile[STRING_MAX - 10] = 0;
            if (nWrite && !oneSlice)
                sprintf(&outFile[strlen(outFile)], "-%03d", p);

            /* open output file if not open yet */
            if (!fp_outfile) {
                if (!getenv("IMOD_NO_IMAGE_BACKUP") && imodBackupFile(outFile))
                    fprintf(stderr, "WARNING: error renaming existing %s to %s~",
                            outFile, outFile);

                if ((fp_outfile = fopen (outFile, "wb")) == 0)
                    exitError("Could not open output file %s", outFile);
            }

            /* write image data and close file */
            kst = oneSlice ? oneSlice : 1;
            knd = oneSlice ? oneSlice : nz;
            for (k=kst; k<=knd; k++) {
                // Create a slice and copy into it
                sl = sliceCreate(nx, ny, SLICE_MODE_FLOAT);
                if (!sl)
                    exitError("Creating slice for output");
                for (j=0; j<ny; j++)
                    for (i=0; i<nx; i++)
                        sl->data.f[i + j * nx] = u[i+1][j+1][k];

                // Convert if necessary and write slice
                if (sliceMode != SLICE_MODE_FLOAT)
                    if (sliceNewMode(sl, sliceMode) < 0)
                        exitError("Converting slice to short");
                if (mrc_write_slice(sl->data.b, fp_outfile, &header,
                                    oneSlice ? nzout : k - 1, 'Z'))
                    exitError("Writing slice %d", k);

                // If doing one slice, accumulate mmm
                if (oneSlice) {
                    sliceMMM(sl);
                    minout = B3DMIN(minout, sl->min);
                    maxout = B3DMAX(maxout, sl->max);
                    sumout += sl->mean;
                }
                sliceFree(sl);
            }

            nzout++;

            // Close the file if not doing one slice or end of run
            if (!oneSlice || p == pmax) {

                // Adjust size and mmm if did one slice
                if (oneSlice) {
                    header.nz = nzout;
                    header.mz = (header.mz * nzout) / nz;
                    header.zlen = (header.zlen * nzout) / nz;
                    header.amin = minout;
                    header.amax = maxout;
                    header.amean = sumout / nzout;
                }

                // Write the MRC header.
                if (mrc_head_write(fp_outfile, &header))
                    exitError("Writing header");

                fclose(fp_outfile);
                printf("output image %s successfully written\n\n", outFile);
                fp_outfile = NULL;
            }
        }
        fflush(stdout);
    } /* for */

    t				= time(NULL);
    finishlocalTime = localtime(&t);
    printf ("program finished at: %s\n", asctime(finishlocalTime));

    end				= clock();
    time_spent		= (double)(end-begin)/ CLOCKS_PER_SEC;
    printf ("Total Processing time: %f\n", time_spent);
    
    hours			= time_spent/3600;
    minutes			= ((int)time_spent/60)%60;
    seconds			= (int)time_spent %60;	
    printf ("Total processing time: %02d hours,%02d minutes, %02d seconds\n", hours, minutes, seconds);  
    /* ---- disallocate storage ---- */
    free(outputFilesPath);
    free_f3tensor (u, 0,nx+1,0,ny+1,0,nz+1);
    exit(0);
    
}



