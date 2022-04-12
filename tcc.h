/*  tc-cylinder: high-accurate terrain corrections
 *
 *  Copyright (C) 2022 S. Olgun and A. Ustun

 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<math.h>

#define R2D (180./M_PI)
#define R2S ((180./M_PI)*3600)
#define aE  (6378137.)
#define fE  (1/298.257222101)
#define ee1 (2*fE-fE*fE)
#define ee2 (ee1/(1.-ee1))
#define cE  (aE*(1.-fE)/(1.-ee1))
#define F2RND(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))
#define SI1N(f) ((double) (f > 1.0 ? (1.0) : (f+0.0)))

typedef struct 
{
    int zN;          /* number of zones */
    int cN;          /* total number of compartments */
    int    *compN;   /* total number of points inside each compertment */ 
    int    *sectorN; /* compartment numbers in zones */
    double *mean;   /* mean topographic height or density inside each compartement */ 
    double *zoneR;   /* raidius of zones */
    double tc;
    double rms;
    double bb;
    double ba;
    double bt;
} CylinderTmpl;

typedef struct 
{
    int nx;          /* number of x grids */
    int ny;          /* number of y grids */
    double xinc;     /* x increement */ 
    double yinc;     /* y increement */ 
    double xTL;      /* Top-Left real x coordinate */ 
    double yTL;      /* Top-Left real y coordinate */ 
    short int *H;    /* z coordinates */
} hgtGrid;

typedef struct 
{
    int nx;          /* number of x grids */
    int ny;          /* number of y grids */
    double xinc;     /* x increement */ 
    double yinc;     /* y increement */ 
    double xTL;      /* Top-Left real x coordinate */ 
    double yTL;      /* Top-Left real y coordinate */ 
    float *H;        /* z coordinates */
} Grid2D;

typedef struct
{
	int N;     /* zone number */
	double *R; /* radious of zones */
} ZoneTmpl;

typedef struct
{
    int idx;
	int iTL;
    int nx;
    int ny;
	double xTL;
	double yTL;
} SubFrame;
void help();
void initHammerTmpl(CylinderTmpl *);
void initHBowieNewell(CylinderTmpl *);
void initHBowieRapp(CylinderTmpl *);
void initTmpl(CylinderTmpl *);
void DestroyTemplate(CylinderTmpl *);
double RGauss(double);
int sphdirsol(double,double, double,double *,double *);
int sphinvsol(double,double, double,double *,double *);
#pragma omp declare simd
void tcfromTmpl(CylinderTmpl *,double,int);
void delTmpl(CylinderTmpl *);
int readhgt(Grid2D *,char *,double,double,double,double,double,double);
void InfoGrid(const Grid2D *);
int minmaxGrid(const Grid2D *g,float *,float *);
void setZone(ZoneTmpl *,int,double []);
void defZone(ZoneTmpl *);
SubFrame *subGridFrame(const Grid2D*,double,double,double);
void Grid2Tmpl(const Grid2D *,SubFrame *,double,double,int,int,CylinderTmpl *);
