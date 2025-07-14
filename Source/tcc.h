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
#include"gmt.h"
#include<time.h>

#define R2D (180/M_PI)
#define R2S ((180/M_PI)*3600)
#define R0 (6371000.0)
#define aE  (6378137.)
#define fE  (1/298.257222101)
#define ee1 (2*fE-fE*fE)
#define ee2 (ee1/(1.-ee1))
#define cE  (aE/(1.-fE))
#define F2RND(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))
#define SI1N(f) ((double) (f > 1.0 ? (1.0) : (f+0.0)))
#define NG  (6.67430e-6) /* m**3 kg**-1 s**-2*/
#define TMD (2670) /* kg m**-3*/
#define SWD (1024) /* kg m**-3*/
#define ATD (1.20) /* kg m**-3*/
#define TPG (2.0*M_PI*NG*TMD)  /* Topographical mGal/m */
#define SPG (2.0*M_PI*NG*SWD)  /* Marine mGal/m */
#define APG (2.0*M_PI*NG*ATD)  /* Atmospherical mGal/m */

typedef struct
{
    double lat;
    double lon;
    double H;
} GeodPos;

typedef struct 
{
    int zN;          /* number of zones */
    int cN;          /* total number of compartments */
    int    *compN;   /* total number of points inside each compertment */ 
    int    *sectorN; /* compartment numbers in zones */
    double *mean;    /* mean topographic height or density inside each compartement */ 
    double *zoneR;   /* radius of zones */
    double *zmR;     /* mean radius of each zone */
    double *zDR;     /* outer and inner radius diff (zone depth) */
    double rms;
    double mH;       /* Average of template heights representing all compartments */
    double bb;       /* Boullard B (spherical correction) */
    double ba;       /* Boullard A (planar plate) correction */
    double tc;       /* Boullard C (terrain correction) */
    double bt;       /* Boullard totatl correction */
} CylinderTmpl;

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
void initTmplTest(CylinderTmpl *);
void initTmpl(CylinderTmpl *);
void initTmplnew(CylinderTmpl *);
void initComp(CylinderTmpl *);
void Destroy_Template(CylinderTmpl *);
double RGauss(double);
int sphdirsol(double,double, double,double *,double *);
int sphinvsol(double,double, double,double *,double *);
void lltoazeq(const double,const double, double *,double *);
void tcfromTmpl(CylinderTmpl *,GeodPos *,int);
void delTmpl(CylinderTmpl *);
void InfoGrid(const Grid2D *);
int minmaxGrid(const Grid2D *g,float *,float *);
void Grid2Tmpl(struct GMT_GRID *,GeodPos *,const int,const int,CylinderTmpl *,int);
SubFrame *subGridFrame(struct GMT_GRID *,GeodPos *,double);
char* ctimer();
void print_info_GRID(struct GMT_GRID *);
struct GMT_GRID* read_GMT_Grid2D(void *,char *,int);
struct GMT_GRID* GMT_Grid_Resamp(void *,struct GMT_GRID *,GeodPos *,double,char *,int);

