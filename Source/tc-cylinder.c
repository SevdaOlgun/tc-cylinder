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

 ***********************************************************************************************
 * tc-cylinder is a C language based open source program which calculate 
 * high-accurate terrain corrections. It also provides 
 * complete Bouguer reductions with spherical form. 

 * The user can specify arbitrary DEM resolutions, the inner region radius,
 * the interpolation region and segmentation type properties. 
 * Please use help argument (-h) before the usage.

 * GMT 6.0.0 (General Mapping Tools) must be installed due to the use of its some tools. 
 * GMT's library should be used when compiling (an example of compiling: 
 * icc -o tc-cylinder tc-cylinder.c tcc.c -lm -I/usr/include/gmt/ -L/pub/opt/gmt/lib/ -lgmt).
 *
 ************************************************************************************************* 
 *  Authers: Sevda Olgun (Kocaeli University),
 *           Aydin Ustun (Hacettepe University)
 *  e-mail : sevda.olgun@kocaeli.edu.tr
 *  Created: March 2022 	    v1.0 Program Created
 *  Created: October 2024 	v2.0 Program Modified
 */

#include<unistd.h>
#include<omp.h>
#include "gmt.h"
#include "tcc.h"
#include<stdio.h>
#define OPTIONS "G:I:T:Z:r:p:hV"

FILE* read_data(int,char *,int);
int main (int argc, char *argv[])
{
    /* command line option flags */
    int Gflag = 0;
    int Iflag = 0;
    int Vflag = 0;
    int rflag = 0;
    /* arguments and variables */
    char Topt='S';  /* default template type S (detailed segmentation) */
    char *gfn;      /* Coarse (outter) grid file name */
    char *ifn;      /* Fine (inner) grid file name */
    char *ires;     /* interpolation resolution of innermost grid */
    int iZ = 56;    /* Inner grid radius index (see templates) */
    int irZ=3;      /* Interpolation zone*/
    char  line[256];
    char  c;

    while((c=getopt(argc,argv,OPTIONS))!=-1)
        switch(c)
        {
            case 'G':
                Gflag = 1;
                gfn=optarg;
                break;
            case 'I':
                Iflag=1;
                ifn=optarg;
                break;
            case 'T':
                Topt=optarg[0];
                break;
            case 'Z':
                iZ=atoi(optarg);
                break;
            case 'r':
                rflag=1;
                ires=optarg;
                break;
            case 'p':
                irZ=atoi(optarg);
                break;
            case 'V':
                Vflag=1;
                break;
            case 'h':
                help();
                exit(EXIT_SUCCESS);
            default:
                exit(EXIT_FAILURE);
        }

    if(Vflag)
        fprintf(stderr,"%s %s START\n",ctimer(),"tcc-omp");

    /* create a cylinder template and initialize */
    CylinderTmpl *T = (CylinderTmpl*) malloc(sizeof(CylinderTmpl));
    char Tname[50];
    switch(Topt) /* select type of cylinder segmentation */
    {
        case 'H':
            initHammerTmpl(T);
            strcpy(Tname,"Hammer"); 
            break;
        case 'R':
            initHBowieRapp(T);
            strcpy(Tname,"Bowie-Rapp"); 
            break;
        case 'N':
            initHBowieNewell(T); 
            strcpy(Tname,"Bowie-Newell"); 
            break;
        case 'T':
            initTmplTest(T); 
            strcpy(Tname,"Test"); 
            break;
        case 'S':
            initTmpl(T); /* default (finest segmentation) */
            strcpy(Tname,"Olgun-Ustun"); 
            break;
	case 'K':
            initTmplnew(T); /* default (finest segmentation) */
            strcpy(Tname,"Olgun-Ustun_new"); 
            break;
        default:
            exit(EXIT_FAILURE);
    }
    if(Vflag)
        fprintf(stderr,"%s Template type is %s\n",ctimer(),Tname);

    if(Gflag)
    {
        void *API;  /* The API control structure */
        /* Initialize the GMT session for zero pad[0] = {0,0} */
        API = GMT_Create_Session ("tcc",0U,0,NULL);
        /* GMT grid (NetCDF) structure to import outer grid */
        struct GMT_GRID *G = read_GMT_Grid2D(API,gfn,Vflag); 

        /* Tabulated geographical file pointer (longitude, latitude, height) */
        FILE *sites = (FILE*) read_data(argc-optind,argv[optind],Vflag);
        /* Observation site pointer */
        GeodPos *P=malloc(sizeof(GeodPos));
        P->lon=0.0;P->lat=0.0;P->H=0.0;

        if(Iflag)
        {
            /* GMT grid (NetCDF) structure to hold inner grid */
            struct GMT_GRID *I = read_GMT_Grid2D(API,ifn,Vflag);

            if(rflag)
            {
                /* GMT grid (NetCDF) structure to generate innermost 
                 * grid */
                /* Resampling (ires) inner grid wrt computation point,
                 * its zone radius  */
                struct GMT_GRID *C = NULL;
                while(fgets(line,256,sites)!=NULL)
                {
                    if(sscanf(line,"%lf%lf%lf\n",&P->lon,&P->lat,&P->H)==3)
                    {
                        C = GMT_Grid_Resamp(API,I,P,T->zoneR[irZ],ires,Vflag);
                        /* Fill template from GMT_GRID (innermost) */
                        Grid2Tmpl(C,P,0,irZ,T,Vflag);
                        if(Vflag)
                            fprintf(stderr,"%s Innermost zone evaluated\n",ctimer());

                        Grid2Tmpl(I,P,irZ,iZ,T,Vflag);
                        if(Vflag)
                            fprintf(stderr,"%s Inner zone evaluated\n",ctimer());

                        Grid2Tmpl(G,P,iZ,T->zN-1,T,Vflag);
                        if(Vflag)
                            fprintf(stderr,"%s Outer zone evaluated\n",ctimer());

                        tcfromTmpl(T,P,Vflag);
                       
		    	//printf("%.15g %.15g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n",
                        //        P->lon,P->lat,P->H,
                        //        T->ba,T->bb,T->tc,T->bt,T->mH,T->rms);
			printf("%.15g %.15g %10.3f %10.3f %8.3f %8.3f %10.3f\n",
                                P->lon,P->lat,P->H,
                                T->ba,T->bb,T->tc,T->bt);
                       
                        /* delete template variables accumlating in */
                        delTmpl(T);
                        ires[strlen(ires)]='s';
                        ires[strlen(ires)]='\0';
                    } 
                }
            }
            else
            {
                while(fgets(line,256,sites)!=NULL)
                {
                    if(sscanf(line,"%lf%lf%lf\n",&P->lon,&P->lat,&P->H)==3)
                    {
                        Grid2Tmpl(I,P,0,iZ,T,Vflag);
                        if(Vflag)
                            fprintf(stderr,"%s Inner zone evaluated\n",ctimer());

                        Grid2Tmpl(G,P,iZ,T->zN-1,T,Vflag);
                        if(Vflag)
                            fprintf(stderr,"%s Outer zone evaluated\n",ctimer());
                        /* Fill template from grid2D table (for inner)*/
                        tcfromTmpl(T,P,Vflag);
                        printf("%.15g %.15g %10.3f %10.3f %8.3f %8.3f %10.3f %8.3f %8.3f\n",
                                P->lon,P->lat,P->H,
                                T->ba,T->bb,T->tc,T->bt,T->mH,T->rms);
                        /* delete template variables accumlating in */
                        delTmpl(T);
                    } 
                }
            }
        }
        else
        {
            while(fgets(line,256,sites)!=NULL)
            {
                if(sscanf(line,"%lf%lf%lf\n",&P->lon,&P->lat,&P->H)==3)
                {
                    /* Fill template from grid2D table (for inner)*/
                    Grid2Tmpl(G,P,0,T->zN-1,T,Vflag);
                    if(Vflag)
                        fprintf(stderr,"%s Outer zone evaluated\n",ctimer());
                    tcfromTmpl(T,P,Vflag);
                    printf("%.15g %.15g %.7g %.7g %.7g %.7g %.7g %.7g %.7g\n",
                            P->lon,P->lat,P->H,
                            T->ba,T->bb,T->tc,T->bt,T->mH,T->rms);
                    /* delete template variables accumlating in */
                    delTmpl(T);
                } 
            }
        }
        free(P);
        GMT_Destroy_Session(API);
    }
    else
    {
        fprintf(stderr,"-G<grid_name> is required\n");
        exit(EXIT_FAILURE);
    }
    Destroy_Template(T);
    return EXIT_SUCCESS;
}

FILE* read_data(int check,char *name,int Vflag)
{
    FILE* df;
    switch(check)
    {
        case 0:
            df=stdin;
            if(Vflag)
                fprintf(stderr,"%s Input coordinates will be read from stdin\n",ctimer());
            break;
        case 1:
            if((df=fopen(name,"r"))==NULL)
            {
                fprintf(stderr,"%s not found!!!\n", name);
                exit(EXIT_FAILURE);
            }
            if(Vflag)
                fprintf(stderr,"%s Input coordinates will be read from external file, %s\n",ctimer(),name);
            break;
        default:
            exit(EXIT_FAILURE);
    }
    return df;
}
