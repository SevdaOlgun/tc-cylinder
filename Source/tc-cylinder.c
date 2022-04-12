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
 *  Created: March 2022 	v1.0 Program Created
 */

#include<unistd.h>
#include<omp.h>
#include "gmt.h"
#include "tcc.h"
#include<stdio.h>
#define OPTIONS "G:I:T:Z:r:p:hV"
int main (int argc, char *argv[])
{
/* command line option flags */
    int Iflag = 0;
    int Vflag = 0;
    int rflag = 0;
    double hgt1 = 0;
    /* arguments and variables */
    char Topt='S';  /* default template type S (detailed segmentation) */
    char *gfn;      /* Coarse (outter) grid file name */
    char *ifn;      /* Fine (inner) grid file name */
    char *ires;     /* interpolation resolution of innermost grid */
    int iZ = 56;    /* Inner grid radius index (see templates) */
    int irZ=3;      /* Interpolation zone*/
    double glon=.0,glat=.0,Hmsl=.0;
    char  line[256];
    char  c;
    while((c=getopt(argc,argv,OPTIONS))!=-1)
            switch(c)
        {
            case 'G':
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
    void *API;                  /* The API control structure */
    struct GMT_GRID *G = NULL;  /* Structure to hold outer grid */
    struct GMT_GRID *I = NULL;  /* Structure to hold inner grid */
    struct GMT_GRID *C = NULL;  /* Structure to hold inner grid */

    char args[128] = {""};           /* str for module cmd args */
    char input[GMT_VF_LEN] = {""};    /* virtual input filename */
    char output[GMT_VF_LEN] = {""};   /* virtual output filename */

    int i=0,row=0,col=0,node=0;

    /* Initialize the GMT session for pad[2] = {2,2} */
    API = GMT_Create_Session ("tcc",0U,0,NULL);
    /* Read in our data table to memory */
    G = GMT_Read_Data(API,
            GMT_IS_GRID,
            GMT_IS_FILE,
            GMT_IS_SURFACE,
            GMT_READ_NORMAL, 
            NULL,gfn,NULL);
    /*... create a Grid2D from G */
    Grid2D *topoG = malloc(sizeof(Grid2D));
    float *pH = NULL;
    topoG->xTL  = G->header->wesn[0]/R2D;
    topoG->yTL  = G->header->wesn[3]/R2D;
    topoG->nx   = G->header->n_columns;
    topoG->ny   = G->header->n_rows;
    topoG->xinc = G->header->inc[0]/R2D;
    topoG->yinc = G->header->inc[1]/R2D;
    topoG->H=malloc(topoG->nx*topoG->ny*sizeof(float));
    topoG->H=G->data;
    //InfoGrid(topoG);

    /* create a cylinder template and initialize */
    CylinderTmpl *T = malloc(sizeof(CylinderTmpl));
    switch(Topt)
    {
        case 'H':
            initHammerTmpl(T);
            break;
        case 'R':
            initHBowieRapp(T);
            break;
        case 'N':
            initHBowieNewell(T); 
            break;
        case 'S':
            initTmpl(T);
            break;
        default:
            exit(EXIT_FAILURE);
    }

    Grid2D *topoI = malloc(sizeof(Grid2D));
    Grid2D *topoC = malloc(sizeof(Grid2D));

    if(Iflag)
    {
        /*... create a hgtGrid from I */
        I = GMT_Read_Data(API,
                GMT_IS_GRID,
                GMT_IS_FILE,
                GMT_IS_SURFACE,
                GMT_READ_NORMAL, 
                NULL,ifn,NULL);
        topoI->xTL  = I->header->wesn[0]/R2D;
        topoI->yTL  = I->header->wesn[3]/R2D;
        topoI->nx   = I->header->n_columns;
        topoI->ny   = I->header->n_rows;
        topoI->xinc = I->header->inc[0]/R2D;
        topoI->yinc = I->header->inc[1]/R2D;
        topoI->H=malloc(topoI->nx*topoI->ny*sizeof(float));
        topoI->H = I->data;
     }

    /* option analysis and file opening ****************/
    FILE *sites;
    switch(argc-optind)
    {
        case 0:
            sites=stdin;
            break;
        case 1:
            if((sites=fopen(argv[optind],"r"))==NULL)
            {
                fprintf(stderr,"%s not found!!!\n", argv[optind]);
                exit(EXIT_FAILURE);
            }
            break;
        default:
            exit(EXIT_FAILURE);
    }

    while(fgets(line,256,sites)!=NULL)
    {
                if(sscanf(line,"%lf%lf%lf\n",&glon,&glat,&Hmsl)==3)
                {
                    glat/=R2D;
                    glon/=R2D;
                    /* extract a subgrid surrounding the template for computation point */
                    /* Innermost template for 1-arc second grid resolution interpolated
                    from original grid file*/
                    if(Iflag)
                    {
                        if(rflag)
                        {
                            SubFrame *irsub = subGridFrame(topoI,glat,glon,T->zoneR[irZ]);
                            SubFrame *Isub = subGridFrame(topoI,glat,glon,T->zoneR[iZ]);
                            GMT_Open_VirtualFile(API,
                                GMT_IS_GRID,
                                GMT_IS_SURFACE,
                                GMT_IN,
                                I,input);
                                /* Virtual file to hold the resulting grid in UTM*/
                            GMT_Open_VirtualFile(API,
                                GMT_IS_GRID,
                                GMT_IS_SURFACE,
                                GMT_OUT,
                                NULL,output);
                                /* Prepare the module arguments */
                            sprintf(args,"-R%.14g/%.14g/%.14g/%.14g -I%s %s -G%s",
                                R2D*irsub->xTL,
                                R2D*(irsub->xTL+topoI->xinc*(irsub->nx-1)),
                                R2D*(irsub->yTL-topoI->yinc*(irsub->ny-1)),
                                R2D*irsub->yTL,
                                ires,
                                input,output);
                                /* Call the grdproject module */
                            GMT_Call_Module (API,"grdsample",GMT_MODULE_CMD,args);
                            /* Obtain the grid from the virtual file */
                            C = GMT_Read_VirtualFile (API,output);
                            /* Close the virtual files */
                            GMT_Close_VirtualFile(API,input);
                            GMT_Close_VirtualFile(API,output);
                            topoC->xTL  = C->header->wesn[0]/R2D;
                            topoC->yTL  = C->header->wesn[3]/R2D;
                            topoC->nx   = C->header->mx;
                            topoC->ny   = C->header->my;
                            topoC->xinc = C->header->inc[0]/R2D;
                            topoC->yinc = C->header->inc[1]/R2D;
                            topoC->H=malloc(topoC->nx*topoC->ny*sizeof(float));
                            topoC->H = C->data;
                            irsub->iTL=C->header->pad[3]*topoC->nx+C->header->pad[0];
                            irsub->nx=C->header->n_columns;
                            irsub->ny=C->header->n_rows;
                            irsub->xTL=topoC->xTL;
                            irsub->yTL=topoC->yTL;
                            Grid2Tmpl(topoC,irsub,glat,glon,0,irZ,T);
                            Grid2Tmpl(topoI,Isub,glat,glon,irZ,iZ,T);
                        }
                        else
                        {
                            SubFrame *Isub = subGridFrame(topoI,glat,glon,T->zoneR[iZ]);
                            Grid2Tmpl(topoI,Isub,glat,glon,0,iZ,T);
                        }
                        SubFrame *Gsub = subGridFrame(topoG,glat,glon,T->zoneR[T->zN-1]);
                        Grid2Tmpl(topoG,Gsub,glat,glon,iZ,T->zN-1,T);
                    }
                    else
                    {
                        if(rflag)
                        {
                                SubFrame *irsub = subGridFrame(topoG,glat,glon,T->zoneR[irZ]);
                                SubFrame *Gsub = subGridFrame(topoG,glat,glon,T->zoneR[T->zN-1]);
                                GMT_Open_VirtualFile(API,
                                GMT_IS_GRID,
                                GMT_IS_SURFACE,
                                GMT_IN,
                                G,input);
                                GMT_Open_VirtualFile(API,
                                    GMT_IS_GRID,
                                    GMT_IS_SURFACE,
                                    GMT_OUT,
                                    NULL,output);
                                    /* Prepare the module arguments */
                                    sprintf(args,"-R%.14g/%.14g/%.14g/%.14g -I%s %s -G%s",
                                    R2D*irsub->xTL,
                                    R2D*(irsub->xTL+topoG->xinc*(irsub->nx-1)),
                                    R2D*(irsub->yTL-topoG->yinc*(irsub->ny-1)),
                                    R2D*irsub->yTL,
                                    ires,
                                    input,output);
                                    /* Call the grdproject module */
                                GMT_Call_Module (API,"grdsample",GMT_MODULE_CMD,args);
                                /* Obtain the grid from the virtual file */
                                C = GMT_Read_VirtualFile (API,output);
                                /* Close the virtual files */
                                GMT_Close_VirtualFile(API,input);
                                GMT_Close_VirtualFile(API,output);
                                topoC->xTL  = C->header->wesn[0]/R2D;
                                topoC->yTL  = C->header->wesn[3]/R2D;
                                topoC->nx   = C->header->mx;
                                topoC->ny   = C->header->my;
                                topoC->xinc = C->header->inc[0]/R2D;
                                topoC->yinc = C->header->inc[1]/R2D;
                                topoC->H=malloc(topoC->nx*topoC->ny*sizeof(float));
                                topoC->H = C->data;
                                irsub->iTL=C->header->pad[3]*topoC->nx+C->header->pad[0];
                                irsub->nx=C->header->n_columns;
                                irsub->ny=C->header->n_rows;
                                irsub->xTL=topoC->xTL;
                                irsub->yTL=topoC->yTL;
                                Grid2Tmpl(topoC,irsub,glat,glon,0,irZ,T);
                                Grid2Tmpl(topoG,Gsub,glat,glon,irZ,T->zN-1,T);
                        }
                        else
                        {
                        SubFrame *Gsub = subGridFrame(topoG,glat,glon,T->zoneR[T->zN-1]);
                        Grid2Tmpl(topoG,Gsub,glat,glon,0,T->zN-1,T);
                        free(Gsub);
                        }
                    }
                    /* calculate terrain correction for computation point */
                    tcfromTmpl(T,Hmsl,Vflag);
                    printf("%.12g %.12g %.5g %.8g %.8g %.8g %.8g\n", glon*R2D,glat*R2D,Hmsl,T->ba,T->bb,T->tc,T->bt);
                    /* freeing subgrid memory and deleting template */
                delTmpl(T);
                }
         }
                  free(topoI);
                  free(topoG);
                  free(topoC);

                  GMT_Destroy_Session(API);
                  return EXIT_SUCCESS;
}
