/*  tc-cylinder : high-accurate terrain corrections 
 *  corrected: 18.02.2023 
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
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include"tcc.h"
int readgrd(Grid2D *g,char *name,double w,double e,double s,double n,double dx,double dy)
{
    /* Open grid file in binary for reading */
    FILE *gin=fopen(name,"rb");
        if(!gin){
        fprintf(stderr,"File I/O error: %s\n",name);
        exit(EXIT_FAILURE);
    }

    /*strcpy(g->name,name);*/
    g->xTL=w;
    g->yTL=n;
    g->xinc=dx;
    g->yinc=dy;
    g->ny=F2RND((n-s)/dy+1);
    g->nx=F2RND((e-w)/dx+1);
    g->H=malloc(g->nx*g->ny*sizeof(double));
    fread(g->H,sizeof(double),g->nx*g->ny,gin);
    fclose(gin);

    return EXIT_SUCCESS;
}

void InfoGrid(const Grid2D *g)
{
    float min=.0,max=.0;
    minmaxGrid(g,&min,&max);
    printf("Grid\n");
    printf(" size: %d %d\n",g->nx,g->ny);
    printf(" TL coordinates: %.12g %.12g\n",g->xTL,g->yTL);
    printf(" increments: %.12g %.12g\n",g->xinc,g->yinc);
    printf(" z bounds: %.12g %.12g\n",min,max);
}

int minmaxGrid(const Grid2D *g,float *min,float *max)
{
    float *pH=g->H;
    float *eH=g->H+g->nx*g->ny;
    *min=*max=g->H[0];
    for(;pH<eH;pH++)
    {
        if(*pH<*min) *min=*pH;
        if(*pH>*max) *max=*pH;
    }
    return EXIT_SUCCESS;
}

void initHammerTmpl(CylinderTmpl *T)
{
    int i=0;
    T->zN=19;
    T->tc=.0;
    T->rms=.0;
    T->ba=.0;
    T->bb=.0;
    T->bt=.0;

    T->zoneR   = malloc(T->zN*sizeof(double));
    T->sectorN = malloc((T->zN-1)*sizeof(int));

    /* Template sturucte and its compartment assignments */
    T->zoneR[ 0]=     2; T->sectorN[ 0]= 8;
    T->zoneR[ 1]=   100; T->sectorN[ 1]= 8;
    T->zoneR[ 2]=   200; T->sectorN[ 2]= 8;
    T->zoneR[ 3]=   300; T->sectorN[ 3]= 8;
    T->zoneR[ 4]=   400; T->sectorN[ 4]= 8;
    T->zoneR[ 5]=   500; T->sectorN[ 5]= 8;
    T->zoneR[ 6]=  1000; T->sectorN[ 6]= 8;
    T->zoneR[ 7]=  1500; T->sectorN[ 7]= 8;
    T->zoneR[ 8]=  2000; T->sectorN[ 8]= 8;
    T->zoneR[ 9]=  3000; T->sectorN[ 9]= 8;
    T->zoneR[10]=  4000; T->sectorN[10]= 8;
    T->zoneR[11]=  6000; T->sectorN[11]= 8;
    T->zoneR[12]=  8000; T->sectorN[12]= 8;
    T->zoneR[13]= 11000; T->sectorN[13]= 8;
    T->zoneR[14]= 15000; T->sectorN[14]= 8;
    T->zoneR[15]= 20000; T->sectorN[15]= 8;
    T->zoneR[16]= 26000; T->sectorN[16]= 8;
    T->zoneR[17]= 33000; T->sectorN[17]= 8;
    T->zoneR[18]= 42000;

    /* find total number of compartments */
    T->cN=0;
    for(i=0;i<T->zN-1;i++)
        T->cN += T->sectorN[i];

    /* allocate memory for compartment variables */
    T->compN=malloc(T->cN*sizeof(int));
    T->mean=malloc(T->cN*sizeof(double));
    for(i=0;i<T->cN;i++)
    {
        T->compN[i]=0;
        T->mean[i]=.0;
    }
}


void initHBowieNewell(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->zN=18;
    Tmpl->tc=.0;
    Tmpl->rms=.0;

    Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));

    Tmpl->zoneR[ 0]=     2; Tmpl->sectorN[ 0]= 4;
    Tmpl->zoneR[ 1]=  16.6; Tmpl->sectorN[ 1]= 6;
    Tmpl->zoneR[ 2]=  53.3; Tmpl->sectorN[ 2]= 6;
    Tmpl->zoneR[ 3]= 170.1; Tmpl->sectorN[ 3]= 8;
    Tmpl->zoneR[ 4]= 390.1; Tmpl->sectorN[ 4]= 8;
    Tmpl->zoneR[ 5]= 894.9; Tmpl->sectorN[ 5]=12;
    Tmpl->zoneR[ 6]=  1530; Tmpl->sectorN[ 6]=12;
    Tmpl->zoneR[ 7]=  2615; Tmpl->sectorN[ 7]=12;
    Tmpl->zoneR[ 8]=  3500; Tmpl->sectorN[ 8]=16;
    Tmpl->zoneR[ 9]=  6653; Tmpl->sectorN[ 9]=16;
    Tmpl->zoneR[10]=  9903; Tmpl->sectorN[10]=16;
    Tmpl->zoneR[11]= 14742; Tmpl->sectorN[11]=16;
    Tmpl->zoneR[12]= 21944; Tmpl->sectorN[12]=20;
    Tmpl->zoneR[13]= 33000; Tmpl->sectorN[13]=20;
    Tmpl->zoneR[14]= 50000; Tmpl->sectorN[14]=20;
    Tmpl->zoneR[15]= 75000; Tmpl->sectorN[15]=20;
    Tmpl->zoneR[16]=110000; Tmpl->sectorN[16]=20;
    Tmpl->zoneR[17]=166740;
    
    Tmpl->cN=0;
    for(i=0;i<Tmpl->zN-1;i++)
        Tmpl->cN += Tmpl->sectorN[i];

    Tmpl->compN=malloc(Tmpl->cN*sizeof(int));
    Tmpl->mean=malloc(Tmpl->cN*sizeof(double));
    for(i=0;i<Tmpl->cN;i++)
    {
        Tmpl->compN[i]=0;
        Tmpl->mean[i]=.0;
    }
}


void initHBowieRapp(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->zN=15;
    Tmpl->tc=.0;
    Tmpl->rms=.0;

    Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));

    Tmpl->zoneR[ 0]=     2; Tmpl->sectorN[ 0]= 4;
    Tmpl->zoneR[ 1]=    68; Tmpl->sectorN[ 1]= 4;
    Tmpl->zoneR[ 2]=   230; Tmpl->sectorN[ 2]= 6;
    Tmpl->zoneR[ 3]=   590; Tmpl->sectorN[ 3]= 8;
    Tmpl->zoneR[ 4]=  1280; Tmpl->sectorN[ 4]=10;
    Tmpl->zoneR[ 5]=  2290; Tmpl->sectorN[ 5]=12;
    Tmpl->zoneR[ 6]=  3500; Tmpl->sectorN[ 6]=16;
    Tmpl->zoneR[ 7]=  5240; Tmpl->sectorN[ 7]=20;
    Tmpl->zoneR[ 8]=  8440; Tmpl->sectorN[ 8]=16;
    Tmpl->zoneR[ 9]= 12400; Tmpl->sectorN[ 9]=20;
    Tmpl->zoneR[10]= 18800; Tmpl->sectorN[10]=24;
    Tmpl->zoneR[11]= 28800; Tmpl->sectorN[11]=14;
    Tmpl->zoneR[12]= 58800; Tmpl->sectorN[12]=16;
    Tmpl->zoneR[13]= 99000; Tmpl->sectorN[13]=28;
    Tmpl->zoneR[14]=166700;
    

    Tmpl->cN=0;
    for(i=0;i<Tmpl->zN-1;i++)
        Tmpl->cN += Tmpl->sectorN[i];

    Tmpl->compN=malloc(Tmpl->cN*sizeof(int));
    Tmpl->mean=malloc(Tmpl->cN*sizeof(double));
    for(i=0;i<Tmpl->cN;i++)
    {
        Tmpl->compN[i]=0;
        Tmpl->mean[i]=.0;
    }
}

void initTmpl(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->zN=134;
    Tmpl->tc=.0;
    Tmpl->rms=.0;

    Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));
  
    Tmpl->zoneR[0]=2; Tmpl->sectorN[0]=4;
    Tmpl->zoneR[1]=20; Tmpl->sectorN[1]=6;
    
    for (i=2;i<32;i=i+1)
     {
        Tmpl->zoneR[i]=(i-1)*50; Tmpl->sectorN[i]=6*((i-1)*50)/50;
     }

     for (i=32;i<43;i=i+1)
     {
         Tmpl->zoneR[i]=1500+(i-31)*100; Tmpl->sectorN[i]=6*(1500+(i-31)*100)/100;
     }

     for (i=43;i<56;i=i+1)
     {
         Tmpl->zoneR[i]=2600+(i-42)*200; Tmpl->sectorN[i]=6*(2400+(i-42)*200)/250;
     }

     for (i=56;i<66;i=i+1)
     {
         Tmpl->zoneR[i]=5200+(i-55)*550; Tmpl->sectorN[i]=6*(5000+(i-55)*550)/750;
     }

     for (i=66;i<77;i=i+1)
     {
         Tmpl->zoneR[i]=10700+(i-65)*1000; Tmpl->sectorN[i]=6*(10000+(i-65)*1000)/1000;
     }
   
     for (i=77;i<97;i=i+1)
     {
         Tmpl->zoneR[i]=21700+(i-76)*1500; Tmpl->sectorN[i]=6*(20000+(i-76)*1500)/1500;
     }

    for (i=97;i<125;i=i+1)
    {
        Tmpl->zoneR[i]=51700+(i-96)*2500; Tmpl->sectorN[i]=6*(50000+(i-96)*2500)/2500;
    }

    for (i=125;i<134;i=i+1)
    {
        Tmpl->zoneR[i]=121700+(i-124)*5000; Tmpl->sectorN[i]=6*(120000+(i-124)*5000)/5000;
    }
    
    Tmpl->cN=0;
    for(i=0;i<Tmpl->zN-1;i++)
    Tmpl->cN += Tmpl->sectorN[i];

    Tmpl->compN=malloc(Tmpl->cN*sizeof(int));
    Tmpl->mean=malloc(Tmpl->cN*sizeof(double));
    for(i=0;i<Tmpl->cN;i++)
    {
        Tmpl->compN[i]=0;
        Tmpl->mean[i]=.0;
    }
}

void DestroyTemplate(CylinderTmpl *tmp)
{
    free(tmp->zoneR);
    free(tmp->sectorN);
    free(tmp->compN);
    free(tmp->mean);
}

SubFrame *subGridFrame(const Grid2D *grd,double lat,double lon,double radi)
{
    SubFrame *sgrd = malloc(sizeof(SubFrame));

    sgrd->idx=0;
    int i=0;
    int py=0,px=0; /* grid coordinates of computation point */
    int sy=0,sx=0; /* bounding box size (half) for Hammer template */ 

    double RG = RGauss(lat);
    double delta2 = .0; 
    double Dlamda = .0; 

    px=F2RND((lon-grd->xTL)/grd->xinc); /* nearest pixel from w */
    py=F2RND((grd->yTL-lat)/grd->yinc); /* nearest pixel from n */

    /* Find longitude difference (Dlamda) of east point 
       with Az = pi/2 and S = radius */
    sphdirsol(M_PI/2,M_PI/2-lat,radi/RG,&delta2,&Dlamda);

    sx=(int)(Dlamda/grd->xinc)+1; /* width in pixel along e-w */
    sy=(int)(radi/RG/grd->yinc)+1;/* height in pixel along n-s */

    /* compute TL coordinates, size and increment of subgrid */
    sgrd->nx=2*sx+1;
    sgrd->ny=2*sy+1;
    sgrd->xTL=grd->xTL+(px-sx)*grd->xinc;
    sgrd->yTL=grd->yTL-(py-sy)*grd->yinc;
    sgrd->iTL=(py-sy)*grd->nx+(px-sx); /* orginal grd index of subgrid TL*/
    return sgrd;
}

double RGauss(double lat)
{
    double cl=cos(lat);
    
    return cE/sqrt(1+ee2*cl*cl);
}

int sphdirsol(double alfa, double delta1, double s, double *delta2, double *Dlamda)
{
    double ca2 = cos(alfa/2);
    double sa2 = sin(alfa/2);

    double z1  = ca2*cos((s-delta1)/2);
    double z2  = ca2*sin((s-delta1)/2);
    double n1  = sa2*cos((s+delta1)/2);
    double n2  = sa2*sin((s+delta1)/2);

    *Dlamda = atan(z1/n1)+atan(z2/n2);
    *delta2 = 2*atan(sqrt((z2*z2+n2*n2)/(z1*z1+n1*n1)));
    return EXIT_SUCCESS;
}
int sphinvsol(double Dlamda, double delta1, double delta2,double *s, double *alfa)
{
    double Dm = (delta1-delta2)/2.;
    double Dp = (delta1+delta2)/2.;
    double ca = cos(Dlamda/2);
    double sa = sin(Dlamda/2);
    double z1  = ca*cos(Dm);
    double z2  = ca*sin(Dm);
    double n1  = sa*cos(Dp);
    double n2  = sa*sin(Dp);

    *alfa = atan(z1/n1)-atan(z2/n2);
    if(*alfa<0.) *alfa += 2*M_PI;
    *s = 2*atan(sqrt((z2*z2+n2*n2)/(z1*z1+n1*n1)));
    return EXIT_SUCCESS;
}

void Grid2Tmpl(const Grid2D *G,SubFrame *F,double lat,double lon,int ri,int ro,CylinderTmpl *T)
{
    int i=0,j=0,k=0,l=0;
    int SN=0,tSN=0,tSN0=0;
    float *spH,*epH;
    double x=0,y=0,s=0,a=0;
    double RG = RGauss(lat);
    double delta=.0,Dalfa= .0;
    int iw = 0; /* index for subframe row */

    /* count the number of sectors for inner grid */
    for(i=0;i<ro;i++)
        tSN0+=T->sectorN[i];
    for(i=0;i<F->ny;i++) /* row loop: top to bottom  */
    {
        y=F->yTL-i*G->yinc;              /* latitude of row */
        x=F->xTL;        /* longitude first element of row */
        iw=F->iTL+i*G->nx;
        spH=&G->H[iw];       /* first pointer address of row */
        epH=&G->H[iw+F->nx]; /* last pointer address of row */
        while(spH<epH)            /* loop over row */
        {
            /* compute spherical polar coordinates */
            sphinvsol(x-lon,M_PI/2-lat,M_PI/2-y,&s,&a);
            s*=RG;  /* s for Gauss sphere */
            if(s>T->zoneR[ro]) goto kout;
            /* find zone and sector and place the height into the cell */
            tSN=tSN0-1;
	    for(k=ro-1;k>=ri;k--) /* loop from inner zone to outer zone */
            {
                SN=T->sectorN[k];
                if(s>T->zoneR[k])
                {
                    Dalfa=2*M_PI/SN;
                    for(l=SN-1;l>=0;l--) /* loop over zone */
                    {
                        if(a>=l*Dalfa) /* find sector */
                        {
                            T->compN[tSN]+=1;
			    delta=*spH - T->mean[tSN];
                            T->mean[tSN]+=delta/T->compN[tSN];
                            goto kout;

                        }
                        else
				tSN--;
                    }
                }
                else
                    tSN-=SN;
	     }	
            kout:
            x+=G->xinc; /* latitude increment */
            spH++;
        }
    }
}
    
void tcfromTmpl(CylinderTmpl *T,double hgt,int v)
{
        int i=0;
        int k=0;
        int l=0;
        int SN = 0;    /* Number of compartment*/
        int tSN= 0;    /* Total number of compartment*/
        double a1 =.0; /* inner zone radius    */
        double a2 =.0; /* auter zone radius    */
        double mR =.0; /* Mean radius of zone  */
        double dA =.0; /* compartment angle    */
        double mA =.0; /* mean azimut of compartment */
        double DH =.0; /* height difference between mean compartment and comp. point heights */
        double dg =.0; /* attraction of compartment  */
        double Pn=.0; //
	double dg2=0.0;
        int ka=1; //coefficient of negatif or pozitif tc
        /* attraction of compartments */
        T->ba=-0.1119*hgt;
        T->bb=-1.464139*pow(10,-3)*hgt-3.533047*pow(10,-7)*pow(hgt,2)+1.002709*pow(10,-13)*pow(hgt,3)+3.002407*pow(10,-18)*pow(hgt,4);
            
            for(k=0;k<T->zN-1;k++)
            {
             a1 = T->zoneR[k];
             a2 = T->zoneR[k+1];
             SN = T->sectorN[k];
             mR = (a1+a2)/2;
             dA = 2*M_PI/SN;
                 
             for(l=0;l<SN;l++)
             {
		     if(T->compN[i]!=0)
		     {
			     tSN++;
                     	     DH = T->mean[i]-hgt;
                             T->rms+=DH*DH;
                             dg2=0.0;
			     if(T->mean[i]>=0)
			     {
				     Pn = sqrt(pow((6371000.0+hgt),2)+mR*mR)-6371000.0;
				     if(T->mean[i] > Pn)
					     DH=-2*Pn+hgt+T->mean[i];
				     else if (Pn > T->mean[i] > hgt)
					     ka=-1;
				     else
					     ka=1;
				     dg=ka*fabs(0.1119*((a2-a1)-sqrt((a2*a2)+DH*DH)+sqrt((a1*a1)+DH*DH))/SN);
			     }
			     else
			     {
				     dg2=0.1119*((a2-a1)-sqrt((a2*a2)+hgt*hgt)+sqrt((a1*a1)+hgt*hgt))/SN;
				     dg=fabs(0.0687*(sqrt(a1*a1+DH*DH)-sqrt(a2*a2+DH*DH)-sqrt(a1*a1+hgt*hgt)+sqrt(a2*a2+hgt*hgt))/SN);
			     }
			     T->tc+=dg+dg2;
		     }
		     i++;
             }
         }
         T->rms=sqrt(T->rms/tSN);
         T->bt=T->ba+T->bb+T->tc;
}

void delTmpl(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->tc=.0;
    Tmpl->rms=.0;

    for(i=0;i<Tmpl->cN;i++)
    {
        Tmpl->compN[i]=0;
        Tmpl->mean[i]=.0;
    }
}

void help()
{
	fprintf(stderr,"tc-cylinder is a C language based open source program which calculate \n");
	fprintf(stderr,"high-accurate terrain corrections. It also provides complete Bouguer reductions \n");
	fprintf(stderr,"with spherical form.\n\n");
	fprintf(stderr,"The algorithm is compatible with DEM input data in both the NetCDF and GeoTIFF formats.\n");
	fprintf(stderr,"The user can specify arbitrary DEM resolutions, the inner region radius,\n");
	fprintf(stderr,"the interpolation region and segmentation type properties. However,\n");
	fprintf(stderr,"the usage pattern highlighted in Example 1 is recommended thanks to \n");
	fprintf(stderr,"its high speed, efficiency, and accuracy.\n\n");
	fprintf(stderr,"GMT 6.0.0 (General Mapping Tools) must be installed due to the use of its some tools.\n");
	fprintf(stderr,"GMT's library should be used when compiling.\n\n");
	fprintf(stderr,"an example of compiling:\n");
	fprintf(stderr,"gcc -o tc-cylinder tc-cylinder.c tcc.c -lm -I/usr/include/gmt/ -L/pub/opt/gmt/lib/ -lgmt\n\n");
	fprintf(stderr,"REQUIRED INPUTS AND ARGUMENTS \n\n");
	fprintf(stderr,"calculation points\n");
        fprintf(stderr,"\t tc-cylinder reads calculation points (longtitude, latitude, height) from standard input.\n\n");
	fprintf(stderr,"-G \n");
        fprintf(stderr,"\t This argument refer to outer or whole computation area DEM.\n" );
        fprintf(stderr,"\t If -I argument (inner grid) also used -Ggridfile will be used for only \n");
        fprintf(stderr,"\t outer zones by tc-cylinder. Any resolution can be chosen as a preference.\n\n");
	fprintf(stderr,"OPTIONAL INPUTS AND ARGUMENTS\n\n");
	fprintf(stderr,"-I \n");
	fprintf(stderr,"\t Specifies the name of the Digital Elevation Model (DEM) used for the interior zone.\n");
	fprintf(stderr,"\t The format should be 2D NetCDF. Any resolution can be chosen as a preference,\n");
	fprintf(stderr,"\t however it should be more precise than the model used with the -G variable.\n");
	fprintf(stderr,"\t The boundaries of this file must be the same as the model corresponding to \n");
	fprintf(stderr,"\t the outer region specified by -G.\n\n");
	fprintf(stderr,"-Z\n");
	fprintf(stderr,"\t The user defines the zone within which the distinction between inner and outer zones\n");
	fprintf(stderr,"\t will be formed using this argument. The default value is 56, which corresponds to \n");
	fprintf(stderr,"\t a distance of 5000 meters from the computation point.\n\n");
	fprintf(stderr,"-r\n");
	fprintf(stderr,"\t The resolution at which the interpolation region should be applied can be determined.\n");
	fprintf(stderr,"\t s is the unit for seconds, m is the unit for minutes, and degrees should be written\n");
	fprintf(stderr,"\t without a unit. 0.5s is the default value.\n\n");
	fprintf(stderr,"-p\n");
	fprintf(stderr,"\t It is possible to specify the zone up to which the interpolation will be applied.\n");
	fprintf(stderr,"\t By default, the value is 3.\n");
	fprintf(stderr,"-T \n");
	fprintf(stderr,"\t This argument specifies the format of the partitions. By default, tc-cylinder uses\n");
	fprintf(stderr,"\t its own detailed partitioning template with S. It is also used with R, H, and N for\n");
	fprintf(stderr,"\t the traditional stencils Hayford Bowie, Hammer, and Niethammer, respectively.\n\n");
	fprintf(stderr,"-h \n");
	fprintf(stderr,"\t help\n\n");
	fprintf(stderr,"OUTPUT FILE \n\n");
	fprintf(stderr,"The following are the outputs and their units in column order:\n\n");
	fprintf(stderr,"Longtitude[degree], Latitude[degree], Height[m], Bouguer Plate Correction[mGal],\n");
	fprintf(stderr,"Spherical Correction[mGal], Terrain Correction[mGal], Complete Bouguer Correction[mGal]\n\n");
	fprintf(stderr,"EXAMPLES AND TEST DATA USAGE\n\n");
	fprintf(stderr,"Example 1\n\n");
	fprintf(stderr,"If 1s resolution is desired in the inner region and 15s resolution in the outer \n");
	fprintf(stderr,"when calculating the correction values for the calculation points on a profile \n");
	fprintf(stderr,"in the Everest region, the following form can be used. Up to zone 56 is considered \n");
	fprintf(stderr,"the inner zone. The interpolation region is preferred up to the third 3th zone at \n");
	fprintf(stderr,"a resolution of 0.5 second. The partitioning template is based on \n");
	fprintf(stderr,"the tc-cylinder template.\n\n");
	fprintf(stderr,"./tc-cylinder pois.lfH -GDEM15s.grd -IDEM1s.grd -TS -Z56 -p3 -r0.5s\n\n");
	fprintf(stderr,"Example 2 \n\n");
	fprintf(stderr,"In contrast to Example 1, if a resolution of 1 second is preferred for the entire region;\n\n");
	fprintf(stderr,"./tc-cylinder pois.lfH -GDEM1s.grd -TS -Z56 -p3 -r0.5s\n\n");
	fprintf(stderr,"Example 3 \n\n");
	fprintf(stderr,"If 3 arc second resolution is desired in the inner region and 15 arc second resolution\n");
	fprintf(stderr,"in the outer when calculating the correction values for the calculation points\n");
	fprintf(stderr,"on a profile in the Everest region, the following form can be used.\n");
	fprintf(stderr,"Up to zone 70 is considered the inner zone. The interpolation region is preferred\n");
	fprintf(stderr,"up to the 10th zone at a resolution of 1 arc second.\n");
	fprintf(stderr,"The partitioning template is based on the tc-cylinder template.\n\n");
	fprintf(stderr,"./tc_cylinder testpois.lfH -GtestDEM15s.grd -ItestDEM3s.grd -TS -Z70 -p10 -r1s \n");
	fprintf(stderr,"Test Data Usage\n\n");
	fprintf(stderr,"If 15 arc second resolution DEM is desired for whole area\n");
	fprintf(stderr,"includes a profile in the Everest region, the following form can be used for test.\n");
	fprintf(stderr,"The interpolation region is preferred\n");
	fprintf(stderr,"up to the 20th zone at a resolution of 0.5 arc second.\n");
	fprintf(stderr,"The partitioning template is based on the tc-cylinder template.\n\n");
	fprintf(stderr,"./tc_cylinder testpois.lfH -GtestDEM15s.grd -TS -p20 -r0.5s > testoutput.txt \n");
}
