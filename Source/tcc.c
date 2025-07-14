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
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include"tcc.h"
#include"gmt.h"
#include<omp.h>

double dAttrofMass(double DH,double a1,double a2,int SN,double MG)
{
    double n1 = DH/a1;
    double n2 = DH/a2;
    return MG*(a2*(1.0-sqrt(1.0+n2*n2))-a1*(1.0-sqrt(1.0+n1*n1)))/SN;
}

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
    //g->H=malloc(g->nx*g->ny*sizeof(float));
    posix_memalign((void **)&g->H,64,g->nx*g->ny*sizeof(float));
    fread(g->H,sizeof(float),g->nx*g->ny,gin);
    fclose(gin);

    return EXIT_SUCCESS;
}

void InfoGrid(const Grid2D *g)
{
    float min=.0,max=.0;
    minmaxGrid(g,&min,&max);
    printf("Grid\n");
    printf(" size: %d %d\n",g->nx,g->ny);
    printf(" TL coordinates: %.15g %.15g\n",g->xTL*R2D,g->yTL*R2D);
    printf(" increments: %.15g %.15g\n",g->xinc*R2D,g->yinc*R2D);
    printf(" z bounds: %.15g %.15g\n",min,max);
}

int minmaxGrid(const Grid2D *g,float *min,float *max)
{
    int i=0;
    *min=*max=g->H[0];
    for(i=0;i<g->nx*g->ny;i++)
    {
        if(g->H[i]<*min) *min = g->H[i];
        if(g->H[i]>*max) *max = g->H[i];
    }
    return EXIT_SUCCESS;
}

void initHammerTmpl(CylinderTmpl *T)
{
    T->zN=19;
    T->tc=.0;
    T->rms=.0;
    T->mH =.0;
    T->ba=.0;
    T->bb=.0;
    T->bt=.0;

    //T->zoneR   = malloc(T->zN*sizeof(double));
    //T->sectorN = malloc((T->zN-1)*sizeof(int));
    posix_memalign((void **)&T->zoneR,64,T->zN*sizeof(double));
    posix_memalign((void **)&T->sectorN,64,(T->zN-1)*sizeof(int));
    posix_memalign((void **)&T->zmR,64,(T->zN-1)*sizeof(double));
    posix_memalign((void **)&T->zDR,64,(T->zN-1)*sizeof(double));

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

    initComp(T);
}


void initHBowieNewell(CylinderTmpl *Tmpl)
{
    Tmpl->zN=18;
    Tmpl->tc=.0;
    Tmpl->rms=.0;
    Tmpl->mH =.0;

    //Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    //Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zoneR,64,Tmpl->zN*sizeof(double));
    posix_memalign((void **)&Tmpl->sectorN,64,(Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zmR,64,(Tmpl->zN-1)*sizeof(double));
    posix_memalign((void **)&Tmpl->zDR,64,(Tmpl->zN-1)*sizeof(double));

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

    Tmpl->zoneR[Tmpl->zN-1] = 1.5/R2D*R0;

    initComp(Tmpl);
}


void initHBowieRapp(CylinderTmpl *Tmpl)
{
    Tmpl->zN=15;
    Tmpl->tc=.0;
    Tmpl->rms=.0;
    Tmpl->mH =.0;

    //Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    //Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zoneR,64,Tmpl->zN*sizeof(double));
    posix_memalign((void **)&Tmpl->sectorN,64,(Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zmR,64,(Tmpl->zN-1)*sizeof(double));
    posix_memalign((void **)&Tmpl->zDR,64,(Tmpl->zN-1)*sizeof(double));

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

    Tmpl->zoneR[Tmpl->zN-1] = 1.5/R2D*R0;
    initComp(Tmpl);
}

void initTmpl(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->zN=134;
    Tmpl->tc=.0;
    Tmpl->rms=.0;
    Tmpl->mH =.0;

    posix_memalign((void **)&Tmpl->zoneR,64,Tmpl->zN*sizeof(double));
    posix_memalign((void **)&Tmpl->sectorN,64,(Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zmR,64,(Tmpl->zN-1)*sizeof(double));
    posix_memalign((void **)&Tmpl->zDR,64,(Tmpl->zN-1)*sizeof(double));

    Tmpl->zoneR[0]=2; Tmpl->sectorN[0]=4;
    Tmpl->zoneR[1]=20; Tmpl->sectorN[1]=6;

    for (i=2;i<32;i++)
    {
        Tmpl->zoneR[i]=(i-1)*50; Tmpl->sectorN[i]=6*((i-1)*50)/50;
    }

    for (i=32;i<43;i++)
    {
        Tmpl->zoneR[i]=1500+(i-31)*100; Tmpl->sectorN[i]=6*(1500+(i-31)*100)/100;
    }

    for (i=43;i<56;i++)
    {
        Tmpl->zoneR[i]=2600+(i-42)*200; Tmpl->sectorN[i]=6*(2400+(i-42)*200)/300;
    }

    for (i=56;i<66;i++)
    {
        Tmpl->zoneR[i]=5200+(i-55)*550; Tmpl->sectorN[i]=6*(5000+(i-55)*550)/750;
    }

    for (i=66;i<77;i++)
    {
        Tmpl->zoneR[i]=10700+(i-65)*1000; Tmpl->sectorN[i]=6*(10000+(i-65)*1000)/1000;
    }

    for (i=77;i<97;i++)
    {
        Tmpl->zoneR[i]=21700+(i-76)*1500; Tmpl->sectorN[i]=6*(20000+(i-76)*1500)/1500;
    }

    for (i=97;i<125;i++)
    {
        Tmpl->zoneR[i]=51700+(i-96)*2500; Tmpl->sectorN[i]=6*(50000+(i-96)*2500)/2500;
    }

    for (i=125;i<134;i++)
    {
        Tmpl->zoneR[i]=121700+(i-124)*5000; Tmpl->sectorN[i]=6*(120000+(i-124)*5000)/5000;
    }
    Tmpl->zoneR[Tmpl->zN-1] = 1.5/R2D*R0;

    initComp(Tmpl);
}


void initTmplnew(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->zN=309;
    Tmpl->tc=.0;
    Tmpl->rms=.0;
    Tmpl->mH =.0;

    //Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    //Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zoneR,64,Tmpl->zN*sizeof(double));
    posix_memalign((void **)&Tmpl->sectorN,64,(Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zmR,64,(Tmpl->zN-1)*sizeof(double));
    posix_memalign((void **)&Tmpl->zDR,64,(Tmpl->zN-1)*sizeof(double));

    Tmpl->zoneR[0]=2; Tmpl->sectorN[0]=2;
    Tmpl->zoneR[1]=10; Tmpl->sectorN[1]=4;

    for (i=2;i<101;i++)
    {
        Tmpl->zoneR[i]=(i-1)*25; Tmpl->sectorN[i]=6*((i-1)*25)/30;
    }

    for (i=101;i<121;i++)
    {
        Tmpl->zoneR[i]=2500+(i-100)*50; Tmpl->sectorN[i]=6*(2500+(i-100)*50)/50;
    }

    for (i=121;i<146;i++)
    {
        Tmpl->zoneR[i]=5000+(i-120)*200; Tmpl->sectorN[i]=6*(5000+(i-120)*200)/200;
    }

    for (i=146;i<226;i++)
    {
        Tmpl->zoneR[i]=10000+(i-145)*500; Tmpl->sectorN[i]=6*(10000+(i-145)*500)/600;
    }

    for (i=226;i<276;i++)
    {
        Tmpl->zoneR[i]=50000+(i-225)*1000; Tmpl->sectorN[i]=6*(50000+(i-225)*1000)/1000;
    }

    for (i=276;i<309;i++)
    {
        Tmpl->zoneR[i]=100000+(i-275)*2000; Tmpl->sectorN[i]=6*(100000+(i-275)*2000)/2000;
    }
    Tmpl->zoneR[Tmpl->zN-1] = 1.5/R2D*R0;

    initComp(Tmpl);
}

void initTmplTest(CylinderTmpl *Tmpl)
{
    int i=0;
    Tmpl->zN=18;
    Tmpl->tc=.0;
    Tmpl->rms=.0;
    Tmpl->mH =.0;

    //Tmpl->zoneR   = malloc(Tmpl->zN*sizeof(double));
    //Tmpl->sectorN = malloc((Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zoneR,64,Tmpl->zN*sizeof(double));
    posix_memalign((void **)&Tmpl->sectorN,64,(Tmpl->zN-1)*sizeof(int));
    posix_memalign((void **)&Tmpl->zmR,64,(Tmpl->zN-1)*sizeof(double));
    posix_memalign((void **)&Tmpl->zDR,64,(Tmpl->zN-1)*sizeof(double));


    for(i=0;i<Tmpl->zN-1;i++)
    {
        Tmpl->zoneR[i]  = i*1.0e+4;
        Tmpl->sectorN[i] = 12;
        Tmpl->cN += Tmpl->sectorN[i];
    }
    Tmpl->zoneR[0]=0.2;
    Tmpl->zoneR[Tmpl->zN-1] = 1.5/R2D*R0;

    initComp(Tmpl);
}

void initComp(CylinderTmpl *Tmpl)
{
    int i=0;

    Tmpl->cN=0;
    for(i=0;i<Tmpl->zN;i++)
        Tmpl->zoneR[i] /= R0;
    for(i=0;i<Tmpl->zN-1;i++)
    {
        Tmpl->cN += Tmpl->sectorN[i];
        Tmpl->zmR[i] = (Tmpl->zoneR[i+1]+Tmpl->zoneR[i])/2;
        Tmpl->zDR[i] = Tmpl->zoneR[i+1]-Tmpl->zoneR[i];
        //printf("%6d %4d %.15g %.15g %.15g\n",i,Tmpl->sectorN[i],
        //        Tmpl->zoneR[i]*R0,Tmpl->zDR[i]*R0,Tmpl->zmR[i]*R0);
    }
    //printf("%6d %4d %.15g %.15g %.15g\n",i,Tmpl->sectorN[i],
    //        Tmpl->zoneR[i]*R0,Tmpl->zDR[i]*R0,Tmpl->zmR[i]*R0);

    const int size = Tmpl->cN;
    //Tmpl->compN=malloc(size*sizeof(int));
    //Tmpl->mean=malloc(size*sizeof(double));
    posix_memalign((void **)&Tmpl->compN,64,size*sizeof(int));
    posix_memalign((void **)&Tmpl->mean,64,size*sizeof(double));
    for(i=0;i<size;i++)
    {
        Tmpl->compN[i]=0;
        Tmpl->mean[i]=.0;
    }
}

void Destroy_Template(CylinderTmpl *T)
{
    free(T->zoneR);
    free(T->sectorN);
    free(T->compN);
    free(T->mean);
    free(T->zmR);
    free(T->zDR);
}

SubFrame *subGridFrame(struct GMT_GRID *G,GeodPos *P,double psi)
{
    /*
     * Determines tightly bounded zone frame
     * Verified formulation A. Ustun 29.03.2023
     */
    SubFrame *sgrd = (SubFrame*) malloc(sizeof(SubFrame));

    sgrd->idx=0;
    int py=0,px=0; /* grid coordinates of computation point */
    int sy=0,sx=0; /* bounding box size (half) for template */ 

    double Dlamda = .0; 
    
    /* nearest pixel from w */
    px=F2RND((P->lon-G->header->wesn[0])/G->header->inc[0]); 
    /* nearest pixel from n */
    py=F2RND((G->header->wesn[3]-P->lat)/G->header->inc[1]); 

    /* Find longitude difference (Dlamda) of east point 
       with Az = pi/2 and S = radius */
    Dlamda = atan(tan(psi)/cos(P->lat/R2D));

    /* width in pixel along e-w */
    sx = (int)(Dlamda*R2D/G->header->inc[0])+1;
    /* height in pixel along n-s */
    sy = (int)(psi*R2D/G->header->inc[1])+1;

    /* compute TL coordinates, size and increment of subgrid */
    sgrd->nx = 2*sx+1;
    sgrd->ny = 2*sy+1;
    sgrd->xTL= G->header->wesn[0]+(px-sx)*G->header->inc[0];
    sgrd->yTL= G->header->wesn[3]-(py-sy)*G->header->inc[1];
    /* orginal grd index of subgrid TL*/
    sgrd->iTL= (py-sy)*G->header->n_columns+(px-sx);
    return sgrd;
}

double RGauss(double lat)
{
    double cl=cos(lat);

    return cE/(1.0+ee2*cl*cl);
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
    if(*alfa<0) *alfa += 2*M_PI;
    *s = 2*atan(sqrt((z2*z2+n2*n2)/(z1*z1+n1*n1)));
    return EXIT_SUCCESS;
}

void lltoazeq(const double lon,const double lat, double *azi,double *rad)
{
    double sf = sin(*rad);
    double s0 = sin(lat);
    double cf = cos(*rad);
    double c0 = cos(lat);
    double cl = cos(*azi-lon);
    double sl = sin(*azi-lon);

    *rad = acos(sf*s0+cf*c0*cl);
    *azi = atan2(sl,(c0*sf/cf-s0*cl));
    if(*azi<0.0) *azi += 2.0*M_PI;
}

void Grid2Tmpl(struct GMT_GRID *G,GeodPos *P,const int ri,const int ro,CylinderTmpl *T,int Vflag)
{
    SubFrame *F = subGridFrame(G,P,T->zoneR[ro]);
    int iw = 0; /* index for subframe row */
    int i=0,j=0,k=0,l=0;
    int tSN=0,tSN0=0;
    double x,y,s,a;
    const double rp = RGauss(P->lat/R2D) + P->H;
    double Dalfa= .0;
    const double *zon  = T->zoneR;
    int    *comp = T->compN;
    double *mean = T->mean;
    const float  *pH   = G->data;
    const int    *sct  = T->sectorN;
    const int Fi  = F->iTL;
    const int Fny = F->ny;
    const int Fnx = F->nx;
    const int Gnx = G->header->n_columns;
    const double Fx  = F->xTL+G->header->registration*G->header->inc[0]/2;
    const double Fy  = F->yTL-G->header->registration*G->header->inc[1]/2;
    const double Gix = G->header->inc[0];
    const double Giy = G->header->inc[1];
    const double lat = P->lat;
    const double lon = P->lon;

    int tid=0;  /* thread id */
    int nth=0;  /* # of threads */
    const int chk=10; /* loop iteration over chunk size */

    /* count the number of sectors for inner grid */
    for(i=0;i<ro;i++)
        tSN0 += sct[i];
    if(Vflag)
        fprintf(stderr,"%s [%d] = %.9g -- [%d] = %.9g m; # of compartments in zone (tSN0) = %d\n",
                ctimer(),ri,zon[ri]*rp,ro,zon[ro]*rp,tSN0);

    nth = omp_get_max_threads();
    if(Vflag)
        fprintf(stderr,"%s ***** Parallel processing with %d threads *****\n",ctimer(),nth);

    //omp_set_num_threads(1);
#pragma omp parallel default(none) shared(comp,mean,nth) firstprivate(tSN0,zon,sct,pH,ri,ro,Fny,Fnx,Fi,Giy,Gix,Gnx,Fx,Fy,lon,lat) private(tid,i,j,k,l,tSN,iw,x,y,s,a,Dalfa)
    {
        tid = omp_get_thread_num();
        //if(tid==0)
        //    nth = omp_get_max_threads();
        //printf("%dth thread of %d threads\n",tid,nth);
#pragma omp for schedule (static, 10)
        for(i=0;i<Fny;i++) /* row loop: top to bottom  */
        {
            y  = Fy-i*Giy;  /* latitude of row */
            x  = Fx;        /* longitude first element of row */
            iw = Fi+i*Gnx;
            for(j=iw;j<iw+Fnx;j++,x+=Gix)  /* loop over row */
            {
                /* compute spherical polar coordinates */
                a = x/R2D; s = y/R2D;
                lltoazeq(lon/R2D,lat/R2D,&a,&s);
                //s *= rp;  /* s for Gauss sphere */
                if(s > zon[ro]) continue;
                /* find zone and sector and place the height into the cell */
                tSN = tSN0-1;
                for(k=ro-1;k>=ri;k--) /* loop from inner zone to outer zone */
                {
                    if(s > zon[k])
                    {
                        Dalfa=2*M_PI/sct[k];
                        for(l=sct[k]-1;l>=0;l--) /* loop over zone */
                        {
                            if(a >= l*Dalfa) /* find sector */
                            {
#pragma omp atomic
                                comp[tSN]++;
#pragma omp atomic
                                mean[tSN] += pH[j];
                                k=ri-1;
                                break;
                            }
                            else
                                tSN--;
                        }
                    }
                    else
                        tSN -= sct[k];
                }	
            }
        }
    }
    free(F);
}

void tcfromTmpl(CylinderTmpl *T,GeodPos *P,int vflag)
{
    int i=0;
    int k=0;
    int l=0;
    int tSN=0,SN = 0;    /* Number of compartment*/
    double a1 =.0; /* inner zone radius    */
    double a2 =.0; /* auter zone radius    */
    double mR =.0; /* Mean radius of zone  */
    double DH =.0; /* height difference between mean compartment and comp. point heights */
    double dg =.0; /* attraction of compartment  */
    double dH =.0; //
    double *pM = T->mean;
    int    *pN = T->compN;
    double lat = P->lat/R2D;
    double Hp  = P->H;
    double rp  = RGauss(lat)+Hp;

    /* attraction of compartments */
    T->ba = -TPG*Hp;
    T->bb = Hp*(1.464139e-3+Hp*(-3.533047e-7+Hp*(1.002709e-13+Hp*3.002407e-18)));

    for(k=0;k<T->zN-1;k++)
    {
        a1 = T->zoneR[k]*rp;
        a2 = a1 + T->zDR[k]*rp;
        mR = T->zmR[k]*rp;
        dH = mR*mR/2.0/rp;

        SN = T->sectorN[k];
        for(l=0;l<SN;l++,pM++,pN++,i++)
        {
            if(*pN!=0)
            {
                *pM /= *pN;
                tSN++;
                DH     = *pM-Hp;
                T->mH += *pM;
                T->rms+= *pM**pM;

                if(DH==0.0)
                    dg = 0.0;
		else if(Hp==0.0) /* Obs. point on sea */
		{
			if(DH>0.0) //integrand point in land
			{
				//dg = dAttrofMass(DH,a1,a2,SN,TPG-APG);
				if(DH>dH) //integrand noktasi kuresel levhanin uzerinde hem neg hem poz degerler olmali
					dg = dAttrofMass(DH-dH,a1,a2,SN,TPG-APG)-dAttrofMass(dH,a1,a2,SN,TPG-APG);
				else //bouguer levhasinin �icinde kaliyor neg 
					dg = -dAttrofMass(DH,a1,a2,SN,TPG-APG);
			}
			else //integrand noktasi denizin icinde
				dg = dAttrofMass(*pM,a1,a2,SN,TPG-SPG);
		}
                else if(Hp>0.0) /* Obs. point on land */
		{
			if(DH>0.0) //integrand noktasi P nin yukarısında
			{
				if(DH>dH) //DHnin bir kismi levhanin icinde bir ustunde
					dg = dAttrofMass(DH-dH,a1,a2,SN,TPG-APG)-dAttrofMass(dH,a1,a2,SN,TPG-APG);
				else // DHnin tamami levhanin icinde
					dg = -dAttrofMass(DH,a1,a2,SN,TPG-APG);
			}
			else // integrand noktasi P nin asagisinda
			{
				if(*pM>=0.0) //integrandkarada
					dg = dAttrofMass(DH,a1,a2,SN,TPG-APG);
				else //integrand denizde
				       	dg = dAttrofMass(Hp,a1,a2,SN,TPG-APG)
                                		+dAttrofMass(Hp-*pM,a1,a2,SN,TPG-SPG)-dAttrofMass(Hp,a1,a2,SN,TPG-SPG); // water and gap effect 
			}
		}
                else /* Obs. point in below sea level */
                    fprintf(stderr,"Observation point inside sea!!!!\n");
            }
            else
                DH = dg = 0.0;
            T->tc += dg;
            if(vflag)
                printf("%d %d %d %d %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",i,k,l,SN,*pN,
                        a1,a2,*pM,DH,dH,dg,T->tc);
        }
    }
    //if(tSN!=T->cN)
    //    fprintf(stderr,"There are empty compartments: filled = %d, expected = %d\n",tSN,T->cN);
    T->mH  /= tSN;
    T->rms /= tSN;
    T->rms=sqrt(T->rms - T->mH*T->mH);
    T->bt = T->ba+T->bb+T->tc;
}

void delTmpl(CylinderTmpl *T)
{
    int i=0;
    T->tc=.0;
    T->rms=.0;
    T->mH=.0;
    T->ba=.0;
    T->bb=.0;
    T->bt=.0;

    for(i=0;i<T->cN;i++)
    {
        T->compN[i]=0;
        T->mean[i]=.0;
    }
}

void help()
{
    fprintf(stderr,"tc-cylinder is a C language based open source program which calculate \n");
    fprintf(stderr,"high-accurate terrain corrections. It also provides complete Bouguer reductions \n");
    fprintf(stderr,"with spherical form.\n\n");
    fprintf(stderr,"The user can specify arbitrary DEM resolutions, the inner region radius,\n");
    fprintf(stderr,"the interpolation region and segmentation type properties. However,\n");
    fprintf(stderr,"the usage pattern highlighted in Example 1 is recommended thanks to \n");
    fprintf(stderr,"its high speed, efficiency, and accuracy.\n\n");
    fprintf(stderr,"GMT 6.0.0 (General Mapping Tools) must be installed due to the use of its some tools.\n");
    fprintf(stderr,"GMT's library should be used when compiling.\n\n");
    fprintf(stderr,"Makefile can be use for compiling:\n");
    fprintf(stderr,"OpenMP should be used when compiling\n\n");
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
    fprintf(stderr,"\t with a unit. 0.5s is the default value.\n\n");
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

char *ctimer()
{
    char *c = malloc(80*sizeof(char));

    time_t t;
    struct tm *tmp;
    time(&t);

    tmp = localtime( &t );

    // using strftime to display time
    strftime(c,80,"%Y-%m-%dT%H:%M:%S ::", tmp);

    return c;
}

struct GMT_GRID* read_GMT_Grid2D(void *API,char *fn,int Vflag)
{
    struct GMT_GRID *G =  GMT_Read_Data(API,
            GMT_IS_GRID,
            GMT_IS_FILE,
            GMT_IS_SURFACE,
            GMT_READ_NORMAL, 
            NULL,fn,NULL);
    if(Vflag)
    {
        fprintf(stderr,"----------------------------------------------------------------------------------------------\n");
        fprintf(stderr,"%s Grid file (%s) has been successfuly read\n",ctimer(),fn);
        fprintf(stderr,"%s %s registration\n",ctimer(),(G->header->registration==0)?"Gridline":"Pixel");
        fprintf(stderr,"%s xmin = %.15g xmax = %.15g xinc = %.15g n_xcols = %d\n",ctimer(),G->header->wesn[0],G->header->wesn[1],G->header->inc[0],G->header->n_columns);
        fprintf(stderr,"%s ymin = %.15g ymax = %.15g yinc = %.15g n_yrows = %d\n",ctimer(),G->header->wesn[2],G->header->wesn[3],G->header->inc[1],G->header->n_rows);
        fprintf(stderr,"----------------------------------------------------------------------------------------------\n");
    }
    return G;
}

struct GMT_GRID* GMT_Grid_Resamp(void *API,struct GMT_GRID *G,GeodPos *P,double Ri,char *ires,int Vflag)
{
    SubFrame *F = subGridFrame(G,P,Ri);

    char args[256] = {""};         /* str for module cmd args */
    char inp[GMT_VF_LEN] = {""};   /* virtual input  filename */
    char out[GMT_VF_LEN] = {""};   /* virtual output filename */

    /* Virtual file to hold the input grid*/
    GMT_Open_VirtualFile(API,
            GMT_IS_GRID,
            GMT_IS_SURFACE,
            GMT_IN,
            G,inp);
    //fprintf(stderr,"%s Name of virtual input grid: (%s)\n",ctimer(),inp);
    /* Virtual file to hold the resulting grid*/
    GMT_Open_VirtualFile(API,
            GMT_IS_GRID,
            GMT_IS_SURFACE,
            GMT_OUT,
            NULL,out);
    /* Prepare the module arguments */
    sprintf(args,"-R%.14g/%.14g/%.14g/%.14g -I%s %s -G%s",
            F->xTL,
            (F->xTL+G->header->inc[0]*(F->nx-1)),
            (F->yTL-G->header->inc[1]*(F->ny-1)),
            F->yTL,
            ires,
            inp,out);
    //printf("%s\n",args);
    /* Call the grdproject module */
    GMT_Call_Module (API,"grdsample",GMT_MODULE_CMD,args);
    /* Obtain the grid from the virtual file */
    struct GMT_GRID *Z = GMT_Read_VirtualFile (API,out);
    if(ires[strlen(ires)-1]=='s')
    {
        ires[strlen(ires)-1]='\0';
        Z->header->inc[0] = Z->header->inc[1] = atof(ires)/3600.0;
    }
    else if(ires[strlen(ires)-1]=='m')
    {
        ires[strlen(ires)-1]='\0';
        Z->header->inc[0] = Z->header->inc[1] = atof(ires)/60.0;
    }
    else
    {
        ires[strlen(ires)-1]='\0';
        Z->header->inc[0] = Z->header->inc[1] = atof(ires);
    }
    /* Close the virtual files */
    GMT_Close_VirtualFile(API,inp);
    GMT_Close_VirtualFile(API,out);
    if(Vflag)
    {
        fprintf(stderr,"----------------------------------------------------------------------------------------------\n");
        fprintf(stderr,"%s Inner most grid file has been successfuly resampled\n",ctimer());
        fprintf(stderr,"%s %s registration\n",ctimer(),(Z->header->registration==0)?"Gridline":"Pixel");
        fprintf(stderr,"%s xmin = %.15g xmax = %.15g xinc = %.15g n_xcols = %d\n",ctimer(),Z->header->wesn[0],Z->header->wesn[1],Z->header->inc[0],Z->header->n_columns);
        fprintf(stderr,"%s ymin = %.15g ymax = %.15g yinc = %.15g n_yrows = %d\n",ctimer(),Z->header->wesn[2],Z->header->wesn[3],Z->header->inc[1],Z->header->n_rows);
        fprintf(stderr,"----------------------------------------------------------------------------------------------\n");
    }
    free(F);

    return Z;
}

void print_info_GRID(struct GMT_GRID *G)
{
    printf("%s\n",G->header->title);
    printf("Registration: %s\n",(G->header->registration==1)?"Gridline":"Pixel");
    printf("xmin = %.15g xmax = %.15g xinc = %.15g n_xcols = %d\n",
            G->header->wesn[0],G->header->wesn[1],
            G->header->inc[0],G->header->n_columns);
    printf("ymin = %.15g ymax = %.15g yinc = %.15g n_yrows = %d\n",
            G->header->wesn[2],G->header->wesn[3],
            G->header->inc[1],G->header->n_rows);
}

