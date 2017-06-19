#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datastream.h"
#include "format.h"
#include "trajdata.h"
#include "charmm.h"
#include "amber.h"

void
clearTrajData(TrajectoryData *t)
{
    t->frames=0;
    t->nAtoms=0;
    t->startTime=0.0;
    t->startFrame=0;
    t->timeStep=0.0;
    t->periodicBoundariesQ=0;
    t->fixedAtoms=0;
    t->fluctChargesQ=0;
    t->degreesOfFreedom=0;
    t->fourdimensQ=0;
    t->velocitiesQ=0;
    t->versionnum=0;
    *t->version=0;
    t->title=0;
    t->titleLength=-1;
    t->titleLines=0;
    t->fixedIndex=0;
} //End clearTrajData()



void
showTrajData(TrajectoryData *t)
{
    int i;

    printf("== Trajectory information ====================\n");

    if (t->titleLines>0 && t->title!=(char **)0)
    {
        printf("title:\n");
        for (i=0; i<t->titleLines; i++)
        {
            printf("  >%s<\n",t->title[i]);
        }
    }
    printf("atoms                  : %d\n",t->nAtoms);
    printf("frames                 : %d\n",t->frames);
    printf("start time             : %f\n",t->startTime);
    printf("start frame            : %d\n",t->startFrame);
    printf("time step              : %1.15f\n",t->timeStep);
    printf("degrees of freedom     : %d\n",t->degreesOfFreedom);
    printf("periodic boundaries?   : %d\n",t->periodicBoundariesQ);
    printf("fixed atoms            : %d\n",t->fixedAtoms);
    printf("fluctuating charges?   : %d\n",t->fluctChargesQ);
    printf("four dimensional?      : %d\n",t->fourdimensQ);
    printf("velocities?            : %d\n",t->velocitiesQ);
    if (*t->version!=0)
        printf("version                : %s\n",t->version);
    if (t->versionnum>0)
        printf("version number         : %d\n",t->versionnum);
} //End showTrajData()



int
readTrajData(DataStream *io, Format *fmt, TrajectoryData *td)
{
    int retval = 0; 

    switch (fmt->trajFormat)
    {
    case CHARMM:
        retval = readHeaderCHARMM(io,fmt,td);
        break;

    case Amber:
      retval = readHeaderAmber(io,fmt,td);
      break;
      
    case NOFORMAT:
      fprintf(stderr,"Unknown format. Cannot read data.\n");
      break;
      
    default:
      fprintf(stderr,"Invalid format. Cannot read data.\n");
      break;
    }

    return retval;

} //End readTrajData()


int
writeHeader(DataStream *io, Format *fmt, TrajectoryData *t)
{
    int retval = 0;

    switch (fmt->trajFormat)
    {
    case CHARMM:
        retval = writeHeaderCHARMM(io,fmt,t);
        break;

    case Amber:
        retval = writeHeaderAmber(io,fmt,t);
        break;

    case NOFORMAT:
        fprintf(stderr,"Unknown format. Cannot write data.\n");
        break;

    default:
        fprintf(stderr,"Invalid format. Cannot write data.\n");
        break;
    }

} //End writeHeader()



Frame *
allocateFrame(int n)
{
    int i;
    Frame *f;

    if ((f=(Frame *)malloc(sizeof(Frame)))==0)
    {
        fprintf(stderr,"cannot allocate storage for frame\n");
        return NULL;
    }

    f->nAtoms=n;

    if ((f->x=(double *)malloc(sizeof(double)*n))==0)
    {
        fprintf(stderr,"cannot allocate storage for frame coordinate\n");
        return NULL;
    }

    if ((f->y=(double *)malloc(sizeof(double)*n))==0)
    {
        fprintf(stderr,"cannot allocate storage for frame coordinate\n");
        return NULL;
    }

    if ((f->z=(double *)malloc(sizeof(double)*n))==0)
    {
        fprintf(stderr,"cannot allocate storage for frame coordinate\n");
        return NULL;
    }

    clearFrame(f);

    return f;

} //End allocateFrame()



void
clearFrame(Frame *f)
{
    int i;

    for (i=0; i<f->nAtoms; f->x[i++]=0.0);
    for (i=0; i<f->nAtoms; f->y[i++]=0.0);
    for (i=0; i<f->nAtoms; f->z[i++]=0.0);

    f->boxX=0.0;
    f->boxY=0.0;
    f->boxZ=0.0;
    f->angleA=90.0;
    f->angleB=90.0;
    f->angleC=90.0;
} //End clearFrame()



void
deallocateFrame(Frame *f)
{
    free(f->x);
    free(f->y);
    free(f->z);
    free(f);
} //End deallocateFrame()



void
copyFrameSelection(Frame *from, Frame *to, int *inx, int n, int *tinx, int tn)
{
    int i;
    int ix;

    if (to->nAtoms!=n)
    {
        fprintf(stderr,"size of destination frame does not match\n");
        return;
    }

    for(i=0; i<n; i++)
    {
	if (tn==0)
	  ix=inx[i]-1;
	else 
	  ix=tinx[inx[i]-1]-1;

        if (ix<0 || ix>=from->nAtoms)
        {
            fprintf(stderr,"index %d out of range during copy\n",ix);
            return;
        }

        to->x[i]=from->x[ix];
	to->y[i]=from->y[ix];
	to->z[i]=from->z[ix];
    }

    to->boxX=from->boxX;
    to->boxY=from->boxY;
    to->boxZ=from->boxZ;
    to->angleA=from->angleA;
    to->angleB=from->angleB;
    to->angleC=from->angleC;
} //End copyFrameSelection()


void
copyAvgFrame(Frame *from, Frame *to, int navg, int *nsel, int **afrom, int **ato) {
    int i,j,k;
    
    for (i=0; i<navg; i++) {
      double cx=0.0;
      double cy=0.0;
      double cz=0.0;
      int nn=0;
      for (j=0; j<nsel[i]; j++) {
        for (k=afrom[i][j]; k<=ato[i][j]; k++) {
          int ix=k-1;
          cx+=from->x[ix];
          cy+=from->y[ix];
          cz+=from->z[ix];
          nn++;
        }
      }
      cx/=(double)nn; 
      cy/=(double)nn; 
      cz/=(double)nn; 
      to->x[i]=cx;
      to->y[i]=cy;
      to->z[i]=cz;
    }

    to->boxX=from->boxX;
    to->boxY=from->boxY;
    to->boxZ=from->boxZ;
    to->angleA=from->angleA;
    to->angleB=from->angleB;
    to->angleC=from->angleC;
  
}


void
copyFrame(Frame *from, Frame *to, int *tinx, int tn)
{
    int i;
    int ix;

    if (to->nAtoms<from->nAtoms)
    {
        fprintf(stderr,"size of source and destination frames should match\n");
        return;
    }

    for(i=0; i<from->nAtoms; i++)
    {
      if (tn==0)
         ix=i;
      else 
	 ix=tinx[i]-1;

      to->x[i]=from->x[ix];
      to->y[i]=from->y[ix];
      to->z[i]=from->z[ix];
    }
    
    for(i=from->nAtoms; i<to->nAtoms; i++) 
    {
      to->x[i]=0.0;
      to->y[i]=0.0;
      to->z[i]=0.0;
    }

    to->boxX=from->boxX;
    to->boxY=from->boxY;
    to->boxZ=from->boxZ;
    to->angleA=from->angleA;
    to->angleB=from->angleB;
    to->angleC=from->angleC;

} //End copyFrame()

void
addFrame(Frame *from, Frame *to)
{
    int i;

    if (to->nAtoms!=from->nAtoms)
    {
        fprintf(stderr,"size of source and destination frames should match\n");
        return;
    }

    for(i=0; i<from->nAtoms; i++)
    {
      to->x[i]+=from->x[i];
      to->y[i]+=from->y[i];
      to->z[i]+=from->z[i];
    }
    
    to->boxX+=from->boxX;
    to->boxY+=from->boxY;
    to->boxZ+=from->boxZ;
    to->angleA+=from->angleA;
    to->angleB+=from->angleB;
    to->angleC+=from->angleC;
} //End addFrame()

void
scaleFrame(Frame *f, double factor) 
{
     int i;
     for (i=0; i<f->nAtoms; i++) 
     {
        f->x[i]*=factor;
        f->y[i]*=factor;
        f->z[i]*=factor;
     }
     
     f->boxX*=factor;
     f->boxY*=factor;
     f->boxZ*=factor;
     f->angleA*=factor;
     f->angleB*=factor;
     f->angleC*=factor;
} //End scaleFrame() 


void
removedriftFrame(Frame *f)
{
     double cx,cy,cz;
     int i;

     cx=cy=cz=0.0;

     for (i=0; i<f->nAtoms; i++) {
       cx+=f->x[i];
       cy+=f->y[i];
       cz+=f->z[i];
     }

     cx/=(double)f->nAtoms;  
     cy/=(double)f->nAtoms;  
     cz/=(double)f->nAtoms;  

     for (i=0; i<f->nAtoms; i++) {
       f->x[i]-=cx;
       f->y[i]-=cy;
       f->z[i]-=cz;
     }
} 

void
mergeFrame(Frame *to, Frame *merge, int offset)
{
    int i;

    for(i=0; i<merge->nAtoms && i+offset<to->nAtoms; i++)
    {
      to->x[i+offset]=merge->x[i];
      to->y[i+offset]=merge->y[i];
      to->z[i+offset]=merge->z[i];
    }
}


void
showFrame(Frame *f, int inx)
{
    int i;

    fprintf(stdout,"== Frame data %8d =========================\n",inx);
    fprintf(stdout,"periodic box: %1.5f x %1.5f x %1.5f ( %1.5f x %1.5f x %1.5f )\n",
            f->boxX,f->boxY,f->boxZ,f->angleA,f->angleB,f->angleC);
    for (i=0; i<f->nAtoms; i++)
        fprintf(stdout,"coor %7d %14.6f %14.6f %14.6f\n",i+1,f->x[i],f->y[i],f->z[i]);
} //End showFrame()



int
readFrame(DataStream *io, Format *fmt, TrajectoryData *td, Frame *f)
{
    int retval = 0;

    switch (fmt->trajFormat)
    {
    case NOFORMAT:
        fprintf(stderr,"input format undefined\n");
        return -1;

    case XYZ:
        //return readFrameXYZ(io,fmt,td,f);
        break;

    case PDB:
        //return readFramePDB(io,fmt,td,f);
        break;

    case CHARMM:
        retval = readFrameCHARMM(io,fmt,td,f);
        break;

    case Amber:
        retval = readFrameAmber(io,fmt,td,f);
        break;

    default:
        fprintf(stderr,"input format unknown\n");
    }

    return retval;

} //End readFrame()



/* write frame to output stream */
int
writeFrame(DataStream *io, Format *fmt, TrajectoryData *td, Frame *f)
{
    int retval = 0;

    switch (fmt->trajFormat)
    {
    case NOFORMAT:
        fprintf(stderr,"output format undefined\n");
        return -1;

    case XYZ:
        //writeFrameXYZ(io,fmt,td,f);
        break;

    case PDB:
        //writeFramePDB(io,fmt,td,f);
        break;

    case CHARMM:
        retval = writeFrameCHARMM(io,fmt,td,f);
        break;

    case Amber:
        retval = writeFrameAmber(io,fmt,td,f);
        break;

    default:
        fprintf(stderr,"output format unknown\n");
        return retval;
    }

    return retval;

} //End writeFrame()


/* set trajectory data from user input */
void 
setTrajDataFromInput(TrajectoryData *t, TrajectoryData *d) 
{
  if (d->frames>0)       
    t->frames=d->frames;
  if (d->nAtoms>0)
    t->nAtoms=d->nAtoms;
  if (d->startTime>0.0)
    t->startTime=d->startTime;
  if (d->startFrame>0) 
    t->startFrame=d->startFrame;
  if (d->timeStep>0.0)
    t->timeStep=d->timeStep;
  if (d->periodicBoundariesQ>0) 
    t->periodicBoundariesQ=d->periodicBoundariesQ;
  if (t->timeStep>0.0 && t->startFrame>0) 
    t->startTime=d->startFrame*d->timeStep;
} //End setTrajDataFromInput

void
parseTrajDataString(char  *str, TrajectoryData *d) 
{
    char *cptr;
    int equal;
    char fstr[128];
    char tstr[128];
    int i,fi,ti;

    cptr=str;

    do
    {
      equal=0;
        fi=0;
        ti=0;
        for (i=0; cptr[i]!=0 && cptr[i]!=','; i++)
        {
            if (cptr[i]=='=')
                equal=i;
            else
                if (equal==0)
                    fstr[fi++]=cptr[i];
                else
                    tstr[ti++]=cptr[i];
        }
        fstr[fi]=0;
        tstr[ti]=0;

        if (!strcmp(fstr,"atoms"))
        {
	  d->nAtoms=atoi(tstr);
        }
        else if (!strcmp(fstr,"frames"))
        {
	  d->frames=atoi(tstr);
        }
        else if (!strcmp(fstr,"delta") || !strcmp(fstr,"timestep"))
        {
	  d->timeStep=atof(tstr);
        }
        else if (!strcmp(fstr,"starttime")) 
        {
	  d->startTime=atof(tstr);
        }
        else if (!strcmp(fstr,"startframe"))
        {
	  d->startFrame=atoi(tstr);
        }
        else if (!strcmp(fstr,"periodic"))
        {
	  d->periodicBoundariesQ=atoi(tstr);
        }

        for (cptr+=i; *cptr==','; cptr++);
    }
    while(*cptr!=0);
} //End parseTrajDataString()


long
calcFrameSize(Format *fmt, TrajectoryData *td)
{
    long retval = 0;

    switch (fmt->trajFormat)
    {
        case NOFORMAT:
            fprintf(stderr,"output format undefined\n");
            return retval;

        case XYZ:
            //retval = calcFrameSizeXYZ(fmt,td);
            break;

        case PDB:
            //retval = calcFrameSizePDB(fmt,td);
            break;

        case CHARMM:
            retval = calcFrameSizeCHARMM(fmt,td);
            break;

        case Amber:
            retval = calcFrameSizeAmber(fmt,td);
            break;

        default:
            fprintf(stderr,"output format unknown\n");
            return retval;
    }

    return retval;
    
} //End calcFrameSize()


long
calcDataSize(Format *fmt, TrajectoryData *td)
{
    long retval = 0;

    switch (fmt->trajFormat)
    {
    case NOFORMAT:
        fprintf(stderr,"output format undefined\n");
        return retval;

    case XYZ:
        //retval = calcDataSizeXYZ(fmt,td);
        break;

    case PDB:
        //retval = calcDataSizePDB(fmt,td);
        break;

    case CHARMM:
        retval = calcDataSizeCHARMM(fmt,td);
        break;

    case Amber:
        retval = calcDataSizeAmber(fmt,td);
        break;

    default:
        fprintf(stderr,"output format unknown\n");
        return retval;
    }

    return retval;

} //End calcDataSize()

