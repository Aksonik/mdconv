#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datastream.h"
#include "format.h"
#include "trajdata.h"
#include "charmm.h"

int
testHeaderCHARMM(DataStream *io, Format *fmt)
{
    char buffer[1024];
    int len;
    int i;

    rewindBuffer(io);
    int n=readData(io,buffer,256);

    if (!strncmp(buffer+4,"CORD",4) || !strncmp(buffer+4,"DROC",4))
    {
        fmt->fortranHeaderSize=4;
    }
    else if (!strncmp(buffer+8,"CORD",4) || !strncmp(buffer+8,"DROC",4))
    {
        fmt->fortranHeaderSize=8;
    }
    else
    {
        return 0;
    }

    if (buffer[0] != 0)
    {
        fmt->endian=LITTLEENDIAN;
    }
    else
    {
        fmt->endian=BIGENDIAN;
    }

    convertIntegerFromBuffer(&len,1,buffer,fmt->fortranHeaderSize,fmt->endian);
    if (len==84)
    {
        fmt->intSize=4;
    }
    else if (len==164)
    {
        fmt->intSize=8;
    }
    else
    {
        fprintf(stderr,"unexpected CHARMM icontrol size: %d\n",len);
        return -1;
    }

    fmt->realSizeXYZ=4;
    fmt->realSizeOther=8;

    return 1;

} //End testHeaderCHARMM()



int
readHeaderCHARMM(DataStream *io, Format *fmt, TrajectoryData *t)
{
    int icontrol[20];
    int len;
    int j,it;
    char *buffer;
    int buffersize;
    int tlenSize;
    int actualtitlelength;
    int copytitlelength;

    buffersize=8192;
    if ((buffer=(char *)malloc(buffersize))==0)
    {
        fprintf(stderr,"cannot allocate storage\n");
        return -1;
    }

    // CORD and icontrol
    if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -2;
    convertIntegerFromBuffer(icontrol,20,buffer+4,fmt->intSize,fmt->endian);

    t->frames=icontrol[0];
    t->degreesOfFreedom=icontrol[7];
    t->fixedAtoms=icontrol[8];
    t->periodicBoundariesQ=icontrol[10];
    t->velocitiesQ=icontrol[4];
    t->fourdimensQ=icontrol[11];
    t->fluctChargesQ=icontrol[12];
    t->versionnum=icontrol[19];

    if (icontrol[2]==0) {
      icontrol[2]=1;
    }

    convertRealFromBuffer(&t->timeStep,1,buffer+4+fmt->intSize*9,fmt->realSizeXYZ,myEndian());
    t->timeStep*=4.88882129E-02;
    t->timeStep*=(double)icontrol[2];

    it=(int)(t->timeStep*100000.0+0.5);
    t->timeStep=(double)it/100000.0;

    t->startTime=(icontrol[1]/icontrol[2])*t->timeStep;

    if (t->velocitiesQ)
        fprintf(stderr,"WARNING: Code to handle velocities not yet implemented\n");

    if (t->fourdimensQ)
        fprintf(stderr,"WARNING: Code to handle 4D data not yet implemented\n");

    if (t->fluctChargesQ)
        fprintf(stderr,"WARNING: Code to handle fluctuating charges not yet implemented\n");

    //Title
    if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -3;
    tlenSize=4;
    convertIntegerFromBuffer(&t->titleLines,1,buffer,tlenSize,fmt->endian);

    if (t->titleLines>100)
    {
        fprintf(stderr,"too many title lines (%d)\n",t->titleLines);
        return -4;
    }

    t->title=(char **)malloc(sizeof(char *)*t->titleLines);
    if (t->title==0)
    {
        fprintf(stderr,"cannot allocate storage for title\n");
        return -5;
    }

    t->titleLength=80;
    actualtitlelength=(len-tlenSize)/t->titleLines;
    copytitlelength=(actualtitlelength>80)?80:actualtitlelength;
    for (j=0; j<t->titleLines; j++)
    {
        t->title[j]=(char *)malloc(t->titleLength+1);
        if (t->title[j]==0)
        {
            fprintf(stderr,"cannot allocate storage for title line\n");
            return -6;
        }
        strncpy(t->title[j],buffer+tlenSize+j*actualtitlelength,copytitlelength);
        *(t->title[j]+copytitlelength)=0;
    }

    //NATOM
    if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -7;
    convertIntegerFromBuffer(&t->nAtoms,1,buffer,fmt->intSize,fmt->endian);
    //fixed atoms
    if (t->fixedAtoms>0)
    {
        buffersize=t->nAtoms*8;
        if ((buffer=(char *)realloc(buffer,buffersize))==0)
        {
            fprintf(stderr,"cannot reallocate storage\n");
            return -8;
        }

        if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -9;
        t->fixedIndex=(int *)malloc(sizeof(int)*(t->nAtoms-t->fixedAtoms));
        convertIntegerFromBuffer(t->fixedIndex,t->nAtoms-t->fixedAtoms,buffer,fmt->intSize,fmt->endian);

        t->refFrame=allocateFrame(t->nAtoms);
    }

    free(buffer);

    return 1;

} //End readHeaderCHARMM()



int
readFrameCHARMM(DataStream *io, Format *fmt,  TrajectoryData *t, Frame *f)
{
    char *buffer;
    int buffersize;
    int len,i;
    double *tmparr;
    double v[6];

    buffersize=t->nAtoms*sizeof(double);
    if (t->periodicBoundariesQ && buffersize<48)
      buffersize=48;
    
    if ((buffer=(char *)malloc(buffersize))==0)
    {
        fprintf(stderr,"cannot allocate storage\n");
        return -1;
    }
    
    if (t->periodicBoundariesQ)
    {
        if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -1;
        convertRealFromBuffer(v,6,buffer,fmt->realSizeOther,fmt->endian);

        f->boxX=sqrt(v[0]*v[0]+v[1]*v[1]+v[3]*v[3]);
	f->boxY=sqrt(v[1]*v[1]+v[2]*v[2]+v[4]*v[4]);
	f->boxZ=sqrt(v[3]*v[3]+v[4]*v[4]+v[5]*v[5]);

	f->angleA=acos((v[4]*(v[2]+v[5])+v[1]*v[3])/(f->boxY*f->boxZ))*180.0/M_PI;
	f->angleB=acos((v[3]*(v[0]+v[5])+v[1]*v[4])/(f->boxX*f->boxZ))*180.0/M_PI;
	f->angleC=acos((v[1]*(v[0]+v[2])+v[3]*v[4])/(f->boxX*f->boxY))*180.0/M_PI;
    }

    if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -1;
    len/=fmt->realSizeXYZ;
    if (len == t->nAtoms)
    {
        convertRealFromBuffer(f->x,t->nAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
        if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -1;
        convertRealFromBuffer(f->y,t->nAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
        if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -1;
        convertRealFromBuffer(f->z,t->nAtoms,buffer,fmt->realSizeXYZ,fmt->endian);

        if (t->fixedAtoms>0)
        {
            memcpy(t->refFrame->x,f->x,sizeof(double)*t->nAtoms);
            memcpy(t->refFrame->y,f->y,sizeof(double)*t->nAtoms);
            memcpy(t->refFrame->z,f->z,sizeof(double)*t->nAtoms);
        }
    }
    else if (len == t->nAtoms-t->fixedAtoms)
    {
        memcpy(f->x,t->refFrame->x,sizeof(double)*t->nAtoms);
        memcpy(f->y,t->refFrame->y,sizeof(double)*t->nAtoms);
        memcpy(f->z,t->refFrame->z,sizeof(double)*t->nAtoms);

        if ((tmparr=(double *)malloc(sizeof(double)*(t->nAtoms-t->fixedAtoms)))==0)
        {
            fprintf(stderr,"cannot allocate storage\n");
            return -1;
        }

        convertRealFromBuffer(tmparr,t->nAtoms-t->fixedAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
        for (i=0; i<t->nAtoms-t->fixedAtoms; i++) f->x[t->fixedIndex[i]-1]=tmparr[i];

        if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -1;
        convertRealFromBuffer(tmparr,t->nAtoms-t->fixedAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
        for (i=0; i<t->nAtoms-t->fixedAtoms; i++) f->y[t->fixedIndex[i]-1]=tmparr[i];

        if ((len=readFortran(io,fmt,buffer,buffersize))<0) return -1;
        convertRealFromBuffer(tmparr,t->nAtoms-t->fixedAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
        for (i=0; i<t->nAtoms-t->fixedAtoms; i++) f->z[t->fixedIndex[i]-1]=tmparr[i];

        free(tmparr);
    }

    free(buffer);

    return 1;

} //End readFrameCHARMM()



int
writeHeaderCHARMM(DataStream *io, Format *fmt, TrajectoryData *t)
{
    int icontrol[20];
    int len;
    int i,j;
    char *buffer;
    int buffersize;
    double tstep;
    int tlines;
    int tlen;
    char *ttit;
    int tlenSize;

    buffersize=(t->nAtoms*16<8192)?8192:t->nAtoms*16;
    if ((buffer=(char *)malloc(buffersize))==0)
    {
        fprintf(stderr,"cannot allocate storage\n");
        return -1;
    }

    // CORD and icontrol
    strcpy(buffer,"CORD");
    for (i=0; i<20; icontrol[i++]=0);
    icontrol[0]=t->frames;
    if (t->timeStep<0.00001) {
      icontrol[1]=0;
    } else {
      icontrol[1]=(int)(t->startTime/t->timeStep+0.5);
    }
    icontrol[2]=1;
    icontrol[3]=t->frames;
    icontrol[7]=(t->nAtoms*3)-6;
    icontrol[10]=t->periodicBoundariesQ;
    icontrol[19]=27;

    convertIntegerToBuffer(icontrol,20,buffer+4,fmt->intSize,fmt->endian);

    tstep=t->timeStep/4.88882129E-02;
    convertRealToBuffer(&tstep,1,buffer+4+9*fmt->intSize,fmt->realSizeXYZ,fmt->endian);

    if ((len=writeFortran(io,fmt,buffer,fmt->intSize*20+4))<0) return -1;

    tlenSize=4;
    tlines=(t->titleLines==0)?1:t->titleLines;
    convertIntegerToBuffer(&tlines,1,buffer,tlenSize,fmt->endian);

    //Title
    if (t->titleLines==0)
    {
        tlen=80;
        sprintf(buffer+tlenSize,"* Generated by MDCONV");
        for (j=strlen(buffer+tlenSize); j<tlen; j++)
            *(buffer+tlenSize+j)=' ';
    }
    else
    {
        tlen=t->titleLength;
        for (j=0; j<t->titleLines; j++)
            strncpy(buffer+tlenSize+j*t->titleLength,t->title[j],t->titleLength);
    }

    if ((len=writeFortran(io,fmt,buffer,tlenSize+tlen*tlines))<0) return -1;

    //NATOM
    convertIntegerToBuffer(&t->nAtoms,1,buffer,fmt->intSize,fmt->endian);
    if ((len=writeFortran(io,fmt,buffer,fmt->intSize))<0) return -1;

    free(buffer);

    return 1;

} //End writeHeaderCHARMM()



int
writeFrameCHARMM(DataStream *io, Format *fmt, TrajectoryData *t, Frame *f)
{
    char *buffer;
    int buffersize;
    int len,i;
    double v[6];

    buffersize=(f->nAtoms<6?6:f->nAtoms)*sizeof(double);
    if ((buffer=(char *)malloc(buffersize))==0)
    {
        fprintf(stderr,"cannot allocate storage\n");
        return -1;
    }

    if (t->periodicBoundariesQ)
    {
        if (fabs(f->angleA-90.0)<0.01 && fabs(f->angleB-90.0)<0.01 && fabs(f->angleC-90.0)<0.01) {
          v[0]=f->boxX;
          v[1]=0.0;
          v[2]=f->boxY;
          v[3]=0.0;
          v[4]=0.0;
          v[5]=f->boxZ;
        } else if (fabs(f->angleA-109.471)<0.01 && fabs(f->angleB-109.471)<0.01 && 
                   fabs(f->angleC-109.471)<0.01 && 
                   fabs(f->boxX-f->boxY)<0.01 && fabs(f->boxX-f->boxZ)<0.01) {
          v[0]=f->boxX*sqrt(1.0-2.0/9.0/3.0);
          v[1]=-f->boxX/3.0/sqrt(3.0);
          v[2]=f->boxX*sqrt(1.0-2.0/9.0/3.0);
          v[3]=-f->boxX/3.0/sqrt(3.0);
          v[4]=-f->boxX/3.0/sqrt(3.0);
          v[5]=f->boxX*sqrt(1.0-2.0/9.0/3.0);
        } else {
          v[0]=v[1]=v[2]=v[3]=v[4]=v[5]=0.0;
        }
	
        convertRealToBuffer(v,6,buffer,fmt->realSizeOther,fmt->endian);
        if ((len=writeFortran(io,fmt,buffer,fmt->realSizeOther*6))<0) return -1;
    }

    convertRealToBuffer(f->x,f->nAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
    if ((len=writeFortran(io,fmt,buffer,fmt->realSizeXYZ*t->nAtoms))<0) return -1;
    convertRealToBuffer(f->y,f->nAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
    if ((len=writeFortran(io,fmt,buffer,fmt->realSizeXYZ*t->nAtoms))<0) return -1;
    convertRealToBuffer(f->z,f->nAtoms,buffer,fmt->realSizeXYZ,fmt->endian);
    if ((len=writeFortran(io,fmt,buffer,fmt->realSizeXYZ*t->nAtoms))<0) return -1;
    
    free(buffer);

    return 1;

} //End writeFrameCHARMM()


long
calcFrameSizeCHARMM(Format *fmt, TrajectoryData *td)
{
    long frame_size;

    if( (fmt == NULL) || (td == NULL) )
    {
        return -1;
    }

    // Size of data
    frame_size = 3 * (td->nAtoms - td->fixedAtoms) * fmt->realSizeXYZ;
    frame_size += (3 * 2 * fmt->fortranHeaderSize); // Size of fortran headers

    // Add size to each frame used if pboundaries is on
    if(td->periodicBoundariesQ == 1)
    {
        frame_size += (6 * fmt->realSizeOther) + (2 * fmt->fortranHeaderSize);
    }

    return frame_size;

} //End calcFrameSizeCHARMM()


long
calcDataSizeCHARMM(Format *fmt, TrajectoryData *td)
{
    long data_size;
    long frame_size;

    if( (fmt == NULL) || (td == NULL) )
    {
        return -1;
    }

    frame_size = calcFrameSizeCHARMM(fmt, td);
    data_size = frame_size * td->frames;

    data_size += 4 + (2 * fmt->fortranHeaderSize) +
                (80 * td->titleLines);  // Title
    data_size += 164 + (2 * fmt->fortranHeaderSize); // Header info
    data_size += fmt->intSize + (2 * fmt->fortranHeaderSize); // Num atoms

    return data_size;

} //End calcDataSizeCHARMM()

