#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "amber.h"
#include "datastream.h"
#include "format.h"
#include "trajdata.h"



int
testHeaderAmber(DataStream *io, Format *fmt)
{
    return 0;
} //End testHeaderAmber()



int
readHeaderAmber(DataStream *io, Format *fmt, TrajectoryData *t)
{
  char *buffer;
  int buffersize;
  int len;
  off_t filesize;
  off_t framesize;
  off_t oframes;
  int ncoord,nlines;
  int nextra;

  buffersize=1024;
  if ((buffer=(char *)malloc(buffersize))==0) {
    fprintf(stderr,"cannot allocate storage\n");
    return -1;
  }

  len=readLine(io,buffer,buffersize);

  t->titleLines=1;
  t->title=(char **)malloc(sizeof(char *)*t->titleLines);
  t->titleLength=len;
  t->title[0]=(char *)malloc(t->titleLength+1);
  strncpy(t->title[0],buffer,t->titleLength);
  *(t->title[0]+t->titleLength)=0;

  filesize=getSize(io);

  if (filesize>0) {
    ncoord=t->nAtoms*3;
    nlines=(int)(ncoord/10);
    nextra=ncoord%10;
    framesize=(off_t)(nlines*81+nextra*8+1);

    if (t->periodicBoundariesQ) {
      framesize+=(off_t)(3*8+1);
    }
    oframes=(off_t)(filesize-(off_t)(len+1))/framesize;
    t->frames=(int)(oframes);
  }

  free(buffer);

  return 1;
} //End readHeaderAmber()



int
readFrameAmber(DataStream *io, Format *fmt,  TrajectoryData *t, Frame *f)
{
  char *buffer;
  int buffersize;
  int len,i,j,n;
  int fulllines;

  double *values;

  buffersize=82;
  if ((buffer=(char *)malloc(buffersize))==0) {
    fprintf(stderr,"cannot allocate storage\n");
    return -1;
  }

  if ((values=(double *)malloc(3*t->nAtoms*sizeof(double)))==0) {
    fprintf(stderr,"cannot allocate storage\n");
    return -1;
  }

  n=0;
  fulllines=(int)(t->nAtoms*3/10);
  for (i=0; i<fulllines; i++) {
    len=readLine(io,buffer,buffersize);
    for (j=0; j<10; j++) {
      values[n++]=getFixedNumber(buffer,j,8);
    }
  }

  if ((t->nAtoms*3)%10!=0) {
    len=readLine(io,buffer,buffersize);
    for (j=0; j<(t->nAtoms*3)%10; j++) {
      values[n++]=getFixedNumber(buffer,j,8);
    }
  }

  n=0;
  for (i=0; i<t->nAtoms; i++) {
    f->x[i]=values[n++];
    f->y[i]=values[n++];
    f->z[i]=values[n++];
  }

  if (t->periodicBoundariesQ) {
    len=readLine(io,buffer,buffersize);
    f->boxX=getFixedNumber(buffer,0,8);
    f->boxY=getFixedNumber(buffer,1,8);
    f->boxZ=getFixedNumber(buffer,2,8);
    f->angleA=90.0;
    f->angleB=90.0;
    f->angleC=90.0;
  }

  free(values);
  free(buffer);

  return 1;
} //End readFrameAmber()


int
writeHeaderAmber(DataStream *io, Format *fmt, TrajectoryData *t)
{
    char title[81];
    int i;

    if (t->titleLines>0)
    {
        strncpy(title,t->title[0],80);
    }
    else
    {
        strcpy(title,"generated with MDCONV");
    }

    for (i=0; i<80 && title[i]!=0; i++);
    for (;i<80; title[i]=' ', i++); 
    title[80]='\n';

    writeData(io,title,81);

    return 1;

} //End writeHeaderAmber()



int
writeFrameAmber(DataStream *io, Format *fmt, TrajectoryData *t, Frame *f)
{
    char buffer[82];
    char *bptr;
    int i,nn;

    bptr=buffer;
    nn=0;
    for (i=0; i<f->nAtoms; i++)
    {
        sprintf(bptr,"%8.3f",f->x[i]);
        bptr+=8;
        if (++nn%10==0)
        {
            *bptr='\n';
            writeData(io,buffer,81);
            bptr=buffer;
        }

        sprintf(bptr,"%8.3f",f->y[i]);
        bptr+=8;
        if (++nn%10==0)
        {
            *bptr='\n';
            writeData(io,buffer,81);
            bptr=buffer;
        }      

        sprintf(bptr,"%8.3f",f->z[i]);
        bptr+=8;
        if (++nn%10==0)
        {
            *bptr='\n';
            writeData(io,buffer,81);
            bptr=buffer;
        }      
    }  

    if (nn%10>0)
    {
        *bptr='\n';
        writeData(io,buffer,(nn%10)*8+1);
    }

    if (t->periodicBoundariesQ)
    {
        sprintf(buffer,"%8.3f%8.3f%8.3f\n",f->boxX,f->boxY,f->boxZ);
        writeData(io,buffer,25);
    }

    return 1;
} //End writeFrameAmber()


long
calcFrameSizeAmber(Format *fmt, TrajectoryData *td)
{
    long frame_size;
    long ncoord;

    if( (fmt == NULL) || (td == NULL) )
    {
        return -1;
    }

    ncoord = td->nAtoms * 3;
    frame_size = (long)(ncoord/10)*81 + (ncoord%10)*8+1;

    if(td->periodicBoundariesQ == 1)
    {
        frame_size += (3*8+1);
    }

    return frame_size;

} //End calcFrameSizeAmber()


long
calcDataSizeAmber(Format *fmt, TrajectoryData *td)
{
    long data_size;

    if( (fmt == NULL) || (td == NULL) )
    {
        return -1;
    }

    // header
    data_size=81;
    
    // data
    data_size += td->frames * calcFrameSizeAmber(fmt, td);

    return data_size;

} //End calcDataSizeAmber()

