#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datastream.h"
#include "format.h"
#include "trajdata.h"
#include "charmm.h"
#include "amber.h"

Endian myendian=NOENDIAN;

void
clearFormat(Format *f)
{
    f->trajFormat=NOFORMAT;
    f->realSizeXYZ=0;
    f->realSizeOther=0;
    f->intSize=0;
    f->endian=NOENDIAN;
    f->fortranHeaderSize=0;
} //End clearFormat()



int
identifyFormat(DataStream *io, Format *f)
{
    int ret=0;

    clearFormat(f);

    fillBuffer(io,1024);

    if ((ret=testHeaderCHARMM(io,f))>0)
    {
      f->trajFormat=CHARMM;
    } else if ((ret=testHeaderAmber(io,f))>0) {
      f->trajFormat=Amber;
    }
    else if (ret<0)
    {
        fprintf(stderr,"error reading file\n");
        return -1;
    }
    else
    {
      //fprintf (stderr,"cannot identify trajectory format\n");
    }

    rewindBuffer(io);

    return ret;

} //End identifyFormat()



void
showFormat(Format *f, char *tag)
{
    char tformat[128];
    char tend[128];

    switch(f->trajFormat)
    {
    case NOFORMAT:
        strcpy(tformat,"unknown");
        break;

    case XYZ:
        strcpy(tformat,"XYZ");
        break;

    case PDB:
        strcpy(tformat,"Protein Data Bank");
        break;

    case CHARMM:
        strcpy(tformat,"CHARMM");
        break;

    case Amber:
        strcpy(tformat,"Amber");
        break;

    default:
        strcpy(tformat,"invalid");
        break;
    }

    switch (f->endian)
    {
    case NOENDIAN:
        strcpy(tend,"unknown");
        break;

    case BIGENDIAN:
        strcpy(tend,"big");
        break;

    case LITTLEENDIAN:
        strcpy(tend,"little");
        break;

    default:
        strcpy(tend,"invalid");
        break;
    }

    fprintf(stdout,"== Trajectory %s format ====================\n",tag);
    fprintf(stdout,"%s, Int: %d bytes, Real: %d/%d bytes, %s endian, FORTRAN header: %d bytes\n",
            tformat,f->intSize,f->realSizeXYZ,f->realSizeOther,tend,f->fortranHeaderSize);

} //End showFormat()



Endian
myEndian()
{
    int num;
    char *nptr;
    int i;

    if (myendian==NOENDIAN)
    {
        nptr=(char *)&num;
        nptr[0]=1;
        for (i=1; i<sizeof(int); i++)
            nptr[i]=0;

        if (num == 1)
            myendian=LITTLEENDIAN;
        else
            myendian=BIGENDIAN;
    }
    return myendian;

} //End myEndian()



void
swapBytes(char *b, int size, int n)
{
    int j,k;
    char c,*cptr;

    cptr=b;
    for (j=0; j<n; j++)
    {
        for (k=0; k<size/2; k++)
        {
            c=cptr[k];
            cptr[k]=cptr[size-k-1];
            cptr[size-k-1]=c;
        }
        cptr+=size;
    }

} //End swapBytes()



void
convertIntegerFromBuffer(int *i, int n, char *buf, int size, Endian endian)
{
    int j,value;
    void *tptr;
    char *cptr;

    if (endian != myEndian())
        swapBytes(buf,size,n);

    cptr=buf;

    if (size==sizeof(int))
    {
        for (j=0; j<n; j++)
        {
            i[j]=*(int *)cptr;
            cptr+=size;
        }
    }
    else if (size>sizeof(int))
    {
        if (myEndian()==BIGENDIAN)
            cptr+=size-sizeof(int);
        for (j=0; j<n; j++)
        {
            i[j]=*(int *)cptr;
            cptr+=size;
        }
    }
    else
    {
        tptr=(void *)i;
        if (myEndian()==BIGENDIAN)
            tptr+=sizeof(int)-size;
        for (j=0; j<n; j++)
        {
            i[j]=0;
            memcpy(tptr,(void *)cptr,size);
            cptr+=size;
            tptr+=sizeof(int);
        }
    }
    
} //End convertIntegerFromBuffer()



void
convertIntegerToBuffer(int *i, int n, char *buf, int size, Endian endian)
{
    int j,k,value;
    void *tptr;
    char *cptr;

    cptr=buf;

    if (size==sizeof(int))
    {
        for (j=0; j<n; j++)
        {
            *(int *)cptr=i[j];
            cptr+=size;
        }
    }
    else if (size>sizeof(int))
    {
        if (myEndian()==BIGENDIAN)
        {
            for (j=0; j<n; j++)
            {
                for(k=0; k<size-sizeof(int); k++)
                    *cptr++=0;
                *(int *)cptr=i[j];
                cptr+=sizeof(int);
            }
        }
        else
        {
            for (j=0; j<n; j++)
            {
                *(int *)cptr=i[j];
                cptr+=sizeof(int);
                for(k=0; k<size-sizeof(int); k++)
                    *cptr++=0;
            }
        }
    }
    else
    {
        tptr=(void *)i;
        if (myEndian()==BIGENDIAN)
            tptr+=sizeof(int)-size;
        for (j=0; j<n; j++)
        {
            memcpy((void *)cptr,tptr,size);
            cptr+=size;
            tptr+=sizeof(int);
        }
    }

    if (endian != myEndian())
        swapBytes(buf,size,n);
} //End convertIntegerToBuffer



void
convertRealFromBuffer(double *f, int n, char *buf, int size, Endian endian)
{
    int j;
    char *cptr;
    float tf;

    if (endian != myEndian())
        swapBytes(buf,size,n);

    cptr=buf;

    if (sizeof(double) == size)
    {
        for (j=0; j<n; j++)
        {
            f[j]=*(double *)cptr;
            cptr+=size;
        }
    }
    else if (sizeof(float) == size)
    {
        for (j=0; j<n; j++)
        {
            f[j]=(double)(*(float *)cptr);
            cptr+=size;
        }
    }
    else
    {
        fprintf(stderr,"cannot handle real size %d\n",size);
        return;
    }

} //End convertRealFromBuffer()



void
convertRealToBuffer(double *f, int n, char *buf, int size, Endian endian)
{
    int j;
    char *cptr;
    float tf;

    cptr=buf;

    if (sizeof(double) == size)
    {
        for (j=0; j<n; j++)
        {
            *(double *)cptr=f[j];
            cptr+=size;
        }
    }
    else if (sizeof(float) == size)
    {
        for (j=0; j<n; j++)
        {
            (*(float *)cptr)=(float)f[j];
            cptr+=size;
        }
    }
    else
    {
        fprintf(stderr,"cannot handle real size %d\n",size);
        return;
    }

    if (endian != myEndian())
        swapBytes(buf,size,n);

} //End convertRealToBuffer()


void  dumpBuffer(char *tbuf, int size) {
  int i;
  fprintf(stderr,"buffer: ");
  for (i=0; i<size; i++) {
    fprintf(stderr,"%02x ",(unsigned char)tbuf[i]);
  }
  fprintf(stderr,"\n");
}

int
readFortran(DataStream *io, Format *f, char *buffer, int maxsize)
{
    char tbuf[8];
    int len;

    if (readData(io,tbuf,f->fortranHeaderSize)<=0)
    {
        fprintf(stderr, "debug: marker 1\n");
        return -1;
    }

    //dumpBuffer(tbuf,f->fortranHeaderSize);
    convertIntegerFromBuffer(&len,1,tbuf,f->fortranHeaderSize,f->endian);

    //  fprintf(stderr,"FORTRAN record len: %d, maxsize: %d\n",len,maxsize);

    if (len<=0 || len>maxsize)
    {
        fprintf(stderr, "debug: marker 2 [len %d][maxsize %d][fhsize %d][endian %d]\n", len, maxsize, f->fortranHeaderSize, f->endian);
        return -1;
    }

    if (readData(io,buffer,len)<=0)
    {
        fprintf(stderr, "debug: marker 3\n");
        return -1;
    }

    if (readData(io,tbuf,f->fortranHeaderSize)<=0)
    {
        fprintf(stderr, "debug: marker 4\n");
        return -1;
    }

    return len;

} //End readFortran()



int
writeFortran(DataStream *io, Format *f, char *buffer, int len)
{
    char tbuf[8];

    convertIntegerToBuffer(&len,1,tbuf,f->fortranHeaderSize,f->endian);
    if (writeData(io,tbuf,f->fortranHeaderSize)<=0) return -1;
    if (writeData(io,buffer,len)<=0) return -1;
    if (writeData(io,tbuf,f->fortranHeaderSize)<=0) return -1;

    return len;

} //End writeFortran()



Endian
getEndianFromString(char *str)
{
    if (str==0 || *str==0) return NOENDIAN;
    if (!strcasecmp(str,"big")) return BIGENDIAN;
    if (!strcasecmp(str,"little")) return LITTLEENDIAN;

    fprintf(stderr,"unknown endian string %s\n",str);
    return NOENDIAN;

} //End getEndianFromString()



TrajectoryFormat
getTrajFormatFromString(char *str)
{
    if (str==0 || *str==0) return NOFORMAT;
    if (!strcasecmp(str,"xyz")) return XYZ;
    if (!strcasecmp(str,"pdb")) return PDB;
    if (!strcasecmp(str,"charmm")) return CHARMM;
    if (!strcasecmp(str,"amber")) return Amber;

    fprintf(stderr,"unknown trajectory format string %s\n",str);
    return NOFORMAT;

} //End getTrajFormatFromString()



void
setFormatFromInput(Format *f, Format *d)
{
  if (d->trajFormat!=NOFORMAT) 
    f->trajFormat=d->trajFormat;
  if (d->endian!=NOENDIAN) 
    f->endian=d->endian;
  if (d->intSize>0 && d->intSize<=16) 
    f->intSize=d->intSize;
  if (d->realSizeXYZ>0 && d->realSizeXYZ<=16) 
    f->realSizeXYZ=d->realSizeXYZ;
  if (d->realSizeOther>0 && d->realSizeOther<=16) 
    f->realSizeOther=d->realSizeOther;
  if (d->fortranHeaderSize>0 && d->fortranHeaderSize<=16) 
    f->fortranHeaderSize=d->fortranHeaderSize;

} //End setFormatFromInput()



void
parseFormatString(char  *str, Format *fmt)
{
    char *cptr;
    int equal;
    char fstr[128];
    char tstr[128];
    int i,fi,ti;

    cptr=str;

    do
    {
        fi=0;
        ti=0;
	equal=0;
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

        if (!strcmp(fstr,"format"))
        {
            fmt->trajFormat=getTrajFormatFromString(tstr);
        }
        else if (!strcmp(fstr,"endian"))
        {
            fmt->endian=getEndianFromString(tstr);
        }
        else if (!strcmp(fstr,"intsize"))
        {
            fmt->intSize=atoi(tstr);
        }
        else if (!strcmp(fstr,"realsize") || !strcmp(fstr,"realsizexyz"))
        {
            fmt->realSizeXYZ=atoi(tstr);
        }
        else if (!strcmp(fstr,"realsizeother"))
        {
            fmt->realSizeOther=atoi(tstr);
        }
        else if (!strcmp(fstr,"fortranheadersize"))
        {
            fmt->fortranHeaderSize=atoi(tstr);
        }

        for (cptr+=i; *cptr==','; cptr++);
    }
    while(*cptr!=0);
} //End parseFormatString()

