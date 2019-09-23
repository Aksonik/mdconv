#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "datastream.h"

void
datastream_init(DataStream *io, int maxlen)
{
    io->maxlen = maxlen;
    io->bptr = 0;
    io->blen = 0;
    io->buffer = 0;
    io->fdes = -1;
} //End datastream_init()



DataStream *
openReadFile(char *fname)
{
    DataStream *io;
    io=(DataStream *)malloc(sizeof(DataStream));
    datastream_init(io, 10000);

    if (!strcmp(fname,"-")) {
      io->fdes=0;
    } else {
      if ((io->fdes=open(fname,O_RDONLY))<0)
        return (DataStream *)0;
    }

    io->buffer=(char *)malloc(io->maxlen);

    io->mode = getMode(io);

    return io;
} //End openReadFile()



DataStream *
openWriteFile(char *fname)
{
    DataStream *io;
    io=(DataStream *)malloc(sizeof(DataStream));
    datastream_init(io, 0);

    if ((io->fdes=open(fname,O_CREAT|O_WRONLY|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0)
        return (DataStream *)0;

    io->mode = getMode(io);

    return io;

} //End openWriteFile()



void
closeFile(DataStream *io)
{
    if (io->fdes>0)
        close(io->fdes);

    if (io->buffer!=0)
    {
        free(io->buffer);
        io->buffer=(char *)0;
    }

} //End closeFile()



int
readn(int fdes, char *buf, int n)
{
    int nreadtotal,nread,nleft;

    nreadtotal=0;
    nleft=n;
    do
    {
        nread=read(fdes,(void *)buf+nreadtotal,nleft);
        if (nread<0) return -1;
        if (nread==0) return nreadtotal;
        nreadtotal+=nread;
        nleft-=nread;
    }
    while (nreadtotal<n);
    return nreadtotal;

} //End readn()



int
writen(int fdes, char *buf, int n)
{
    int nwritetotal,nwritten,nleft;

    nwritetotal=0;
    nleft=n;
    do
    {
        nwritten=write(fdes,(void *)buf+nwritetotal,nleft);
        if (nwritten<=0) return -1;
        nwritetotal+=nwritten;
        nleft-=nwritten;
    }
    while (nwritetotal<n);

    return nwritetotal;

} //End writen()



int
fillBuffer(DataStream *io, int len)
{
    if (len>io->maxlen)
    {
        fprintf(stderr,"internal buffer size not large enough (requested: %d, available :%d)\n",
                len,io->maxlen);
        return -1;
    }
    if ((io->blen=readn(io->fdes,(char *)io->buffer,len))<0)
    {
        fprintf(stderr,"Error reading from file: ");
        fprintf(stderr,"%s\n",strerror(errno));
        return -1;
    }

    io->bptr=0;
    return 1;

} //End fillBuffer()



int
rewindBuffer(DataStream *io)
{
    io->bptr=0;
} //End rewindBuffer()



int
readData(DataStream *io, char *buf, int len)
{
    int ncopy;
    int dataleft;

    int reqlen=len;

    do
    {
        if (io->bptr>=io->blen)
            fillBuffer(io,io->maxlen);

        if (io->blen>0)
        {
            dataleft=io->blen-io->bptr;
            ncopy=(dataleft>reqlen)?reqlen:dataleft;
            memcpy(buf+len-reqlen,io->buffer+io->bptr,ncopy);
            io->bptr+=ncopy;
            reqlen-=ncopy;
        }
    }
    while(reqlen>0 && io->blen>0);

    //    fprintf(stderr,"read %d bytes\n",len-reqlen);

    return (len-reqlen);

} //End readData()

int 
readLine(DataStream *io, char *buf, int maxlen) 
{
    int maxncopy;
    int dataleft;
    int i;
    char *c;
    
    int len=maxlen-1;

    int reqlen=len;

    int foundcr=0;

    do
    {
        if (io->bptr>=io->blen)
            fillBuffer(io,io->maxlen);

        if (io->blen>0)
        {
            dataleft=io->blen-io->bptr;
            maxncopy=(dataleft>reqlen)?reqlen:dataleft;

	    c=(char *)(io->buffer+io->bptr);
	    for (i=0; i<maxncopy && c[i]!='\n'; i++);
	    
	    if (i<maxncopy) {
	      memcpy(buf+len-reqlen,(void *)c,i);
	      io->bptr+=i+1;
	      reqlen-=i;
	      foundcr=1;
	    } else {
	      memcpy(buf+len-reqlen,(void *)c,maxncopy);
	      io->bptr+=maxncopy;
	      reqlen-=maxncopy;
	    }
        }
    }
    while(reqlen>0 && io->blen>0 && !foundcr);

    buf[len-reqlen]=0;

    return (len-reqlen);
  
}

int
writeData(DataStream *io, char *buf, int len)
{
    return writen(io->fdes,(char *)buf,len);
} //End writeData()

int
datastream_setfdes(DataStream *io, int fdes)
{
    if( (io == NULL) || (fdes < 0) )
    {
        return -1;
    }

    io->fdes = fdes;
    io->mode = getMode(io);

    return fdes;

} //End datastream_setfdes()

off_t
getSize(DataStream *io) 
{
  struct stat info;

  if (!fstat(io->fdes,&info)) {
    if (S_ISREG(info.st_mode)) {
      return info.st_size;
    } else {
      return (off_t)0;
    }
  } else {
    return (off_t)-1;
  }
} //End getSize()


double getFixedNumber(char *buffer, int inx, int size) {
  char *cptr;
  char *dptr;
  char save;
  double value;

  cptr=buffer+inx*size;
  dptr=buffer+(inx+1)*size;
  
  save=*dptr;
  *dptr=0;
  value=atof(cptr);
  *dptr=save;

  return value;
} //End getFixedNumber()


int
fileSeek(DataStream *io, int len, int whence)
{
    off_t seek_amount;
    int retval = 0;

    if(S_ISREG(io->mode))
    {
        if( (io->bptr + len) >= (io->blen) )
        {
            seek_amount = io->bptr + len - io->blen;

            if(lseek(io->fdes, seek_amount, whence) >= 0)
            {
                if(fillBuffer(io,io->maxlen))
                {
                    retval = 1;
                }
            }
        }
        else
        {
            io->bptr += len;
        }
    }

    return retval;

} //End fileSeek()


mode_t
getMode(DataStream *io)
{
    struct stat info;

    if(fstat(io->fdes, &info) == 0)
    {
        return info.st_mode;
    }
    else
    {
        return (mode_t) 0;
    }

} //End getMode()

