#ifndef _DATASTREAM_H
#define _DATASTREAM_H

typedef struct {
  int  fdes;

  char *buffer;
  int maxlen;

  int blen;  /* length of data in buffer */
  int bptr;  /* offset of beginning of buffer */
  mode_t mode;  /* mode returned from stat of data source */
} DataStream;

DataStream *openReadFile(char *fname);
DataStream *openWriteFile(char *fname);

//DataStream *openSocket();

void closeFile(DataStream *io);

int fillBuffer(DataStream *io, int len);
int rewindBuffer(DataStream *io);
int readData(DataStream *io, char *buf, int len);
int readLine(DataStream *io, char *buf, int maxlen);

int writeData(DataStream *io, char *buf, int len);

void datastream_init(DataStream *io, int maxlen);
int datastream_setfdes(DataStream *io, int fdes);

int readn(int fdes, char *buf, int n);
int writen(int fdes, char *buf, int n);

off_t getSize(DataStream *io);

double getFixedNumber(char *buffer, int inx, int size);

int fileSeek(DataStream *io, int len, int whence);
mode_t getMode(DataStream *io);


#endif
