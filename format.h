#ifndef _FORMAT_H
#define _FORMAT_H

#include "datastream.h"

typedef enum { NOFORMAT, Amber, CHARMM, XYZ, PDB } TrajectoryFormat;
typedef enum { NOENDIAN, BIGENDIAN, LITTLEENDIAN, MIXEDENDIAN } Endian;

typedef struct {
  TrajectoryFormat trajFormat;
  int realSizeXYZ;
  int realSizeOther;
  int intSize;
  Endian endian;
  int fortranHeaderSize;
} Format;

Endian myEndian();
Endian getEndianFromString(char *str);

void clearFormat(Format *fmt);

TrajectoryFormat getTrajFormatFromString(char *str);

void parseFormatString(char *str, Format *fmt);

int identifyFormat(DataStream *io, Format *fmt);
void showFormat(Format *fmt, char *tag);

void swapBytes(char *buf, int size, int n);

void convertIntegerFromBuffer(int *i, int n, char *buf, int size, Endian endian);
void convertRealFromBuffer(double *f, int n, char *buf, int size, Endian endian);
void convertIntegerToBuffer(int *i, int n, char *buf, int size, Endian endian);
void convertRealToBuffer(double *f, int n, char *buf, int size, Endian endian);

int readFortran(DataStream *io, Format *f, char *buffer, int maxsize);
int writeFortran(DataStream *io, Format *f, char *buffer, int len);

void setFormatFromInput(Format *f, Format *d);

#endif
