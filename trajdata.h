#ifndef _TRAJDATA_H
#define _TRAJDATA_H

#include "format.h"

typedef struct {
  int nAtoms;

  double *x;
  double *y;
  double *z; 

  double boxX;
  double boxY;
  double boxZ;
  
  double angleA;
  double angleB;
  double angleC;
} Frame;

typedef struct {
  int frames;
  int nAtoms;
  double startTime;
  int startFrame;
  double timeStep;
  int periodicBoundariesQ;
  int fixedAtoms;
  int velocitiesQ;
  int fourdimensQ;
  int fluctChargesQ;
  int degreesOfFreedom;
  int versionnum;
  char version[256];
  int titleLines;
  int titleLength;
  char **title;
  int *fixedIndex;
  Frame *refFrame;
} TrajectoryData;

void clearTrajData(TrajectoryData *t);
void showTrajData(TrajectoryData *t);

int readTrajData(DataStream *io, Format *fmt, TrajectoryData *t);

int writeHeader(DataStream *io, Format *fmt, TrajectoryData *t);

Frame *allocateFrame(int n);
void clearFrame(Frame  *f);
void deallocateFrame(Frame *f);

void copyFrameSelection(Frame *from, Frame *to, int *inx, int n, int *tinx, int tn);
void copyFrame(Frame *from, Frame *to, int *tinx, int tn);
void addFrame(Frame *from, Frame *to);
void scaleFrame(Frame *f, double factor);
void copyAvgFrame(Frame *from, Frame *to, int navg, int *nsel, int **afrom, int **ato);
void mergeFrame(Frame *to, Frame *merge, int offset);
void removedriftFrame(Frame *f);

int readFrame(DataStream *io, Format *fmt, TrajectoryData *td, Frame *f);
int writeFrame(DataStream *io, Format *fmt, TrajectoryData *td, Frame *f);

void showFrame(Frame *f, int inx);

void setTrajDataFromInput(TrajectoryData *f, TrajectoryData *d);

void parseTrajDataString(char *str, TrajectoryData *d);

long calcFrameSize(Format *fmt, TrajectoryData *td);
long calcDataSize(Format *fmt, TrajectoryData *td);

#endif
