#ifndef _AMBER_H
#define _AMBER_H

#include "datastream.h"
#include "format.h"
#include "trajdata.h"

int testHeaderAmber(DataStream *io, Format *fmt);
int readHeaderAmber(DataStream *io, Format *fmt, TrajectoryData *t);
int readFrameAmber(DataStream *io, Format *fmt, TrajectoryData *t, Frame *f);
int writeHeaderAmber(DataStream *io, Format *fmt, TrajectoryData *t);
int writeFrameAmber(DataStream *io, Format *fmt, TrajectoryData *t, Frame *f);
long calcFrameSizeAmber(Format *fmt, TrajectoryData *td);
long calcDataSizeAmber(Format *fmt, TrajectoryData *td);

#endif
