#ifndef _CHARMM_H
#define _CHARMM_H

#include "datastream.h"
#include "format.h"
#include "trajdata.h"

int testHeaderCHARMM(DataStream *io, Format *fmt);
int readHeaderCHARMM(DataStream *io, Format *fmt, TrajectoryData *t);
int readFrameCHARMM(DataStream *io, Format *fmt, TrajectoryData *t, Frame *f);
int writeHeaderCHARMM(DataStream *io, Format *fmt, TrajectoryData *t);
int writeFrameCHARMM(DataStream *io, Format *fmt, TrajectoryData *t, Frame *f);
long calcFrameSizeCHARMM(Format *fmt, TrajectoryData *td);
long calcDataSizeCHARMM(Format *fmt, TrajectoryData *td);

#endif
