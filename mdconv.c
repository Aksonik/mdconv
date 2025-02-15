/* usage:
   mdconv [options] DCDfile(s)
   options: -info
            -verbose 
            -out DCDfile
            -pdb PDBfile
            -merge DCDfile -offset value
            -translate file
            -frames  from:to[=from:to] -skip num
            -framelist file
            -atoms   from:to[=from:to]
            -atomlist file
            -inpdata atoms=value,frames=value,timestep=value,starttime=value,
                     startframe=value,periodic=[1|0]
            -from|to format=[xyz|pdb|charmm|amber],endian=[big|little],
	             intsize=value,realsize=value,realsizeother=value,
                     fortranheadersize=value
            -unwrap -box xsize ysize zsize
            -minseg  from:to[=from:to]
            -avgatoms file
            -avgframes  
            -lsqfit
            -removedrift
            -contacts A D SIGfile
*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "datastream.h"
#include "format.h"
#include "trajdata.h"
#include "charmm.h"
#include "amber.h"

#include "pdb.h"
#include "vector.h"

#define MAXTRAJFILES 4096

int verbose=0;
int withinSelection(int inx, int *from, int *to, int n, int skip);
int getRangeFromString(char *atoms, int *from, int *to, int maxsel);
void unwrapFrame(Frame *from, Frame *last, Frame *addwrap,double xbox, double ybox, double zbox); 
void minsegFrame(Frame *from, int *minsegfrom, int *minsegto, int nminsegsel, double xbox, double ybox, double zbox); 
void alignFrame(Frame *f, Frame *ref);
int readAveragingList(char *fname, int ***avgfrom, int ***avgto, int **avgnsel);
int readTranslationTable(char *fname, int **tab);

//void frotu_(double *r, double *u);

void usage() {
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"mdconv [options] DCDfile(s)\n");
  fprintf(stderr,"options: -info\n");
  fprintf(stderr,"         -verbose\n");
  fprintf(stderr,"         -out DCDfile\n");
  fprintf(stderr,"         -pdb PDBfile\n");
  fprintf(stderr,"         -merge DCDfile -offset num\n");
  fprintf(stderr,"         -translate file\n");
  fprintf(stderr,"         -frames  from:to[=from:to]  -skip num\n");
  fprintf(stderr,"         -framelist file\n");
  fprintf(stderr,"         -atoms   from:to[=from:to]\n");
  fprintf(stderr,"         -atomlist file\n");
  fprintf(stderr,"         -inpdata atoms=value,frames=value,timestep=value,starttime=value,\n");
  fprintf(stderr,"                  startframe=value,periodic=[1|0]\n");
  fprintf(stderr,"         -from|to format=[xyz|pdb|charmm|amber],endian=[big|little],\n");
  fprintf(stderr,"                  intsize=value,realsize=value,realsizeother=value,\n");
  fprintf(stderr,"                  fortranheadersize=value\n");
  fprintf(stderr,"         -unwrap -box xsize ysize zsize\n");
  fprintf(stderr,"         -minseg from:to[=from:to]\n");
  fprintf(stderr,"         -avgatoms file\n");
  fprintf(stderr,"         -avgframes\n");
  fprintf(stderr,"         -lsqfit\n");
  fprintf(stderr,"         -removedrift\n");
  fprintf(stderr,"         -contacts A D SIGfile\n"); // greg

  exit(0);  
}

int main(int argc, char **argv) {
  int i,j,n;
  int infoQ=0;

  DataStream *inpio,*tinpio,*outio,*mergio;
  Format inpfmt,tinpfmt,outfmt,minpfmt;
  Format fromfmt,tofmt;
  TrajectoryData inptrajdata,tinptrajdata,outtrajdata,minptrajdata;
  TrajectoryData giventrajdata;
  Frame *frame=0;
  Frame *last=0;
  Frame *addwrap=0;
  Frame *sframe=0;
  Frame *mframe=0;   // for merging
  Frame *rframe=0; // for reference (lsqfit)
  Frame *aframe=0; // for averaging

  char ofile[256];
  char *trajfile[MAXTRAJFILES];
  char pdbfile[256];
  char sigfile[256];	// greg
  char transfile[256];
  char mfile[256];
  char avgatomsfile[256];
  char framelistfile[256];
  char atomlistfile[256];
  int mergeoffset=0;

  int ntrajfile=0;
  int it,itt;
  int is,iss;
  int iset,iseg;

  char frames[8192];
  char atoms[8192];
  char minseg[8192];

  int nframes;
  int maxframe;
  int minframe;
  int iframe;
  int totframes;
  int mtotframes;
  int skipnum=0;
  int nwriteframe=0;

  int from[1000000];
  int to[1000000];
  int nsel=0;

  int **avgatomsfrom;
  int **avgatomsto;
  int *avgatomsnsel;
  int navgatoms=0;

  int **framesfrom;
  int **framesto;
  int *framesnsel;
  int nframeseltot=0;

  int **atomsfrom;
  int **atomsto;
  int *atomsnsel;
  int natomseltot=0;

  int minsegfrom[100];
  int minsegto[100];
  int nminsegsel=0;

  int avgframes=0;
  int lsqfit=0;
  int lsqfitfirst=1;

  int *atomlist;
  int natomsel=0;
  int nframesel=0;

  int numAtoms;
  int mnumAtoms;

  int *transtab;
  int ntrans;

  int unwrap,removedrift;
  double xboxsize, yboxsize, zboxsize;
  int wrapped;

// greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg

  int contacts;
  int coni,conj;
  double condx,condy,condz,cond;

  double A;
  double D;

  FILE *confile;
  FILE *pairsfile;

  char resname[10000][4];	// number of residues and residue name characters (+1)
  char segid[10000][6];

  int nsigres=0,nsigval=0;
  char sigres[1000][4];
  float sigval[1000];

// greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg

  *ofile=0;  
  *pdbfile=0;
  *sigfile=0; // greg
  *frames=0;
  *atoms=0;
  *transfile=0;
  *avgatomsfile=0;
  *framelistfile=0;
  *atomlistfile=0;
  *mfile=0;
  *minseg=0;
  transtab=0;

  nsel=0;
  natomsel=0;
  nframesel=0;

  ntrans=0;

  inpio=0;
  outio=0;
  mergio=0;
  tinpio=0;

  unwrap=0;
  xboxsize=-1;
  yboxsize=-1;
  zboxsize=-1;

  removedrift=0;

  contacts=0;	// greg

  clearFormat(&tofmt);
  clearFormat(&fromfmt);
  
  clearTrajData(&giventrajdata);

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-info")) {
      infoQ=1;
    } else if (!strcmp(argv[i],"-verbose")) {
      verbose=1;
    } else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"-help")) {
      usage();
    } else if (!strcmp(argv[i],"-out")) {
      strcpy(ofile,argv[++i]);
    } else if (!strcmp(argv[i],"-pdb")) {
      strcpy(pdbfile,argv[++i]);
    } else if (!strcmp(argv[i],"-avgatoms")) {
      strcpy(avgatomsfile,argv[++i]);
    } else if (!strcmp(argv[i],"-framelist")) {
      strcpy(framelistfile,argv[++i]);
    } else if (!strcmp(argv[i],"-atomlist")) {
      strcpy(atomlistfile,argv[++i]);
    } else if (!strcmp(argv[i],"-avgframes")) {
      avgframes=1;
    } else if (!strcmp(argv[i],"-lsqfit")) {
      fprintf (stderr,"lsqfit not implemented in this version\n");
      exit(1);
      lsqfit=1;
    } else if (!strcmp(argv[i],"-merge")) {
      strcpy(mfile,argv[++i]);
    } else if (!strcmp(argv[i],"-offset")) {
      mergeoffset=atoi(argv[++i]);
    } else if (!strcmp(argv[i],"-inpdata")) {
      parseTrajDataString(argv[++i],&giventrajdata);
    } else if (!strcmp(argv[i],"-frames")) {
      strncpy(frames,argv[++i],8191);
    } else if (!strcmp(argv[i],"-skip")) {
      skipnum=atoi(argv[++i]);
    } else if (!strcmp(argv[i],"-atoms")) {
      strncpy(atoms,argv[++i],8191);
    } else if (!strcmp(argv[i],"-translate")) {
      strcpy(transfile,argv[++i]);
    } else if (!strcmp(argv[i],"-to")) {
      parseFormatString(argv[++i],&tofmt);
    } else if (!strcmp(argv[i],"-from")) {
      parseFormatString(argv[++i],&fromfmt);
    } else if (!strcmp(argv[i],"-unwrap")) {
      unwrap=1;
    } else if (!strcmp(argv[i],"-minseg")) {
      strncpy(minseg,argv[++i],8191);
    } else if (!strcmp(argv[i],"-removedrift")) {
      removedrift=1;
    } else if (!strcmp(argv[i],"-box")) {
      xboxsize=atof(argv[++i]);
      yboxsize=atof(argv[++i]);
      zboxsize=atof(argv[++i]);

// greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg

    } else if (!strcmp(argv[i],"-contacts")) {
      contacts=1;

      A=atoi(argv[++i]);
      D=atoi(argv[++i]);

      confile=fopen("con.dat","w");
      pairsfile=fopen("pairs.dat","w");

      strcpy(sigfile,argv[++i]); // read sigmas from a file
      strcpy(pdbfile,argv[++i]);

      FILE *fptr;
      fptr=fopen(pdbfile,"r");

      PDBEntry *pdb;
      pdb=new PDBEntry[400000];

      int natom=0;

      while(!feof(fptr))
      {
       if (pdb[natom].read(fptr)>0) 
       {
        strcpy(resname[natom],pdb[natom].residueName());
        strcpy(segid[natom],pdb[natom].segmentID());
       }
       natom++;
      }
      fclose(fptr);

// read sigmas from a file:

      FILE *fpsig;

      fpsig=fopen(sigfile,"r");

      int i=0;

      while(!feof(fpsig))
      {
       i++;
//       printf("ndx: %i\n",i);  
       if(i==1)
       {
        if(fscanf(fpsig,"%s",&sigres[nsigres])==1)
        {
         printf("res: %s\n",sigres[nsigres]);
         nsigres++;
        }
       }
       else if(i==2)
       {
        if(fscanf(fpsig,"%f",&sigval[nsigval])==1){
        printf("val: %f\n",sigval[nsigval]);
        i=0;
        nsigval++;
        }
       }
      }
/*
      printf("res first: %s\n",sigres[0]);
      printf("res last %s\n",sigres[nsigres-1]);

      printf("val first: %f\n",sigval[0]);
      printf("val last: %f\n",sigval[nsigval-1]);
*/
      fclose(fpsig);

// greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg

    } else {
      trajfile[ntrajfile]=(char *)malloc(strlen(argv[i])+1);
      strcpy(trajfile[ntrajfile++],argv[i]);
    }
  }

  if (*minseg!=0) {
    nminsegsel=getRangeFromString(minseg,minsegfrom,minsegto,100);
  }

  if (ntrajfile<=0) {
    trajfile[ntrajfile]=(char *)malloc(2);
    strcpy(trajfile[ntrajfile++],"-");
  }

  if (*ofile==0 && !infoQ && !contacts) {		// greg
    fprintf(stderr,"please provide output file\n");
    exit(1);
  }

  if (*ofile!=0) {
    if ((outio=openWriteFile(ofile))==0) {
      fprintf(stderr,"cannot open output file %s for writing\n",ofile);
      exit(1);
    }
    if (verbose) 
      fprintf(stdout,"Trajectory output file: %s\n",ofile);
  }

  if (*mfile!=0) {
    if ((mergio=openReadFile(mfile))==0) {
      fprintf(stderr,"cannot open merge DCD file %s for reading\n",mfile);
      exit(1);
    }
    if (verbose)
      fprintf(stdout,"Merging coordinates from file %s with offset %d\n",mfile,mergeoffset);

    clearFormat(&minpfmt);
    identifyFormat(mergio,&minpfmt);
    clearTrajData(&minptrajdata);
    readTrajData(mergio,&minpfmt,&minptrajdata);

    if (minptrajdata.nAtoms==0) {
      fprintf(stderr,"Number of atoms in merge file is zero!\n");
      exit(1);
    }

    if (minptrajdata.frames==0) {
      fprintf(stderr,"Number of frames iin merge file s zero!\n");
      exit(1);
    }

    mnumAtoms=minptrajdata.nAtoms;
    mframe=allocateFrame(mnumAtoms);

    mtotframes=minptrajdata.frames;
  } 

  if (*transfile!=0) {
    if (!access(transfile,R_OK)) {
      ntrans=readTranslationTable(transfile,&transtab);
    } else {
      fprintf(stderr,"cannot read translation file %s\n",transfile);
      exit(1);
    }
  }

  if (*avgatomsfile!=0) {
    if (!access(avgatomsfile,R_OK)) {
      navgatoms=readAveragingList(avgatomsfile,&avgatomsfrom,&avgatomsto,&avgatomsnsel);
    } else {
      fprintf(stderr,"cannot read atom averaging input file %s\n",avgatomsfile);
      exit(1);
    }
  }

  if (*framelistfile!=0) {
    if (!access(framelistfile,R_OK)) {
      nframeseltot=readAveragingList(framelistfile,&framesfrom,&framesto,&framesnsel);
    } else {
      fprintf(stderr,"cannot read frame list file %s\n",framelistfile);
      exit(1);
    }
  }

  if (*atomlistfile!=0) {
    if (!access(atomlistfile,R_OK)) {
      natomseltot=readAveragingList(atomlistfile,&atomsfrom,&atomsto,&atomsnsel);
      if (verbose) {
        fprintf(stderr,"read %d indices from atom list file\n",natomseltot);
      }
    } else {
      fprintf(stderr,"cannot read atom list file %s\n",atomlistfile);
      exit(1);
    }
  }

  iframe=0;
  for (it=0; it<ntrajfile; it++) {
    if ((inpio=openReadFile(trajfile[it]))==0) {
      fprintf(stderr,"cannot open trajectory file %s for reading\n",
              trajfile[it]);
      exit(1);
    }

    if (verbose) 
      fprintf(stdout,"Reading trajectory from: %s\n",trajfile[it]);

    clearFormat(&inpfmt);
    identifyFormat(inpio,&inpfmt);
    setFormatFromInput(&inpfmt,&fromfmt);
    if (verbose) showFormat(&inpfmt,"input");

    // preset trajdata to allow readHeader to figure out missing information
    
    clearTrajData(&inptrajdata);
    setTrajDataFromInput(&inptrajdata,&giventrajdata);
    readTrajData(inpio,&inpfmt,&inptrajdata);
    //set trajdata from user input
    setTrajDataFromInput(&inptrajdata,&giventrajdata);
    if (infoQ) showTrajData(&inptrajdata);    

    if (inptrajdata.nAtoms==0) {
      fprintf(stderr,"Number of atoms is zero!\n");
      exit(1);
    }

    if (inptrajdata.frames==0) {
      fprintf(stderr,"Number of frames is zero!\n");
      exit(1);
    }

    if (it==0) { 
      numAtoms=inptrajdata.nAtoms;
      frame=allocateFrame(numAtoms);

      if (*atoms!=0 || natomseltot>0) {
        if (*atoms!=0) 
          nsel=getRangeFromString(atoms,from,to,1024);
        if (natomseltot>0) {
          for (is=0; is<natomseltot; is++) {
            for (iss=0; iss<atomsnsel[is]; iss++) {
              from[nsel]=atomsfrom[is][iss];
              to[nsel]=atomsto[is][iss];
              nsel++;
            }
          }
        }

        for (natomsel=0,i=0; i<nsel; i++) 
         for (j=from[i]; j<=to[i]; j++) 
	  if (j>=0 && j<=numAtoms) natomsel++;
     
        if (natomsel>0) {
         if ((atomlist=(int *)malloc(sizeof(int)*natomsel))==0) {
	  fprintf(stderr,"cannot allocate storage for atom list\n");
	  exit(1);
         }
         for (natomsel=0,i=0; i<nsel; i++) 
	  for (j=from[i]; j<=to[i]; j++) 
	    if (j>=0 && j<=numAtoms) atomlist[natomsel++]=j;
        } else {
          fprintf(stderr,"invalid atom selection\n");
          exit(1);
        }
        sframe=allocateFrame(natomsel);
      } else if (ntrans>0) {
	sframe=allocateFrame(numAtoms);
      } else if (navgatoms>0) {
        sframe=allocateFrame(navgatoms);
      } else if (mnumAtoms>0) {
        if (mnumAtoms+mergeoffset>numAtoms) {
          sframe=allocateFrame(mergeoffset+mnumAtoms);
        } else {
          sframe=allocateFrame(numAtoms);
        }
      } else {
        sframe=frame;
      }

      if (lsqfit) 
        rframe=allocateFrame(sframe->nAtoms);
      
      if (avgframes) {
        aframe=allocateFrame(sframe->nAtoms);
        aframe->angleA=0.0;
        aframe->angleB=0.0;
        aframe->angleC=0.0;
      }

      minframe=-1;
      maxframe=-1;
      if (*frames!=0 || nframeseltot>0) {
        if (*frames!=0) 
          nframesel=getRangeFromString(frames,from,to,1024);
        if (nframeseltot>0) 
          for (is=0; is<nframeseltot; is++) {
            for (iss=0; iss<framesnsel[is]; iss++) {
              from[nframesel]=framesfrom[is][iss];
              to[nframesel]=framesto[is][iss];
              nframesel++;
            }
          }
      
        for (i=0; i<nframesel; i++) {
          if (to[i]>maxframe) maxframe=to[i];
          if (from[i]<minframe || minframe<0) minframe=from[i];
        }
      }

      totframes=inptrajdata.frames;
      for (itt=1; itt<ntrajfile; itt++) {
        if ((tinpio=openReadFile(trajfile[itt]))==0) {
          fprintf(stderr,"cannot open trajectory file %s for reading\n",
                  trajfile[itt]);
          exit(1);
        }
        identifyFormat(tinpio,&tinpfmt);
	setFormatFromInput(&tinpfmt,&fromfmt);

	clearTrajData(&tinptrajdata);
	setTrajDataFromInput(&tinptrajdata,&giventrajdata);
        readTrajData(tinpio,&tinpfmt,&tinptrajdata);
	setTrajDataFromInput(&tinptrajdata,&giventrajdata);
        totframes+=tinptrajdata.frames;

        if (tinptrajdata.fixedAtoms>0) 
          deallocateFrame(tinptrajdata.refFrame);
        closeFile(tinpio);
      }
      if (verbose) 
	fprintf(stderr,"total number of frames: %d\n",totframes);

      if (minframe<=0) minframe=1;
      if (maxframe>totframes || maxframe<0) maxframe=totframes; 

      if (nframesel<=1) {
        nframes=(maxframe-minframe)+1;
	nframes=(int)((double)nframes/(double)(skipnum+1)-0.001)+1;
      } else {
        nframes=0;
        for (i=minframe; i<=maxframe; i++) {
          if (withinSelection(i,from,to,nframesel,skipnum))
             nframes++;
        }
      }

      if (verbose)
	fprintf(stderr,"writing %d frames\n",nframes);

      if (outio!=0) {
	outfmt.trajFormat=inpfmt.trajFormat;
	outfmt.endian=myEndian();
	outfmt.intSize=sizeof(int);
	outfmt.realSizeXYZ=4;
        outfmt.realSizeOther=8;
	outfmt.fortranHeaderSize=4;
        setFormatFromInput(&outfmt,&tofmt);
        if (verbose) showFormat(&outfmt,"output");

        memcpy((void *)&outtrajdata,(void *)&inptrajdata,
               sizeof(TrajectoryData));
        outtrajdata.nAtoms=sframe->nAtoms;
        outtrajdata.frames=(avgframes?1:nframes);
        if (skipnum>0) {
          outtrajdata.timeStep=inptrajdata.timeStep*(double)(skipnum+1);
          outtrajdata.startTime=outtrajdata.timeStep;
        }
        outtrajdata.startTime+=outtrajdata.timeStep*(double)(minframe-1);
        writeHeader(outio,&outfmt,&outtrajdata);
      }
    } else {
      if (inptrajdata.nAtoms!=numAtoms) {
        fprintf(stderr,"number of atoms changed (previous: %d, current: %d)\n",
                numAtoms,inptrajdata.nAtoms);
        exit(1);
      }
    }

    if (minframe<0 || minframe<=iframe+inptrajdata.frames) {
      for (i=0; 
           i<inptrajdata.frames && (i+iframe<maxframe || maxframe<0); i++) {
        clearFrame(frame);
        if (readFrame(inpio,&inpfmt,&inptrajdata,frame)<0) {
          fprintf(stderr,"error reading next frame\n");
          exit(1);
        }

        if ((nframesel==0 && ((i+iframe)%(skipnum+1)==0)) 
	    || withinSelection(i+iframe+1,from,to,nframesel,skipnum)) {
          if (natomsel>0) 
    	    copyFrameSelection(frame,sframe,atomlist,natomsel,transtab,ntrans);
	  else if (ntrans>0) 
	    copyFrame(frame,sframe,transtab,ntrans);
          else if (navgatoms>0)
            copyAvgFrame(frame,sframe,navgatoms,avgatomsnsel,avgatomsfrom,avgatomsto);
          else if (mnumAtoms>0) 
            copyFrame(frame,sframe,(int *)0,0);

	  if (unwrap) {
	    if (last==0) {
	      last=allocateFrame(sframe->nAtoms);
	      addwrap=allocateFrame(sframe->nAtoms);
	      memcpy((void *)last->x,(void *)sframe->x,sizeof(double)*sframe->nAtoms);
	      memcpy((void *)last->y,(void *)sframe->y,sizeof(double)*sframe->nAtoms);
	      memcpy((void *)last->z,(void *)sframe->z,sizeof(double)*sframe->nAtoms);
    
              for (iset=0; iset<sframe->nAtoms; addwrap->x[iset++]=0.0);
              for (iset=0; iset<sframe->nAtoms; addwrap->y[iset++]=0.0);
              for (iset=0; iset<sframe->nAtoms; addwrap->z[iset++]=0.0);
            } else 
              unwrapFrame(sframe,last,addwrap,xboxsize,yboxsize,zboxsize);

            last->boxX=sframe->boxX;
            last->boxY=sframe->boxY;
            last->boxZ=sframe->boxZ;
            last->angleA=sframe->angleA;
            last->angleB=sframe->angleB;
            last->angleC=sframe->angleC;
	  } 

// greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg

	  if(contacts)
          {

//           fprintf(stderr,"calculate contacts here!\n");
//           fprintf(stderr,"number of atoms: %i\n",inptrajdata.nAtoms);
//           fprintf(stderr,"number of frames: %i\n",inptrajdata.frames);

           last=allocateFrame(sframe->nAtoms);
           memcpy((void *)last->x,(void *)sframe->x,sizeof(double)*sframe->nAtoms);
           memcpy((void *)last->y,(void *)sframe->y,sizeof(double)*sframe->nAtoms);
           memcpy((void *)last->z,(void *)sframe->z,sizeof(double)*sframe->nAtoms);

//           fprintf(stderr,"coordinates of atom 0: %lf %lf %lf\n",last->x[0],last->y[0],last->z[0]);
////           fprintf(stderr,"%lf\n",sframe->x[0]);  // the same
////           fprintf(stderr,"%i\n",sframe->nAtoms); // the same 
           fprintf(stderr,"frame number: %i\n",i);
//           fprintf(stderr,"box size: %lf %lf %lf\n",sframe->boxX,sframe->boxY,sframe->boxZ);

           int cc=0;

           for(coni=0;coni<inptrajdata.nAtoms-1;coni++)
           {
            for(conj=coni+1;conj<inptrajdata.nAtoms;conj++)
            {
//             fprintf(stderr,"coordinates of atom: %i %lf %lf %lf\n",
//             coni,last->x[coni],last->y[coni],last->z[coni]);
//             fprintf(stderr,"coordinates of atom: %i %i\n",
//             coni,conj);

             condx=fabs(last->x[coni]-last->x[conj]);
             condy=fabs(last->y[coni]-last->y[conj]);
             condz=fabs(last->z[coni]-last->z[conj]);

//             fprintf(stderr,"%f %f\n",condx,fabs(condx));

             condx=condx-(int)(condx/sframe->boxX)*sframe->boxX;
             condy=condy-(int)(condy/sframe->boxY)*sframe->boxY;
             condz=condz-(int)(condz/sframe->boxZ)*sframe->boxZ;

             if(condx>0.5*sframe->boxX)
             {condx=sframe->boxX-condx;
             }
             if(condy>0.5*sframe->boxY)
             {condy=sframe->boxY-condy;
             }
             if(condz>0.5*sframe->boxZ)
             {condz=sframe->boxZ-condz;
             }

             cond=condx*condx+condy*condy+condz*condz;

             double rgA;
             double rgB;

// to read residues properly the input pdb file cannot have anything else than atoms !!!
/*
             if(!strcmp(resname[coni],"CGA")){rgA=18.1273;}
             else if (!strcmp(resname[coni],"CGB")){rgA=20.6692;}
             else if (!strcmp(resname[coni],"CGC")){rgA=23.3446;}
             else if (!strcmp(resname[coni],"CGP")){rgA=23.3446;}
             else if (!strcmp(resname[coni],"CGL")){rgA=26.9569;}
             else if (!strcmp(resname[coni],"CGI")){rgA=47.7310;}

             if(!strcmp(resname[conj],"CGA")){rgB=18.1273;}
             else if (!strcmp(resname[conj],"CGB")){rgB=20.6692;}
             else if (!strcmp(resname[conj],"CGC")){rgB=23.3446;}
             else if (!strcmp(resname[conj],"CGP")){rgB=23.3446;}
             else if (!strcmp(resname[conj],"CGL")){rgB=26.9569;}
             else if (!strcmp(resname[conj],"CGI")){rgB=47.7310;}
*/

             int ir;

             for (ir=0;ir<=nsigres-1;ir++) {
             if (!strcmp(resname[coni],sigres[ir])) {
              rgA=sigval[ir];
//              printf("%i %f %f\n",ir,sigval[ir],rgA);
              }
             }
             for (ir=0;ir<=nsigres-1;ir++) {
             if (!strcmp(resname[conj],sigres[ir])) {
              rgB=sigval[ir];
//              printf("%i %f %f\n",ir,sigval[ir],rgB);
              }
             }

             double rc;

             rc=(A*(rgA+rgB)*0.5+D)*(A*(rgA+rgB)*0.5+D);

//             fprintf(stderr,"squared cutoff: %i %i %i %s %s %f %f %f %f\n",
//             i,coni,conj,resname[coni],resname[conj],rgA,rgB,sqrt(rc),sqrt(cond));

//             fprintf(stderr,"coordinates of atom: %i %lf %lf %lf\n",
//             coni,last->x[coni],last->y[coni],last->z[coni]);
//             fprintf(stderr,"coordinates of atom: %i %lf %lf %lf\n",
//             coni,last->x[conj],last->y[conj],last->z[conj]);


//             fprintf(stderr,"resid i-2: %s\n",resname[coni]);
//             fprintf(stderr,"resid i-1: %s\n",resname[coni]);


//             fprintf(stderr,"resid i: %s\n",resname[coni]);
//             fprintf(stderr,"resid j: %s\n",resname[conj]);


//             fprintf(stderr,"resid j+1: %s\n",resname[conj]);
//             fprintf(stderr,"resid j+2: %s\n",resname[conj]);


             if(cond<=rc)
             {cc=cc+1;
              fprintf(pairsfile,"%i %i %i\n",i+1,coni+1,conj+1);
//              fprintf(pairsfile,"%i %i %i %s %s %12.6f\n",i+1,coni,conj,resname[coni],resname[conj],sqrt(cond));
//              fprintf(pairsfile,"%i %i %i %s %s %s %s %12.6f\n",i+1,coni,conj,resname[coni],resname[conj],segid[coni],segid[conj],sqrt(cond));
             }

            }
           }

//           fprintf(stderr,"tot. num. of con.:%i %i\n",i,cc);
           fprintf(confile,"%i %i\n",i+1,cc);

          }

// greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg greg

          if (nminsegsel>1) {
            minsegFrame(sframe,minsegfrom,minsegto,nminsegsel,xboxsize,yboxsize,zboxsize);
          }

          if (mergio!=0) {
            clearFrame(mframe);
            if (readFrame(mergio,&minpfmt,&minptrajdata,mframe)<0) {
              fprintf(stderr,"error reading next frame from merge file\n");
              exit(1);
            }
	    mergeFrame(sframe,mframe,mergeoffset);
          }

          if (lsqfit) {
            if (lsqfitfirst) {
              copyFrame(sframe,rframe,(int *)0,0);
              lsqfitfirst=0;
            } else {
              alignFrame(sframe,rframe);
            }
          }

          if (infoQ) 
            showFrame(sframe,i+iframe+1);

          if (avgframes) 
            addFrame(sframe,aframe);  

          if (removedrift)
            removedriftFrame(sframe);

          if (outio!=0 && !avgframes) {
	    if (verbose) 
	      fprintf(stderr,"writing frame %d --> %d\n",i+iframe+1,++nwriteframe);
            writeFrame(outio,&outfmt,&outtrajdata,sframe);
	  }
        }
      }
    }
    iframe+=inptrajdata.frames;
  }

  if(contacts)
  {
   fclose(confile); // greg
   fclose(pairsfile); // greg
  }

  // finish up

  if (avgframes) {
    scaleFrame(aframe,1.0/(double)nframes);
    if (verbose) 
      fprintf(stderr,"writing average frame\n");
    writeFrame(outio,&outfmt,&outtrajdata,aframe);
    deallocateFrame(aframe);
  }

  deallocateFrame(frame);
  if (unwrap && last!=0)
    deallocateFrame(last);
  if (unwrap && addwrap!=0)
    deallocateFrame(addwrap);

  if (natomsel>0 || ntrans>0 || navgatoms>0 || mnumAtoms>0) 
    deallocateFrame(sframe);

  if (navgatoms>0) {
    for (i=0; i<navgatoms; i++) {
      free(avgatomsfrom[i]);
      free(avgatomsto[i]);
    }
    free(avgatomsfrom);
    free(avgatomsto);
    free(avgatomsnsel);
  }

  if (nframeseltot>0) {
    for (i=0; i<nframeseltot; i++) {
      free(framesfrom[i]);
      free(framesto[i]);
    }
    free(framesfrom);
    free(framesto);
    free(framesnsel);
  }

  if (natomseltot>0) {
    for (i=0; i<natomseltot; i++) {
      free(atomsfrom[i]);
      free(atomsto[i]);
    }
    free(atomsfrom);
    free(atomsto);
    free(atomsnsel);
  }

  if (inptrajdata.fixedAtoms>0) 
    deallocateFrame(inptrajdata.refFrame);

  closeFile(inpio);

  if (outio!=0) 
    closeFile(outio); 

  if (mergio!=0) 
    closeFile(mergio);

  for (it=0; it<ntrajfile; it++)
    free(trajfile[it]);

  if (ntrans>0) 
    free(transtab);

  exit(0);
  return 0;
}

int getRangeFromString(char *atoms, int *from, int *to, int maxsel) {
  // syntax is from:to[=from:to] ...
  char *cptr;
  int i,fi,ti;
  char fstr[128];
  char tstr[128];
  int colon;
  int nsel=0;
 
  cptr=atoms;
  do { 
    colon=0;
    fi=0;
    ti=0;
    for (i=0; cptr[i]!=0 && cptr[i]!='='; i++) {
      if (cptr[i]==':') 
        colon=i;
      else {
        if (colon==0) 
          fstr[fi++]=cptr[i];
        else 
          tstr[ti++]=cptr[i];
      }
    }
    if (nsel>=maxsel) {
       fprintf(stderr,"maximum number of selections reached\n");
       exit(1);
    }
    fstr[fi]=0;
    tstr[ti]=0;

    from[nsel]=atoi(fstr);
    if (*tstr!=0) 
      to[nsel]=atoi(tstr);
    else 
      to[nsel]=from[nsel];
    nsel++;
     
    for (cptr+=i; *cptr=='='; cptr++);
  } while(*cptr!=0);
  return nsel;
}

int withinSelection(int inx, int *from, int *to, int n, int skip) {
  int i;
  for (i=0; i<n; i++) 
    if (inx>=from[i] && inx<=to[i] && ((inx-from[i])%(skip+1))==0) return 1;
  return 0;
}

int readAveragingList(char *fname, int ***avgfrom, int ***avgto, int **avgnsel) {
  FILE *fptr;
  char **tempbuf;
  int n,i;
  int maxlines=2000000;

  if ((fptr=fopen(fname,"r"))==0) {
    fprintf(stderr,"cannot open file %s for reading\n",fname);
    exit(1);
  }

  n=0;
  if ((tempbuf=(char **)malloc(maxlines*sizeof(char *)))==0) {
    fprintf(stderr,"cannot allocate memory\n");
    exit(1);
  }
 
  for (i=0; i<maxlines; i++) {
    if ((tempbuf[i]=(char *)malloc(1024*sizeof(char)))==0) {
      fprintf(stderr,"cannot allocate memory for tempbuf[i]\n");
      exit(1);
    }
  }

  while(!feof(fptr)) {
    if (fgets(tempbuf[n],1024,fptr)) n++;
  }

  fclose(fptr);

  if (n>0) {
    if ((*avgfrom=(int **)malloc(sizeof(int *)*n))==0) {
      fprintf(stderr,"cannot allocate storage for avgfrom\n");
      exit(1);
    }
    if ((*avgto=(int **)malloc(sizeof(int *)*n))==0) {
      fprintf(stderr,"cannot allocate storage for avgto\n");
      exit(1);
    }
    if ((*avgnsel=(int *)malloc(sizeof(int)*n))==0) {
      fprintf(stderr,"cannot allocate storage for avgnsel\n");
      exit(1);
    }
    for (i=0; i<n; i++) {
      if (((*avgfrom)[i]=(int *)malloc(sizeof(int)*1024))==0) {
         fprintf(stderr,"cannot allocate storage for avgfrom[i]\n");
         exit(1);
      }
      if (((*avgto)[i]=(int *)malloc(sizeof(int)*1024))==0) {
         fprintf(stderr,"cannot allocate storage for avgto[i]\n");
         exit(1);
      }
      (*avgnsel)[i]=getRangeFromString(tempbuf[i],(*avgfrom)[i],(*avgto)[i],1024);
    }
  }

  for (i=0; i<maxlines; i++) {
    free(tempbuf[i]);
  }
  free(tempbuf);

  return n;
}

int readTranslationTable(char *fname, int **tab) {
  FILE *fptr;
  int *tempbuf;
  int n,i;
  int maxlines=2000000;
  
  if ((fptr=fopen(fname,"r"))==0) {
    fprintf(stderr,"cannot open file %s for reading\n",fname);
    exit(1);
  }

  n=0;
  if ((tempbuf=(int *)malloc(maxlines*sizeof(int)))==0) {
    fprintf(stderr,"cannot allocate memory\n");
    exit(1);
  }

  while(!feof(fptr)) {
    if (fscanf(fptr,"%d",&tempbuf[n])) n++;
  }

  fclose(fptr);

  if (n>0) {
    if ((*tab=(int *)malloc(n*sizeof(int)))==0) {
      fprintf(stderr,"cannot allocate memory\n");
      exit(1);
    }
    
    for (i=0; i<n; i++) {
      (*tab)[i]=tempbuf[i];
    }
  } else {
    *tab=0;
  }

  free(tempbuf);
  
  return n;
}  

void
unwrapFrame(Frame *from, Frame *last, Frame *addwrap,double xbox, double ybox, double zbox) {
  int i;

  if (xbox<0)
    xbox=from->boxX;

  if (ybox<0)
    ybox=from->boxY;

  if (zbox<0)
    zbox=from->boxZ;

  for (i=0; i<from->nAtoms; i++) {
    if (last->x[i]-from->x[i]>0.5*xbox) {
      addwrap->x[i]+=1.0;
    } else if (from->x[i]-last->x[i]>0.5*xbox) {
      addwrap->x[i]-=1.0;
    }
    last->x[i]=from->x[i];
    from->x[i]+=addwrap->x[i]*xbox;

    if (last->y[i]-from->y[i]>0.5*ybox) {
      addwrap->y[i]+=1.0;
    } else if (from->y[i]-last->y[i]>0.5*ybox) {
      addwrap->y[i]-=1.0;
    }
    last->y[i]=from->y[i];
    from->y[i]+=addwrap->y[i]*ybox;

    if (last->z[i]-from->z[i]>0.5*zbox) {
      addwrap->z[i]+=1.0;
    } else if (from->z[i]-last->z[i]>0.5*zbox) {
      addwrap->z[i]-=1.0;
    }
    last->z[i]=from->z[i];
    from->z[i]+=addwrap->z[i]*zbox;
  }
}

void 
minsegFrame(Frame *from, int *minsegfrom, int *minsegto, int nminsegsel, double xbox, double ybox, double zbox) {
  int i,iseg,iset;
  double avgx[100];
  double avgy[100];
  double avgz[100];
  int nset; 

  if (xbox<0)
    xbox=from->boxX;

  if (ybox<0)
    ybox=from->boxY;

  if (zbox<0)
    zbox=from->boxZ;

//  fprintf(stderr,"box size %lf %lf %lf\n",xbox,ybox,zbox);

  for (iseg=0; iseg<nminsegsel; iseg++) {
    avgx[iseg]=0.0;
    avgy[iseg]=0.0;
    avgz[iseg]=0.0;
    nset=0;
    for (iset=minsegfrom[iseg]; iset<=minsegto[iseg]; iset++) {
      avgx[iseg]+=from->x[iset-1];
      avgy[iseg]+=from->y[iset-1];
      avgz[iseg]+=from->z[iset-1];
      nset++;
    }
    avgx[iseg]/=(double)nset;
    avgy[iseg]/=(double)nset;
    avgz[iseg]/=(double)nset;
//    fprintf(stderr,"avg %d: %lf %lf %lf %d %d %d\n",iseg,avgx[iseg],avgy[iseg],avgz[iseg],nset,minsegfrom[iseg],minsegto[iseg]);
  }

  for (iseg=1; iseg<nminsegsel; iseg++) {
    double dx=avgx[iseg]-avgx[0];
    double dy=avgy[iseg]-avgy[0];
    double dz=avgz[iseg]-avgz[0];
    
    int idelx=(dx>0)?(int)(dx/xbox+0.5):(int)(dx/xbox-0.5);
    int idely=(dy>0)?(int)(dy/ybox+0.5):(int)(dy/ybox-0.5);
    int idelz=(dz>0)?(int)(dz/zbox+0.5):(int)(dz/zbox-0.5);

    if (idelx!=0 || idely!=0 || idelz!=0) {
      for (iset=minsegfrom[iseg]; iset<=minsegto[iseg]; iset++) {
        from->x[iset-1]-=(double)idelx*xbox;
        from->y[iset-1]-=(double)idely*ybox;
        from->z[iset-1]-=(double)idelz*zbox;
      }
    }
  }
}

void
alignFrame(Frame *f, Frame *ref) {
  int i;
  
  double crx,cry,crz;
  double ccx,ccy,ccz;

  for (i=0; i<f->nAtoms; i++) {
    crx+=ref->x[i];
    cry+=ref->y[i];
    crz+=ref->z[i];
 
    ccx+=f->x[i];
    ccy+=f->y[i];
    ccz+=f->z[i];
  }
  crx/=(double)f->nAtoms;
  cry/=(double)f->nAtoms;
  crz/=(double)f->nAtoms;

  ccx/=(double)f->nAtoms;
  ccy/=(double)f->nAtoms;
  ccz/=(double)f->nAtoms;

  double r[9];
  double u[9];

  for (i=0; i<9; i++) 
    r[i]=0.0;

  double rx,ry,rz;
  double cx,cy,cz;
  for (i=0; i<f->nAtoms; i++) {
    rx=ref->x[i]-crx;
    ry=ref->y[i]-cry;
    rz=ref->z[i]-crz;
  
    f->x[i]-=ccx;
    f->y[i]-=ccy;
    f->z[i]-=ccz;

    cx=f->x[i];
    cy=f->y[i];
    cz=f->z[i];

    r[0]+=cx*rx;
    r[1]+=cx*ry;
    r[2]+=cx*rz; 

    r[3]+=cy*rx;
    r[4]+=cy*ry;
    r[5]+=cy*rz; 

    r[6]+=cz*rx;
    r[7]+=cz*ry;
    r[8]+=cz*rz; 
  }
  //frotu_(r,u);
  
  for (i=0; i<f->nAtoms; i++) {
    cx=u[0]*f->x[i]+u[3]*f->y[i]+u[6]*f->z[i]+crx;  
    cy=u[1]*f->x[i]+u[4]*f->y[i]+u[7]*f->z[i]+cry;  
    cz=u[2]*f->x[i]+u[5]*f->y[i]+u[8]*f->z[i]+crz;  
    f->x[i]=cx;
    f->y[i]=cy;
    f->z[i]=cz;
  } 
}

