//
//
// Very Simple PFILE interface.
// 
//
//
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu



#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <values.h>
#include <math.h>
/*#include <ieeefp.h> */
#include <float.h>
#include <assert.h>
#include "error.h"
#include "spi.h"
#include "general.h"

#define PF_DEBUGLEVEL 0
#define PF_GENINDEX 1



// ======================================================================
// ============           SPI         ===================================
// ======================================================================

SPI::SPI(const char *const fname, bool swap)
{
  if (fname == NULL)
    error("SPI::SPI Can't open NULL file");
  
  inf = fopen(fname,"r");
  if (inf == NULL)
    error("SPI::SPI can't open file (%s)",fname);
  in_streamp = new InFtrLabStream_PFile(PF_DEBUGLEVEL, "", inf, 
				        PF_GENINDEX,swap);
  
  local_fname = copyToNewStr(fname);
  buf_size = 0;
}

SPI::~SPI()
{
  delete in_streamp;
  fclose(inf);
  delete [] local_fname;
}

size_t SPI::read_ftrs(const size_t pos,
		      float*& fb)
{
  const size_t n_frames = in_streamp->num_frames(pos);
  if (n_frames == SIZET_BAD) {
    error("SPI::read_ftrs: Couldn't find number of frames "
	  "at sentence %lu in input pfile.\n",(unsigned long) pos);
  }
  ftr_buf.growByNIfNeeded(2,n_frames*in_streamp->num_ftrs());
  const SegID seg_id = in_streamp->set_pos(pos, 0);
  if (seg_id == SEGID_BAD) {
    error("SPI::read_ftrs, Couldn't seek to start of sentence %lu "
	  "in input pfile.",(unsigned long) pos);
  }
  const size_t n_read =
    in_streamp->read_ftrs(n_frames, ftr_buf.ptr);
  if (n_read != n_frames) {
    error("SPI::read_ftrs At sentence %lu in input pfile, "
	  "only read %lu frames when should have read %lu.\n",
	  (unsigned long) pos, 
	  (unsigned long) n_read, (unsigned long) n_frames);
  }
  fb = ftr_buf.ptr;
  return n_read;
}


size_t SPI::read_labs(const size_t pos,
		      UInt32*& lb)
{
  const size_t n_frames = in_streamp->num_frames(pos);
  if (n_frames == SIZET_BAD) {
    error("SPI::read_labs: Couldn't find number of frames "
	  "at sentence %lu in input pfile.\n",(unsigned long) pos);
  }
  lab_buf.growByNIfNeeded(2,n_frames*in_streamp->num_labs());
  const SegID seg_id = in_streamp->set_pos(pos, 0);
  if (seg_id == SEGID_BAD) {
    error("SPI::read_labs, Couldn't seek to start of sentence %lu "
	  "in input pfile.",(unsigned long) pos);
  }
  const size_t n_read =
    in_streamp->read_labs(n_frames, lab_buf.ptr);
  if (n_read != n_frames) {
    error("SPI::read_labs At sentence %lu in input pfile, "
	  "only read %lu frames when should have read %lu.\n",
	  (unsigned long) pos, 
	  (unsigned long) n_read, (unsigned long) n_frames);
  }
  lb = lab_buf.ptr;
  return n_read;
}


size_t SPI::read_ftrslabs(const size_t pos,
			  float*& fb,
			  UInt32*& lb)
{
  const size_t n_frames = in_streamp->num_frames(pos);
  if (n_frames == SIZET_BAD) {
    error("SPI::read_labs: Couldn't find number of frames "
	  "at sentence %lu in input pfile.\n",(unsigned long) pos);
  }
  lab_buf.growByNIfNeeded(2,n_frames*in_streamp->num_labs());
  ftr_buf.growByNIfNeeded(2,n_frames*in_streamp->num_ftrs());
  const SegID seg_id = in_streamp->set_pos(pos, 0);
  if (seg_id == SEGID_BAD) {
    error("SPI::read_labs, Couldn't seek to start of sentence %lu "
	  "in input pfile.",(unsigned long) pos);
  }
  const size_t n_read =
    in_streamp->read_ftrslabs(n_frames, ftr_buf.ptr,lab_buf.ptr);
  if (n_read != n_frames) {
    error("SPI::read_labs At sentence %lu in input pfile, "
	  "only read %lu frames when should have read %lu.\n",
	  (unsigned long) pos, 
	  (unsigned long) n_read, (unsigned long) n_frames);
  }
  fb = ftr_buf.ptr;
  lb = lab_buf.ptr;
  return n_read;
}


// ======================================================================
// ============           SPI2         ==================================
// ======================================================================


SPI2::SPI2(const char *const fname,
	   const char *const fname2,
           bool swap1, bool swap2)
{
  if (fname == NULL)
    error("SPI::SPI Can't open NULL file");
  if (fname2 == NULL)
    error("SPI::SPI Can't open NULL second file");

  inf = fopen(fname,"r");
  if (inf == NULL)
    error("SPI::SPI can't open file (%s)",fname);
  in_streamp = new InFtrLabStream_PFile(PF_DEBUGLEVEL, "", inf,
				       PF_GENINDEX,swap1);


  inf2 = fopen(fname2,"r");
  if (inf2 == NULL)
    error("SPI::SPI can't open file (%s)",fname2);
  in_streamp2 = new InFtrLabStream_PFile(PF_DEBUGLEVEL, "", inf2, 
					    PF_GENINDEX,swap2);

  if (in_streamp->num_segs() != in_streamp2->num_segs()) {
    error("ERROR: the second pfile (name %s with %d segments) must have same number of segments as the first pfile has (name %s with %d segments)",
	  fname2,
	  in_streamp2->num_segs(),
	  fname,
	  in_streamp->num_segs()
	  );
  }

  if (in_streamp->num_ftrs() == 0)
    error("ERROR: pfile %s has zero features per frame.",fname);
  if (in_streamp2->num_ftrs() == 0)
    error("ERROR: pfile %s has zero features per frame.",fname2);

  for (int pos=0;pos<(int)in_streamp->num_segs();pos++) {
    if (in_streamp->num_frames(pos)
	!= 
	in_streamp2->num_frames(pos))
      {
	error("ERROR: the second pfile (%s) has %d frames at segment %d but the first one (%s) has %d frames at that segment.",
	      fname2,
	      in_streamp2->num_frames(pos),
	      pos,
	      fname,
	      in_streamp->num_frames(pos));
      }
  }
  
  local_fname = copyToNewStr(fname);
  local_fname2 = copyToNewStr(fname2);
  buf_size = 0;
}

SPI2::~SPI2()
{
  delete in_streamp;
  fclose(inf);
  delete [] local_fname;

  delete in_streamp2;
  fclose(inf2);
  delete [] local_fname2;

}


size_t SPI2::read_ftrs(const size_t pos,
		       float*& fb)
{

  // read in sentence from pfile 1
  const size_t n_frames1 = in_streamp->num_frames(pos);
  if (n_frames1 == SIZET_BAD) {
    error("SPI::read_ftrs: Couldn't find number of frames "
	  "at sentence %lu in input pfile %s.\n",(unsigned long) pos,local_fname);
  }
  ftr_buf1.growByNIfNeeded(2,n_frames1*in_streamp->num_ftrs());
  const SegID seg_id1 = in_streamp->set_pos(pos, 0);
  if (seg_id1 == SEGID_BAD) {
    error("SPI::read_ftrs, Couldn't seek to start of sentence %lu "
	  "in input pfile %s.",(unsigned long) pos,local_fname);
  }
  const size_t n_read1 =
    in_streamp->read_ftrs(n_frames1, ftr_buf1.ptr);
  if (n_read1 != n_frames1) {
    error("SPI::read_ftrs At sentence %lu in input pfile %s, "
	  "only read %lu frames when should have read %lu.\n",
	  (unsigned long) pos, 
	  local_fname,
	  (unsigned long) n_read1, (unsigned long) n_frames1);
  }

  // read in sentence from pfile 2 in to a separate buffer
  const size_t n_frames2 = in_streamp2->num_frames(pos);
  if (n_frames2 == SIZET_BAD) {
    error("SPI::read_ftrs: Couldn't find number of frames "
	  "at sentence %lu in input pfile %s.\n",(unsigned long) pos,local_fname2);
  }
  ftr_buf2.growByNIfNeeded(2,n_frames2*in_streamp2->num_ftrs());
  const SegID seg_id2 = in_streamp2->set_pos(pos, 0);
  if (seg_id2 == SEGID_BAD) {
    error("SPI::read_ftrs, Couldn't seek to start of sentence %lu "
	  "in input pfile %s.",(unsigned long) pos,local_fname2);
  }
  const size_t n_read2 =
    in_streamp2->read_ftrs(n_frames2, ftr_buf2.ptr);
  if (n_read2 != n_frames2) {
    error("SPI::read_ftrs At sentence %lu in input pfile %s, "
	  "only read %lu frames when should have read %lu.\n",
	  (unsigned long) pos, 
	  local_fname2,
	  (unsigned long) n_read2, (unsigned long) n_frames2);
  }


  // now we need to merge ftr_buf1 and ftr_buf2 into ftr_buf

  ftr_buf.growByNIfNeeded(2,n_frames2*
			  (in_streamp->num_ftrs() 
			   +
			   in_streamp2->num_ftrs()
			   ));

  float 
    *p1 = ftr_buf1.ptr,
    *p2 = ftr_buf2.ptr,
    *p  = ftr_buf.ptr;

  // Assumes that n_frames1 is always > 0
  assert ( n_frames1 > 0 );
  int feat = n_frames1; do {
    int i;
    // Note that this assumes that there will always
    // be greater than 0 features.
    i=in_streamp->num_ftrs();
    do {
      *p++ = *p1++;
    } while (--i != 0);
    i=in_streamp2->num_ftrs();
    do {
      *p++ = *p2++;
    } while (--i != 0);
  } while (--feat != 0);

  fb = ftr_buf.ptr;
  return n_read2;
}



// only read labels from first pfile.
size_t SPI2::read_labs(const size_t pos,
		      UInt32*& lb)
{
  const size_t n_frames = in_streamp->num_frames(pos);
  if (n_frames == SIZET_BAD) {
    error("SPI::read_labs: Couldn't find number of frames "
	  "at sentence %lu in input pfile.\n",(unsigned long) pos);
  }
  lab_buf.growByNIfNeeded(2,n_frames*in_streamp->num_labs());
  const SegID seg_id = in_streamp->set_pos(pos, 0);
  if (seg_id == SEGID_BAD) {
    error("SPI::read_labs, Couldn't seek to start of sentence %lu "
	  "in input pfile.",(unsigned long) pos);
  }
  const size_t n_read =
    in_streamp->read_labs(n_frames, lab_buf.ptr);
  if (n_read != n_frames) {
    error("SPI::read_labs At sentence %lu in input pfile, "
	  "only read %lu frames when should have read %lu.\n",
	  (unsigned long) pos, 
	  (unsigned long) n_read, (unsigned long) n_frames);
  }
  lb = lab_buf.ptr;
  return n_read;
}


size_t SPI2::read_ftrslabs(const size_t pos,
			  float*& fb,
			  UInt32*& lb)
{
  SPI2::read_ftrs(pos,fb);
  return SPI2::read_labs(pos,lb);
}





// ======================================================================
// ============           SPO         ===================================
// ======================================================================



SPO::SPO(const char *const fname,
	 const size_t n_ftrs,
	 const size_t n_labs,
         bool swap)
{

  if (!strcmp(fname,"-")) {
    ouf = stdout;
  } else if ((ouf = fopen(fname, "w")) == NULL) {
    error("SPO::SPO Couldn't open output file (%s) for writing.",fname);
  }
  out_streamp = new OutFtrLabStream_PFile(PF_DEBUGLEVEL,
					     "",
					     ouf,
					     n_ftrs,
					     n_labs,
					     PF_GENINDEX,swap);
  local_fname = copyToNewStr(fname);
  pos = 0;
}

void
SPO::write_ftrslabs(const size_t len,
		    float* ftr_buf,
		    UInt32* lab_buf)
{
  out_streamp->write_ftrslabs(len, ftr_buf, lab_buf);
  out_streamp->doneseg((SegID) pos++);
}


SPO::~SPO()
{
  delete out_streamp;
  fclose(ouf);
  delete local_fname;
}


#ifdef MAIN


main(int argc,char *argv[])
{
  if (argc < 3) {
    fprintf(stderr,"ERROR: Need two arguments, [inpfile] [outpfile];")
    return -1;
  }  

  SPI ipf(argv[1]);

  SPO opf(argv[2],ipf.n_ftrs(),ipf.n_labs());

  float* ftr_buf;
  UInt32* lab_buf;

  for (int i=0;(unsigned)i<ipf.n_segs();i++) {
    size_t len = ipf.read_ftrslabs((size_t)i,ftr_buf,lab_buf);
    opf.write_ftrslabs(len,ftr_buf,lab_buf);
  }
}

#endif
