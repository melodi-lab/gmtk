
//static char* rcsid = "$Id$";
// revised pfile code (withouth intvec,fltvec or quicknet) - KK

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "pfile.h"

// This is the verision string that appears in the PFile header

static char* pfile_version0_string = 
    "-pfile_header version 0 size 32768";


// This routine is used to get one unsigned integer argument from a pfile
// header stored in a buffer. It returns 0 on success, else -1.


intv_int32_t
swapb_i32_i32(intv_int32_t val)
{
  intv_uint32_t uval;
  intv_uint32_t res;
  intv_int32_t b0, b1, b2, b3;

  uval = (intv_uint32_t) val;
  b0 = uval >> 24;
  b1 = (uval >> 8) & 0x0000ff00;
  b2 = (uval << 8) & 0x00ff0000;
  b3 = uval << 24;

  res = b0 | b1 | b2 | b3;
  return (intv_int32_t) res;
}


// TODO: this should be changed to i16 at some point.
short
swapb_short_short(short sval) {
 
  short res;
  short s0,s1;

  usval = (unsigned short) sval;  

  s0 = usval >> 8;
  s1 = (usval << 8) & 0x0000ff00;
  res = s0 | s1 ;
  
  return res;
}
	

intv_int32_t
copy_i32_i32(intv_int32_t from) {
  return from;
}

void
copy_i32_vi32(size_t len, intv_int32_t from, intv_int32_t* to)
{
  size_t i;

  for (i=0; i<len; i++)
    *to++ = from;
}

void
swapb_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;
for (i=0; i<len; i++)
    *to++ = swapb_i32_i32(*from++);
}

void
copy_vi32_vi32(size_t len, const intv_int32_t* from, intv_int32_t* to)
{
  size_t i;

  for (i=0; i<len; i++)
  *to++ = *from++;
}

void
copy_f_vf(size_t len, float from, float* to)
{
  size_t i;

  for (i=0; i<len; i++)
    *to++ = from;
}


static int
get_uint(const char* hdr, const char* argname, unsigned int* val)
{
    const char* p;		// Pointer to argument
    int count = 0;		// Number of characters scanned

    // Find argument in header
    p = strstr(hdr, argname);
    if (p==NULL)
	return -1;
    // Go past argument name
    p += strlen(argname);
    // Get value from stfing
    sscanf(p, " %u%n", val, &count);

    // We expect to pass one space, so need >1 characters for success.
    if (count > 1)
	return 0;
    else
	return -1;
}

// This routine is used to find string-valued data in the PFile header.
// Returns pointer to string or null.
static char*
get_str(const char* hdr, const char* argname)
{
    char* p;		// Pointer to argument.

    // Find argument in header.
    p = strstr(hdr, argname);
    if (p==NULL)
	return NULL;
    // Go past argument name.
    p += strlen(argname);
    // Pass over spaces.
    while (isspace(*p))
	p++;
    return p;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// InFtrLabStream_PFile
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

InFtrLabStream_PFile::InFtrLabStream_PFile(int a_debug,
 					   const char* a_filename,
					   FILE* a_file,
					   int a_indexed,
                                           short swap)
  : file(a_file),
    filename(a_filename),
    indexed(a_indexed),
    buffer(NULL),
    sentind(NULL),
    bswap(swap)
{
    if (fseek(file, 0, SEEK_SET))
    {
	error("Failed to seek to start of pfile '%s' header - "
		 "cannot read PFiles from streams.",
	      filename);
    }

    read_header();

    // Allocate frame buffer.
    buffer = new PFile_Val[num_cols];

    // Remember the width, in bytes, of one column
    bytes_in_row = num_cols * sizeof(PFile_Val);

    if (indexed)
    {
	// Allocate space for sentence index
	sentind = new UInt32[total_sents + 1];
	if (sentind_offset != 0)
	    build_index_from_sentind_sect();
	else
	    build_index_from_data_sect();
    }
    // Move to the start of the PFile
    rewind();
}

InFtrLabStream_PFile::~InFtrLabStream_PFile()
{
    if (indexed)
	delete[] sentind;
    delete[] buffer;
}


void
InFtrLabStream_PFile::read_header()
{
    char* header;
    unsigned int ndim;		// Number of dimensions of data section
    unsigned int size;		// Size of data section
    unsigned int rows;		// Number of rows in section
    unsigned int cols;		// Number of columns in section
    long offset;		// offset in data section
    char* p;			// Temporary pointer
    int ec;			// Error code

  // Allocate space for header - do not use stack, as putting big things
  // on the stack is dangerous on SPERT.

   header = new char[PFILE_HEADER_SIZE]; // Store header here

  // Read in pfile header

    if (fread(header, PFILE_HEADER_SIZE, 1, file)!=1)
	error("Failed to read pfile '%s' header.",
		 filename);

  // Check pfile header

    if (strstr(header, pfile_version0_string)==NULL)
	error("Bad PFile header version in '%s'.",
	      filename);

    p = strstr(header, "-data");
    if (p==NULL)
	error("Cannot find pfile -data parameter in header of"
		 " '%s'.", filename);

    sscanf(p, "-data size %u offset %li ndim %u nrow %u ncol %u",
	   &size, &offset, &ndim, &rows, &cols);
    if (offset!=0 || ndim!=2 || (rows*cols)!=size)
	error("Bad or unrecognized pfile header -data args in"
		 " '%s'.", filename);

 // Find feature, label and target details and check okay

    ec = get_uint(header, "-first_feature_column",
			   &first_ftr_col);
    ec |= get_uint(header, "-num_features",
			    &num_ftr_cols);

    first_lab_col = 0;      // Initialize to zero.
    if (get_uint(header, "-num_labels", &num_lab_cols))
        num_lab_cols = 0;       // Set to zero if not found.
    if (num_lab_cols > 0)
    {
        ec |= get_uint(header, "-first_label_column", &first_lab_col);
    }

// A few other bits of information
    ec |= get_uint(header, "-num_frames", &total_frames);
    if (total_frames!=rows)
	error("Inconsistent number of frame ins header for"
		 " pfile '%s'.", filename);
    ec |= get_uint(header, "-num_sentences", &total_sents);
    if (ec)
	error("Problems reading pfile '%s' header.",
		 filename);


 // Check the "format" field.
 // The format value is a string containing 'f's, 'd's or prefixed repeat
 // counts.
 // e.g. ddffffffd or 2d24f

	long i;		// Local counter.

 // First unpack the string.

	char* unpacked_format = new char[cols];
	memset(unpacked_format, '\0', cols); // Set to null.
	char* hdr_format = get_str(header, "-format");
	if (hdr_format==NULL)
	{
	    error("Could not find '-format' field in PFile '%d'.",
		      filename);
	}
	char* upk_ptr = unpacked_format;
	unsigned upk_cnt = 0;
	char c = *hdr_format++;
	while (c=='f' || c=='d' || isdigit(c))
	{
	    if (c=='f' || c=='d')
	    {
		if (upk_cnt>cols)
		{
		    error("Format too long in PFile '%s' header.",
			      filename);
		}
		*upk_ptr++ = c;
		upk_cnt++;
		c = *hdr_format++;
	    }
	    else
	    {
		long rpt;	// No of times to repeat formatting character.

		rpt = strtol((hdr_format-1), &hdr_format, 10);
		upk_cnt += rpt;
		if (upk_cnt > cols)
		{
		    error("Format too long in PFile '%s' header.",
			      filename);
		}
		c = *hdr_format++;
		if (c!='f' && c!='d')
		{
		    error("Format field corrupted in PFile '%s' header.",
			      filename);
		}
		for (i=0; i<rpt; i++)
		{
		    *upk_ptr++ = c;
		}
		c = *hdr_format++;
	    }
	}

// Now check the string.
	// Two 'd's first for sentence and column.
	if (unpacked_format[0]!='d' || unpacked_format[1]!='d')
	{
	    error("Format field corrupted in PFile '%s' header.",
		      filename);
	}
	size_t col;
	for (col=first_ftr_col; col<(first_ftr_col+num_ftr_cols); col++)
	{
	    if (unpacked_format[col]!='f')
	    {
		error("Format field corrupted in PFile '%s' header - "
			  "not all feature columns have format 'f'.",
			  filename);
	    }
	}
	for (col=first_lab_col; col<(first_lab_col+num_lab_cols); col++)
	{
	    if (unpacked_format[col]!='d')
	    {
		error("Format field corrupted in PFile '%s' header - "
			  "not all label columns have format 'd'.",
			  filename);
	    }
	}
	delete[] unpacked_format;
    

// Put some stuff we stored locally in the main structure
    data_offset = (offset*sizeof(PFile_Val) + PFILE_HEADER_SIZE);
    num_cols = cols;

// Get the sentence index information
    sentind_offset = 0;
    
    p = strstr(header, "-sent_table_data");
    if (p!=NULL)
    {
	sscanf(p, "-sent_table_data size %u offset %li ndim %u",
	       &size, &offset, &ndim);
	if (size!=total_sents+1 || ndim!=1)
	{
	    error("Bad or unrecognized header"
		     " -sent_table_data args in PFile '%s'.", 
		     filename);
	}
	sentind_offset = (offset*sizeof(PFile_Val) + PFILE_HEADER_SIZE);
    }

// Some simple consistency checks on header information
    if (first_ftr_col >= num_cols
	|| (first_ftr_col + num_ftr_cols) > num_cols
	|| first_lab_col >= num_cols
	|| (first_lab_col + num_lab_cols) > num_cols
	|| ( (first_lab_col>=first_ftr_col) // Check for overlapped labs & ftr
	     && (first_lab_col<(first_ftr_col+num_ftr_cols)) )
	|| total_sents > total_frames
	)
    {
	error("Inconsistent pfile header values in '%s' "
		 "- probably corrupted PFile.", filename);
    }
    delete[] header;
}

// Build a sentence index using the sentence index section in the PFile

void
InFtrLabStream_PFile::build_index_from_sentind_sect()
{
    if (fseek(file, sentind_offset, SEEK_SET)!=0)
    {
        error("Failed to move to start of sentence index data, "
		  "sentind_offset=%li - probably a corrupted PFile.",
		  sentind_offset);
    }
    long size = sizeof(PFile_Val) * (total_sents + 1);
    if (fread((char *) sentind, size, 1, file)!=1)
    {
	error("Failed to read sent_table_data section in '%s'"
		 " - probably a corrupted PFile.", filename);
    }

    // The index is in big-endian format on file - convert it to the host
    // endianness


    if (bswap)  
      swapb_vi32_vi32(total_sents+1, (intv_int32_t*) sentind,
	  	    (intv_int32_t*) sentind);
      

    // Check that the index of the frame after the last one is the same
    // as the number of frames in the PFile    

    if (sentind[total_sents] != total_frames)
    {
	error("Last sentence index (%lu) does not correspond"
		 " with number of frames (%i) in PFile '%s' - probably a"
		 " corrupted PFile.", (unsigned long) sentind[total_sents],
		 total_frames, filename);
    }
}

// Build a sentence index using the data section in the PFile
// This entails scanning the whole PFile

void
InFtrLabStream_PFile::build_index_from_data_sect()
{
    long last_sentno = -1;	// Number of the last sentence.
    long next_frameno = -1;	// Number of the next frame within sentence.
    long abs_frameno = 0;	// Frame number from beginning of data.

    if (fseek(file, data_offset, SEEK_SET)!=0)
    {
        error("Failed to move to start of data in '%s',"
		  " data_offset=%li - probably a corrupted PFile.",
		  filename, data_offset);
    }

    do
    {
	if (fread((char *) buffer, sizeof(PFile_Val)*num_cols,
		  1, file) != 1)
	{
	    if (feof(file))
	    {
		last_sentno++;
		break;
	    }
	    else
	    {
		error("Failed to read pfile record from '%s',"
			 "last_sentno=%li next_frameno=%li abs_frameno=%li "
			 "filepos=%li - probably a corrupted PFile.",
			 filename,
			 last_sentno, next_frameno, abs_frameno, ftell(file));
	    }
	}
	// Convert from big endian to native


        long sentno;
        if (bswap) 	
            sentno = swapb_i32_i32(buffer[0].l);
        else
            sentno = buffer[0].l;

        if (sentno != last_sentno)
        {
            last_sentno++;

	    if (last_sentno == (long) total_sents)
	    {
		break;
	    }
            else if (sentno != last_sentno)
            {
		error("Non-sequential sentence numbers in "
			 "pfile '%s', "
			 "sentno=%li next_sentno=%li, abs_frameno=%li - "
			 "probably a corrupted PFile.",
			 filename,
			 sentno, last_sentno, abs_frameno);
            }
	    else
	    {
		sentind[sentno] = (UInt32) abs_frameno;
		next_frameno = -1;
	    }
        }

	long frameno;

        if (bswap)
          frameno = swapb_i32_i32(buffer[1].l);
        else
          frameno = buffer[1].l; 


        next_frameno++;
        if (frameno != next_frameno)
        {
	    error("Incorrect frame number in pfile '%s', "
		     "sentno=%li frameno=%li next_frameno=%li abs_frameno=%li"
		     " - probably a corrupted PFile.", filename,
		     sentno, frameno, next_frameno, abs_frameno);
	}
	abs_frameno++;
    } while(1);

    if (last_sentno!=(long) total_sents)
    {
	error("Not enough sentences in pfile '%s', "
		 "header says %lu, file has %li - probably a corrupted "
		 "PFile.", filename, total_sents, last_sentno);
    }
    // Need to add one extra index so we can calculate the length of the last
    // sentence.
    sentind[total_sents] = (UInt32) abs_frameno;

    // Check that the index of the frame after the last one is the same
    // as the number of frames in the PFile
    if ((unsigned long) abs_frameno != (unsigned long) total_frames)
    {
	error("Last sentence index (%lu) does not correspond"
		 " with number of frames (%lu) in PFile '%s' - probably a"
		 " corrupted PFile.", (unsigned long) abs_frameno,
		 (unsigned long) total_frames, filename);
    }
}

// Read one frame from the PFile, setting "buffer", "pfile_sent" and
// "pfile_frame"

inline void
InFtrLabStream_PFile::read_frame()
{
    int ec;			// Return code

    ec = fread((char *) buffer, bytes_in_row, 1, file);
    if (ec!=1)
    {
	if (feof(file))
	{ 
	    pfile_sent = SENT_EOF;
	    pfile_frame = 0;
	}
	else
	{
	    error("Failed to read frame from PFile"
		     " '%s', sent=%li frame=%li row=%li file_offset=%li.",
		     filename, current_sent, current_frame,
		     current_row, ftell(file));
	}
    }
    else
    {

        if (bswap) {
     	   pfile_sent = swapb_i32_i32((intv_int32_t) buffer[0].l);
	   pfile_frame = swapb_i32_i32((intv_int32_t) buffer[1].l);
        }
        else {
           pfile_sent =   buffer[0].l;
           pfile_frame =  buffer[1].l;
        }

    }
}

// Move to the start of the PFile
int
InFtrLabStream_PFile::rewind()
{
    // Move to the start of the data and initialise our own file offset.
    if (fseek(file, data_offset, SEEK_SET)!=0)
    {
	error("Rewind failed to move to start of data in "
		 "'%s', data_offset=%li - probably corrupted PFile.",
		 filename, data_offset);
    }
    current_sent = -1;
    current_frame = 0;
    current_row = 0;
    // Read the first frame from the PFile
    read_frame();
    return 0;			// Should return senence ID
}

size_t
InFtrLabStream_PFile::read_ftrslabs(size_t frames, float* ftrs,
				       UInt32* labs)
{
    size_t count;		// Count of number of frames

    for (count = 0; count < frames; count++)
    {
	if (pfile_sent == current_sent)
	{
	    if (pfile_frame==current_frame)
	    {
		if (ftrs!=NULL)
		{
                    if (bswap) 
        		    swapb_vi32_vi32(num_ftr_cols, 
				   (const intv_int32_t *)(&(buffer[first_ftr_col].f)),(intv_int32_t *)ftrs);
                            else
				copy_vi32_vi32(num_ftr_cols,
                                   (const intv_int32_t *)(&(buffer[first_ftr_col].f)),(intv_int32_t *)ftrs);
		    ftrs += num_ftr_cols;
		}
		if (labs!=NULL)
		{
                    if (bswap)
		      swapb_vi32_vi32(num_lab_cols,
			     (const intv_int32_t *) &(buffer[first_lab_col].l),
			     (intv_int32_t *) labs);
                    else
                      copy_vi32_vi32(num_lab_cols,
                             (const intv_int32_t *) &(buffer[first_lab_col].l),
                             (intv_int32_t *) labs);
		    labs += num_lab_cols;
		}
		read_frame();
		current_frame++;
		current_row++;
	    }
	    else
	    {
		error("Inconsistent frame number in PFile '%s',"
			 " sentence=%li frame=%li - probably corrupted PFile.",
			 filename, current_sent, current_frame);
	    }
	}
	else
	{
	    // Different sentence number - simply return number of frames
	    // so far
	    break;
	}
    }
    return count;
}

SegID
InFtrLabStream_PFile::nextseg()
{
    int ret;			// Return value

    // Skip over existing frames in sentence
    size_t skip_count = 0;	// Number of frames we have skipped
    while (current_sent==pfile_sent)
    {
	    if (pfile_frame==current_frame)
	    {
		current_frame++;
		current_row++;
		skip_count++;
		read_frame();
	    }
	    else
	    {
		error("Inconsistent frame number in PFile '%s',"
			 " sentence=%li frame=%li - probably corrupted PFile.",
			 filename, current_sent, current_frame);
	    }
    }
    if (skip_count!=0)
    {
    }

    // If possible, check that the end of sentence ties up with the index
    if (indexed)
    {
	if (current_row != (long) sentind[current_sent+1])
	    error("Sentence data inconsistent with index in"
		     " PFile '%s' - sentence %li should end at row %li  "
		     "but instead ends at row %li", filename,
		     current_sent,
		     (long) sentind[current_sent+1],
		     current_row);
	ret = SEGID_BAD;
    }
    // If we are at the end of file, exit gracefully
    if (current_sent==(long) (total_sents-1))
    {
	ret = SEGID_BAD;
    }
    else
    {
	current_sent++;
	if (current_sent!=pfile_sent)
	{
	    error("Sentence %li in PFile '%s' has sentence "
		     "number %li - probably corrupted PFile.",
		     current_sent, filename, pfile_sent);
	}
	current_frame = 0;
	if (current_frame!=pfile_frame)
	{
	    error("Sentence %li in PFile '%s' has first frame "
		     "number %li, should be 0 - probably corrupted PFile.",
		     current_sent, filename, pfile_frame);
	}
	ret = 0;		// FIXME - should return proper segment ID
    }
    return ret;
}

SegID
InFtrLabStream_PFile::set_pos(size_t segno, size_t frameno)
{
    int ret;			// Return value

    assert(segno < total_sents); // Check we do not seek past end of file

    if (indexed)
    {
	long this_sent_row;	// The number of the row at the start of sent.
	long next_sent_row;	// The row at the start of next sent.
	long row;		// The number of the row we require
	long offset;		// The position as a file offset

	this_sent_row = sentind[segno];
	row = sentind[segno] + frameno;
	next_sent_row = sentind[segno+1];
	if (row > next_sent_row)
	{
	    error("Seek beyond end of sentence %li.",
		     (unsigned long) segno);
	}
	offset = bytes_in_row * row + data_offset;
	
	if (fseek(file, offset, SEEK_SET)!=0)
	{
	    error("Seek failed in PFile "
		     "'%s', offset=%li - file problem?",
		     filename, offset);
	}
	current_sent = segno;
	current_frame = frameno;
	current_row = row;

	// Read the frame from the PFile
	read_frame();
	ret = 0;		// FIXME - should return sentence ID
    }
    else
    {
	// Tried to seek when not indexed
	ret = SEGID_BAD;
    }
    return ret;
}


int
InFtrLabStream_PFile::get_pos(size_t* segnop, size_t* framenop)
{
    size_t segno;		// The segno value returned.
    size_t frameno;		// The frameno value returned.

    if (current_sent==-1)
    {
	segno = SIZET_BAD;
	frameno = SIZET_BAD;
    }
    else
    {
	segno = (size_t) current_sent;
	frameno = (size_t) current_frame;
	assert(current_frame>=0);
    }
    if (segnop!=NULL)
	*segnop = segno;
    if (framenop!=NULL)
	*framenop = frameno;
    return OK;
}

size_t
InFtrLabStream_PFile::num_frames() {
  
  return total_frames;
}

size_t
InFtrLabStream_PFile::num_frames(unsigned int segno)
{
    size_t num_frames;		// Number of frames returned

    assert(segno<total_sents || segno==ALL);

    if (segno==ALL)
    {
	num_frames = total_frames;
    }
    else
    {
	if (indexed)
	{
	    num_frames = sentind[segno+1] - sentind[segno];
	}
	else
	    num_frames = SIZET_BAD;
    }
    
    return num_frames;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// OutFtrLabStream_PFile
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

OutFtrLabStream_PFile::OutFtrLabStream_PFile(int a_debug,
						   const char* a_filename,
						   FILE* a_file,
						   size_t a_ftrs,
						   size_t a_labs,
						   int a_indexed, 
                                                   short swap)
  : 
    file(a_file),
    filename(a_filename),
    indexed(a_indexed),
    num_ftr_cols(a_ftrs),
    num_lab_cols(a_labs),
    current_sent(0),
    current_frame(0),
    current_row(0),
    index(NULL),
    index_len(0),
    bswap(swap)
{
    int ec;			// Error code.

    // Move to the start of the PFile data section
    ec = fseek(file, PFILE_HEADER_SIZE, SEEK_SET);
    if (ec!=0)
    {
	error("Failed to seek to data section in output PFile "
		   "'%s' - cannot write PFiles to streams.",
		   filename );
    }

    // Allocate a buffer for one frame.
    buffer = new PFile_Val[num_ftr_cols + num_lab_cols + 2];
    // On SPERT, it is safer to malloc big things early.
    // We could do this later, but do not want to crash with a memory
    // error after writing a big PFile.
    header = new char[PFILE_HEADER_SIZE];
    if (indexed!=0)
    {
	index = new UInt32[DEFAULT_INDEX_SIZE];
	index_len = DEFAULT_INDEX_SIZE;
	index[0] = 0;
    }
}

OutFtrLabStream_PFile::~OutFtrLabStream_PFile()
{
    // doneseg() should be called before closing pfile - not doing this
    // is an error.
    if (current_frame!=0)
    {
	error("PFile '%s' closed mid sentence.",
		  filename);
    }
    if (indexed)
    {
	write_index();
    }
    write_header();
    if (index!=NULL)
	delete[] index;
    delete[] header;
    delete[] buffer;
}

void
OutFtrLabStream_PFile::write_ftrslabs(size_t frames, const float* ftrs,
				      const UInt32* labs)
{
    size_t i;			// Local counter.
    const size_t num_ftrs = num_ftr_cols; // Local version of value in object.
    const size_t num_labs = num_lab_cols; // Local version of value in object.
    const size_t cols = num_ftrs + num_labs + 2;
    const size_t bytes_in_frame = sizeof(PFile_Val) * cols;

    // Note we convert all data to big endian when we write it.
    for (i=0; i<frames; i++)
      {
	int ec;			// Return code.
	
        if (bswap) {
	  buffer[0].l = swapb_i32_i32(current_sent);
	  buffer[1].l = swapb_i32_i32(current_frame);
        }
        else {
          buffer[0].l = current_sent;
          buffer[1].l = current_frame;
        }
	if (num_ftrs!=0)
	{
	    if (ftrs!=NULL)
	    {
               if (bswap) {
		swapb_vi32_vi32(num_ftrs, (const intv_int32_t *)ftrs, (intv_int32_t *)&(buffer[2].f));
                }
                else {
                   copy_vi32_vi32(num_ftrs, (const intv_int32_t *)ftrs, (intv_int32_t*)&(buffer[2].f));
                }
		ftrs += num_ftrs;
	    }
	    else
		copy_f_vf(num_ftrs, 0.0f, &(buffer[2].f));
	}
	if (num_labs!=0)
	{
	    if (labs!=NULL)
	    {
                if (bswap) 
		  swapb_vi32_vi32(num_labs, (const intv_int32_t*) labs,
			       &(buffer[2+num_ftrs].l));
                else
                  copy_vi32_vi32(num_labs, (const intv_int32_t*) labs,
                               &(buffer[2+num_ftrs].l));
	    	  labs += num_labs;
	    }
	    else
		copy_i32_vi32(num_labs, 0, &(buffer[2+num_ftrs].l));
	}
	ec = fwrite((char*) buffer, bytes_in_frame, 1, file);

	if (ec!=1)
	{
	    error("Failed to write frame to PFile '%s' - only written %i items",
		       filename,ec); 
	}
	current_frame++;
	current_row++;
    }
}

void
OutFtrLabStream_PFile::doneseg(SegID)
{
    if (current_frame==0)
	error("wrote zero length sentence.");
    current_sent++;
    current_frame = 0;
    // Update the index if necessary.
    if (indexed)
    {
	// If the index is not large enough, make it bigger.
	if ((size_t) current_sent>=index_len)
	{
	    size_t new_index_len = index_len * 2;
	    UInt32* new_index = new UInt32[new_index_len];
	    copy_vi32_vi32(index_len, (intv_int32_t*) index, 
			   (intv_int32_t*) new_index);
	    delete[] index;
	    index_len = new_index_len;
	    index = new_index;
	}
	// Update the index.
	index[current_sent] = current_row;
    }
}

void
OutFtrLabStream_PFile::write_index()
{
    // write it out in original byte order
    if (bswap) 
      swapb_vi32_vi32(current_sent+1, (intv_int32_t*) index,
		   (intv_int32_t*) index);
    fwrite(index, (current_sent+1) * sizeof(UInt32), 1, file);

    // Swap it back just in case we want to use it again

    swapb_vi32_vi32(current_sent+1, (intv_int32_t*) index,
		   (intv_int32_t*) index);
}

void
OutFtrLabStream_PFile::write_header()
{
    int chars = 0;		// Number of characters added to header this
				// printf
    int count = 0;		// Total number of characters in header so far
    char* ptr = NULL;		// Point into header array
    size_t i;			// Local counter
    int ec;			// Error code

    // Unused sections of the header should be filled with \0.
    memset(header, '\0', PFILE_HEADER_SIZE);
    ptr = header;

    // Note - some sprintfs are broken - cannot use return value.

    // The version string.
    sprintf(ptr, "%s\n", pfile_version0_string);
    chars = strlen(ptr);
    count += chars; ptr += chars;

    // "Vertical" information.
    // -num_sentences
    sprintf(ptr, "-num_sentences %lu\n", (unsigned long) current_sent);
    chars = strlen(ptr);
    count += chars; ptr += chars;
    // -num_frames
    sprintf(ptr, "-num_frames %lu\n", (unsigned long) current_row);
    chars = strlen(ptr);
    count += chars; ptr += chars;

    // Feature information.
    sprintf(ptr, "-first_feature_column %lu\n", (unsigned long) 2);
    chars = strlen(ptr);
    count += chars; ptr += chars;
    sprintf(ptr, "-num_features %lu\n", (unsigned long) num_ftr_cols);
    chars = strlen(ptr);
    count += chars; ptr += chars;

    // Label information.
    sprintf(ptr, "-first_label_column %lu\n",
	    (unsigned long)  2+num_ftr_cols);
    chars = strlen(ptr);
    count += chars; ptr += chars;
    sprintf(ptr, "-num_labels %lu\n", (unsigned long) num_lab_cols);
    chars = strlen(ptr);
    count += chars; ptr += chars;

    // The format string.
    sprintf(ptr, "-format dd");
    chars = strlen(ptr);
    count += chars; ptr += chars;
    for (i=0; i<num_ftr_cols; i++)
	*ptr++ = 'f';
    for (i=0; i<num_lab_cols; i++)
	*ptr++ = 'd';
    count += (num_ftr_cols + num_lab_cols);
    *ptr++ = '\n';
    count++;

    // The details of the data sections.
    size_t cols = num_ftr_cols + num_lab_cols + 2;
    size_t data_size = cols * current_row;
    sprintf(ptr, "-data size %lu offset %lu ndim %lu nrow %lu ncol %lu\n",
	    (unsigned long) data_size, (unsigned long) 0, (unsigned long) 2,
	    (unsigned long) current_row, (unsigned long) cols);
    chars = strlen(ptr);
    count += chars; ptr += chars;

    // If necessary, details of the sentence index.
    if (indexed)
    {
	size_t sentind_size = current_sent + 1;
	sprintf(ptr, "-sent_table_data size %lu offset %lu ndim 1\n",
		(unsigned long) sentind_size, (unsigned long) data_size);
	chars = strlen(ptr);
	count += chars; ptr += chars;
    }

    // The end of the header.
    sprintf(ptr, "-end\n");
    chars = strlen(ptr);
    count += chars; ptr += chars;

    assert((unsigned long) count<=PFILE_HEADER_SIZE);

    // Seek to start of file to write header.
    ec = fseek(file, 0L, SEEK_SET);
    if (ec!=0)
    {
	error("Failed to seek to start of PFile '%s' - %s.",
		   filename, strerror(ec));
    }
    ec = fwrite(header, PFILE_HEADER_SIZE, 1, file);
    if (ec!=1)
    {
	error("Failed to write header in PFile '%s' - %s.",
		   filename, strerror(ec));
    }
}

