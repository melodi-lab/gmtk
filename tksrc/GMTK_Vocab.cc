/*-
 * GMTK_Vocab.cc
 *     Vocab
 *
 * Written by Gang Ji <gang@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 *  $Header: /homes/gang/research/programs/cpp/gmtk/RCS/GMTK_NGramCPT.h,v 1.1 20
04/05/13 00:59:36 gang Exp $
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

#include <cstring>

#include "GMTK_Vocab.h"
#include "fileParser.h"
#include "error.h"


/*-
 *-----------------------------------------------------------------------
 * nextPrime
 *      Find the prime number just bigger than twice the input
 *
 * Results:
 *      Return the prime number just bigger than 2*n.
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
unsigned nextPrime(unsigned n);


/*-
 *-----------------------------------------------------------------------
 * Vocab::Vocab
 *      Default constructor
 *
 * Results:
 *      Create the buffer
 *
 * Side Effects:
 *
 *-----------------------------------------------------------------------
 */
Vocab::Vocab() : _size(0), _tableSize(0), _indexTable(NULL), _stringTable(NULL), _unkIndex(~0x0) {
}


Vocab::Vocab(unsigned size) : _size(size), _unkIndex(size) {
	_tableSize = nextPrime(_size << 1 );
	if ( ! (_indexTable = new HashEntry [_tableSize]) )
		error("out of memory");

	if ( ! (_stringTable = new const char * [_size]) )
		error("out of memory");
	::memset(_stringTable, 0, sizeof(char *) * _size);
}


/*-
 * Vocab::~Vocab
 *      Default destructor
 *
 * Results:
 *      None.
 *
 * Side Effects
 *      Clean up the memory
 */
Vocab::~Vocab() {
	delete [] _indexTable;
	delete [] _stringTable;
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::index
 *      Find the index of a word
 *
 * Results:
 *      Return the index of the word, or <unk> index if not found.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
unsigned Vocab::index(const char* word) const {
	const HashEntry *pos = findPos(word);

	if ( pos && pos->key != NULL )
		return pos->wid;

	return _unkIndex;
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::word
 *      Find the string of a word from index
 *
 * Results:
 *      Return the string of the word
 *
 * Side Effects:
 *      Only the char* is returned.  Becareful not to delete [] or modify
 *      from out side.
 *
 *-----------------------------------------------------------------------
 */
const char * Vocab::word(unsigned index) const {
	if ( index >= _size )
		return "<unk>";

	return _stringTable[index];
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::read
 *      Read words from a file.
 *
 * Results:
 *      Read words and indices.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
void Vocab::read(const char *filename) {

  iDataStreamFile ifs(filename, false, false);	// ascii, no cpp

  // start allocating wids at 0 (so the first word in the
  // file has 0, the next has 1, and so on).
  unsigned wid = 0;
  // ????????????????????
  // TODO: why is ifs.prepareNext() being called? why not just call readStr()?????
  //
  while ( ifs.prepareNext() ) {
    char *word = NULL;
    if ( ! ifs.readStr(word) )
      error("Error: cannot read string when reading vocab file '%s'.",filename);
    insert(word, wid++);
    if ( wid > _size )
      error("Error: vocab file file '%s' contains more words than specified vocab size of %d",
	    filename,_size);
    // readStr() allocates its own buffer, so we free it here.
    delete [] word;		
  }
  
  if ( wid < _size )
    error("Error: in vocab with cardinality %d, read only %d words from file '%s'.", _size, wid, filename);
}



/*-
 *-----------------------------------------------------------------------
 * Vocab::read
 *      read in a new vocab specification into the object.
 *
 * Results:
 *      none
 *
 * Side Effects:
 *      The new object 'this' is now instantiated if no error occurs.
 *
 *-----------------------------------------------------------------------
 */
void Vocab::read(iDataStreamFile& is) {
  NamedObject::read(is);
  is.read(_size, "Can't read Vocab's cardinality (i.e., the lexicon size)");

  resize(_size);

  char *vocabFile;
  if ( !is.readStr(vocabFile) )
    error("ERROR: reading file '%s' line %d, Vocab '%s' can't read vocab filename", 
	  is.fileName(), is.lineNo(),name().c_str());
  
  read(vocabFile);

  delete [] vocabFile;
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::resize
 *      Resize the vocab buffer.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      All previous information will be erased.
 *
 *-----------------------------------------------------------------------
 */
void Vocab::resize(unsigned card) {
	_unkIndex = _size = card;

	delete [] _indexTable;
	_tableSize = nextPrime(_size << 1 );
	if ( ! (_indexTable = new HashEntry [_tableSize]) )
		error("trying to resize vocab to be of size %d, but out of memory",_tableSize);

	delete [] _stringTable;
	if ( ! (_stringTable = new const char * [_size]) )
		error("out of memory");
	::memset(_stringTable, 0, sizeof(char *) * _size);
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::insert
 *      Insert a word into the vocabulary.
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      If the word is "<unk>", the <unk> index will be changed.
 *
 *-----------------------------------------------------------------------
 */
void Vocab::insert(const char* key, unsigned wid) {

  if ( wid >= _size )
    error("Error: word id %d exeeds size %d in Vocab object", wid, _size);

  // Insert x as active
  HashEntry *pos = findPos(key);

  // check for duplicates
  if ( pos->key != NULL )		
    error("Error: in vocab inserting word %s more than once", key);

  pos->wid = wid;
  // delete [] pos->key; no need
  pos->key = new char [strlen(key) + 1];
  strcpy(pos->key, key);

  if ( strcmp(key, "<unk>") == 0 )
    _unkIndex = wid;
  
  _stringTable[wid] = pos->key;
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::findPos
 *      Find the hash entry for a word string.
 *
 * Results:
 *      Pointer to the entry for word string.
 *
 * Side Effects:
 *      Be sure the table is not full.
 *
 *-----------------------------------------------------------------------
 */
Vocab::HashEntry* Vocab::findPos(const char* key) const {
	unsigned pos = hash(key);
	unsigned numberOfCollisions = 0;

	while ( _indexTable[pos].key != NULL && strcmp(_indexTable[pos].key, key) ) {
		pos = (pos + (++numberOfCollisions << 1) - 1) % _tableSize;
	}

	return _indexTable + pos;
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::HashEntry::HashEntry
 *      Default constructor
 *
 * Results:
 *      None.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
Vocab::HashEntry::HashEntry(const char* theKey, int theWid) : wid(theWid) {
	key = new char [strlen(theKey) + 1];
	strcpy(key, theKey);
}


/*-
 *-----------------------------------------------------------------------
 * Vocab::HashEntry::operator =
 *      Overloaded assignment operator
 *
 * Results:
 *      The reference to itself.
 *
 * Side Effects:
 *      Deep copy data.
 *
 *-----------------------------------------------------------------------
 */
const Vocab::HashEntry& Vocab::HashEntry::operator = (const Vocab::HashEntry& he) {
	if ( &he != this ) {
		delete [] key;
		key = new char [strlen(he.key) + 1];
		strcpy(key, he.key);
		wid = he.wid;
	}

	return *this;
}


/*-
 *-----------------------------------------------------------------------
 * isPrime
 *      Determin whether a number is prime or not.
 *
 * Results:
 *      Return true if n is a prime number.
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
bool isPrime(unsigned n) {
	if ( n == 2 || n == 3 )
		return true;
	if ( n == 1 || (! (n & 0x1)) )
		return false;

	for ( unsigned i = 3; i * i <= n; i += 2 )
		if ( ! (n % i ) )
			return false;

	return true;
}


/*-
 *-----------------------------------------------------------------------
 * Function
 *      Find a prime number just bigger than 2 * n
 *
 *      TODO: this could be very slow -- do this with a table instead.
 *
 * Results:
 *      Find a prime number just bigger than 2 * n
 *
 * Side Effects:
 *      None.
 *
 *-----------------------------------------------------------------------
 */
unsigned nextPrime(unsigned n) {
	n = (n << 1) + 1;

	while ( ! isPrime(n) )
		n += 2;

	return n;
}


////////////////////////////////////////////////////////////////////
//        Test Driver
////////////////////////////////////////////////////////////////////

#ifdef MAIN


#include <cstdio>


int main() {
	Vocab vocab(13);

	vocab.read("w2.dct");

	printf("vocab size %d\n", vocab.size());

	printf("index of oh is %d\n", vocab.index("oh"));
	printf("index of foo is %d\n", vocab.index("foo"));

	printf("string of 3 is %s\n", vocab.word(3));

	return 0;
}


#endif
