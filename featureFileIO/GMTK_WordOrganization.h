/*  Generated header
 *  File Name : GMTK_WordOrganization.h
 *
 *  Created   : 2003-11-19 10:53:43 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
 *
 * $Header$
 *
*/

#ifndef GMTK_WORDORGANIZATION_H
#define GMTK_WORDORGANIZATION_H

// Big Endian: byte with address 0xXXXXXX00 is in the most significant (big end)
//               of a 32 bit word.
// Little Endian: byte with address 0xXXXXXX00 is in the least significant 
//               (little end) of a 32 bit word.


// ************ ASSUMPTIONS (that we shouldn't make) ***************
//  byte = 8 bits
//  short = 2 bytes
//  long = 4 bytes


enum ByteEndian { BYTE_BIG_ENDIAN,BYTE_LITTLE_ENDIAN,BYTE_UNSUPPORTED_ORG };
enum ShortEndian { SHORT_BIG_ENDIAN,SHORT_LITTLE_ENDIAN,SHORT_UNSUPPORTED_ORG };


ByteEndian getWordOrganization();


#endif
