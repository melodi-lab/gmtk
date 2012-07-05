/*  Generated header
 *  File Name : GMTK_WordOrganization.cc
 *
 *  Created   : 2003-11-19 10:53:43 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
 *
 * $Header$
 *
*/

#include <cstdio>
#include "GMTK_WordOrganization.h"


ByteEndian getWordOrganization() {

  ByteEndian  byteEndian;
  ShortEndian shortEndian;

  union { unsigned long ln; unsigned char chrs[4]; unsigned short shrs[2]; };
  ln = 0xAABBCCDD;

  if (shrs[1] == 0xAABB) {
    shortEndian = SHORT_LITTLE_ENDIAN;
    if (chrs[3] == 0xAA) {
      byteEndian = BYTE_LITTLE_ENDIAN;
    } else if (chrs[2] == 0xAA) {
      byteEndian = BYTE_BIG_ENDIAN;
    } else {
      byteEndian = BYTE_UNSUPPORTED_ORG;
    }
  } else if (shrs[0] == 0xAABB) {
    shortEndian = SHORT_BIG_ENDIAN;
    if (chrs[1] == 0xAA) {
      byteEndian = BYTE_LITTLE_ENDIAN;
   } else if (chrs[0] == 0xAA) {
     byteEndian = BYTE_BIG_ENDIAN;
   } else {
     byteEndian = BYTE_UNSUPPORTED_ORG;
   }
  } else {
    shortEndian = SHORT_UNSUPPORTED_ORG;
    byteEndian = BYTE_UNSUPPORTED_ORG;
  }

  return byteEndian;
}

#ifdef TEST_WO
int main(void) {

  ByteEndian byteEndian = getWordOrganization();

  printf("The word organization on this machine is %s.\n",(byteEndian==BYTE_BIG_ENDIAN)?"BIG ENDIAN":(byteEndian==BYTE_LITTLE_ENDIAN)?"LITTLE_ENDIAN":"UNSUPPORTED");

  return 0;

}
#endif
