/* 
 * A lex file for scanning a GM structure file
 * 
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (c) 2001, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this software
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 *
 *  $Header$
 *
 */

%{
/* need this for the call to atof() below */
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_FileParser.h"

#ifndef _AIX
extern "C" int fileno(FILE*); 
#endif

int debugLexer = 0;

%}

/* General string separated by "" characters, on a single line. */
string  \"[^\n"]+\"

/* White space */
ws      [ \t]+

/* Basic alpha and numeric stuff */
alpha   [A-Za-z]
dig     [0-9]

/* A valid identifier */ 
ident   ({alpha})({alpha}|{dig}|\_|\-)*

/* numeric stuff */
unsigned {dig}+
int     [-+]?{unsigned}
flt1    [-+]?{dig}+\.?([eE][-+]?{dig}+)?
flt2    [-+]?{dig}*\.{dig}+([eE][-+]?{dig}+)?
flt     {flt1}|{flt2}

/* integer ranges */
int_rng {int}:{int}


/* all valid keywords */
keyword GRAPHICAL_MODEL|frame|variable|type|cardinality|switchingparents|conditionalparents|discrete|continuous|hidden|observed|value|nil|using|mapping|collection|DenseCPT|SparseCPT|DeterministicCPT|mixGaussian|gausSwitchMixGaussian|logitSwitchMixGaussian|mlpSwitchMixGaussian|chunk

separator ":"|";"|"{"|"}"|"("|")"|"|"|","

%%

"%"[^\n]*   /* eat up one-line comments */

  /*
     ***************************************
     ** support for cpp line directives    *
     ***************************************
  */

   /* parse CPP file/line directives to get the #include correct */
   /* ^"#"{ws}{int}{ws}{string}{ws}{int} { } */


^"#".* {
            /* eat up cpp any and all possible file/line directives for now. 
               Ultimately parse this */
                ;
}

[ \t]+    /* eat up just whitespace */

\n        { FileParser::tokenInfo.srcLine++; }

{string}     {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::TT_String;
            if (debugLexer)
              printf( "A string: %s\n", yytext);
	    return FileParser::tokenInfo.tokenType;
            }


{int}      {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::TT_Integer;
            FileParser::tokenInfo.int_val = atoi( yytext );
            if (debugLexer)
              printf( "An integer: %s (%d)\n", yytext,
                    atoi( yytext ) );
	    return FileParser::tokenInfo.tokenType;
            }



{flt}       {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::TT_Real;
            FileParser::tokenInfo.float_val = atof( yytext );
            if (debugLexer)
              printf( "A float: %s (%g)\n", yytext,
                    atof( yytext ) );
	    return FileParser::tokenInfo.tokenType;
            }

{separator} {
            FileParser::tokenInfo.tokenStr = yytext;
            if (debugLexer)
              printf( "A separator: %s\n", yytext );	     	    	       
            switch(*yytext) {
                 case ':': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_Colon;
		   return FileParser::tokenInfo.tokenType;
		   break;
                 case ',': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_Comma;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case ';': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_SemiColon;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '{': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_LeftBrace;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '}': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_RightBrace;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '(': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_LeftParen;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case ')': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_RightParen;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '|': 
		   FileParser::tokenInfo.tokenType = FileParser::TT_VirtBar;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 default: { abort(); }
             }
            }

{keyword}   {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::TT_Keyword;
            if (debugLexer)
              printf( "A keyword: %s\n", yytext );
	    return FileParser::tokenInfo.tokenType;
            }

{ident}     { 
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::TT_Identifier;
            if (debugLexer)
	      printf( "An identifier: %s\n", yytext );
	    return FileParser::tokenInfo.tokenType;
            }

.           { 
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::TT_Undefined;
            if (debugLexer)
              printf( "Unrecognized character: %s\n", yytext );
	    return FileParser::tokenInfo.tokenType;
            }

%%

#if 0
int
main( int argc, char **argv )
{
    ++argv, --argc;  /* skip over program name */
    if ( argc > 0 )
            yyin = fopen( argv[0], "r" );
    else
            yyin = stdin;

    while (1) {
       int rc = yylex();
       printf("rc = %d\n",rc);
    }

}

#endif
