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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"

#include "GMTK_FileParser.h"

#define YY_NEVER_INTERACTIVE 1

#ifndef _AIX
extern "C" int fileno(FILE*) throw(); 
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


/* all valid keywords 
 * This set must be consistent with the enum TokenKeyword in FileParser.h 
 * and with kw_table[] in FileParser.cc   
*/


keyword GRAPHICAL_MODEL|frame|variable|type|cardinality|switchingparents|conditionalparents|discrete|continuous|hidden|observed|observation|weight|value|nil|using|mapping|collection|DenseCPT|SparseCPT|DeterministicCPT|mixGaussian|gausSwitchMixGaussian|logitSwitchMixGaussian|mlpSwitchMixGaussian|chunk|elimination_hint

separator ":"|";"|"{"|"}"|"("|")"|"|"|","

%%

"%"[^\n]*   /* eat up one-line comments */


  /*
     ***************************************
     ** support for cpp line directives   **
     ***************************************
  */


^"#line"{ws}{int}{ws}{string}.* {
            if (debugLexer)
              printf( "A CPP line directive with 'line': %s\n", yytext);
            /* get the line number and file name */
            char *ptr = yytext;
            ptr += 6; /* skip over '#line' to first white space */
            while (*ptr && isspace(*ptr)) ptr++;
            /* ptr now points to start of integer. */
            char *p;
            long l = strtol(ptr,&p,0);
            if (p == ptr)
              error("ERROR: tokenizer couldn't obtain line number.");
            ptr = p;
            while (*ptr && isspace(*ptr)) ptr++;
            if (*ptr != '\"')
              error("ERROR: tokenizer expecting file name.");
            p = ++ptr;
            while (*ptr && *ptr != '\"') ptr++;
            if (*ptr != '\"')
              error("ERROR: tokenizer expecting end of file name.");
            *ptr = '\0';
            FileParser::fileNameParsing = p;
            FileParser::tokenInfo.srcLine = l-1;
}

^"#"{ws}{int}{ws}{string}.* {
            if (debugLexer)
              printf( "A CPP line directive: %s\n", yytext);               
            /* get the line number and file name */
            char *ptr = yytext;
            ptr++; /* skip over '#' to first white space */
            while (*ptr && isspace(*ptr)) ptr++;
            /* ptr now points to start of integer. */
            char *p;
            long l = strtol(ptr,&p,0);
            if (p == ptr)
              error("ERROR: tokenizer couldn't obtain line number.");
            ptr = p;
            while (*ptr && isspace(*ptr)) ptr++;
            if (*ptr != '\"')
              error("ERROR: tokenizer expecting file name.");
            p = ++ptr;
            while (*ptr && *ptr != '\"') ptr++;
            if (*ptr != '\"')
              error("ERROR: tokenizer expecting end of file name.");
            *ptr = '\0';
            FileParser::fileNameParsing = p;
            FileParser::tokenInfo.srcLine = l-1;
}

^"#".* {
            if (debugLexer)
              printf("Blank line\n");
            /*  eat up cpp any and all possible file/line directives for now. 
             *  Ultimately parse this 
             */
}

[ \t]+    /* eat up just whitespace */

\n        { FileParser::tokenInfo.srcLine++; 
            /* "." eats up all chars except for newline, so
               we count lines here.
             */
          }

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
            FileParser::tokenInfo.doub_val = atof( yytext );
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
