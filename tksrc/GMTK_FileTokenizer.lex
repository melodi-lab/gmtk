/* 
 *  A lex file for scanning a GM structure file
 *  $Header$
 *
 */

%{
/* need this for the call to atof() below */
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

#include "GMTK_FileParser.h"

extern "C" int fileno(FILE*); 

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
int     {dig}+
flt1    [-+]?{dig}+\.?([eE][-+]?{dig}+)?
flt2    [-+]?{dig}*\.{dig}+([eE][-+]?{dig}+)?
flt     {flt1}|{flt2}

/* integer ranges */
int_rng {int}:{int}


/* all valid keywords */
keyword GRAPHICAL_MODEL|frame|variable|type|disposition|cardinality|switchingparents|conditionalparents|discrete|continuous|hidden|observed|nil|using|mapping|mixGaussian|gausSwitchMixGaussian|logitSwitchMixGaussian|mlpSwitchMixGaussian|chunk

separator ":"|";"|"{"|"}"|"("|")"|"|"

%%


"#"[^\n]*   /* eat up one-line comments */


[ \t]+    /* eat up whitespace */

\n        { FileParser::tokenInfo.srcLine++; }


{int}+      {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::Token_Integer;
            FileParser::tokenInfo.int_val = atoi( yytext );
            if (debugLexer)
              printf( "An integer: %s (%d)\n", yytext,
                    atoi( yytext ) );
	    return FileParser::tokenInfo.tokenType;
            }


{flt}       {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::Token_Real;
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
		   FileParser::tokenInfo.tokenType = FileParser::Token_Colon;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case ';': 
		   FileParser::tokenInfo.tokenType = FileParser::Token_SemiColon;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '{': 
		   FileParser::tokenInfo.tokenType = FileParser::Token_LeftBrace;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '}': 
		   FileParser::tokenInfo.tokenType = FileParser::Token_RightBrace;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '(': 
		   FileParser::tokenInfo.tokenType = FileParser::Token_LeftParen;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case ')': 
		   FileParser::tokenInfo.tokenType = FileParser::Token_RightParen;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 case '|': 
		   FileParser::tokenInfo.tokenType = FileParser::Token_VirtBar;
		   return FileParser::tokenInfo.tokenType;
		   break;
		 default: { abort(); }
             }
            }

{keyword}   {
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::Token_Keyword;
            if (debugLexer)
              printf( "A keyword: %s\n", yytext );
	    return FileParser::tokenInfo.tokenType;
            }

{ident}     { 
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::Token_Identifier;
            if (debugLexer)
	      printf( "An identifier: %s\n", yytext );
	    return FileParser::tokenInfo.tokenType;
            }

.           { 
            FileParser::tokenInfo.tokenStr = yytext;
            FileParser::tokenInfo.tokenType = FileParser::Token_Undefined;
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
