/* 
 *  A lex file for scanning a GM structure file
 *  $Header$
 *
 */

%{
/* need this for the call to atof() below */
#include <math.h>

#include <stdio.h>


#include "GMTK_FileParser.h"


%}

/* General string */
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

/* all valid keywords */
keyword model|frame|variable|disposition|parents|attributes|parameters

separator ":"|";"

%%


"#"[^\n]*   /* eat up one-line comments */

[ \t\n]+    /* eat up whitespace */


{int}+    {
            printf( "An integer: %s (%d)\n", yytext,
                    atoi( yytext ) );
	    return FileParser::TT_INT;
            }


{flt}        {
            printf( "A float: %s (%g)\n", yytext,
                    atof( yytext ) );
	    return FileParser::TT_FLT;
            }

{keyword}   {
            printf( "A keyword: %s\n", yytext );
	    return FileParser::TT_KWD;
            }

{ident}     { printf( "An identifier: %s\n", yytext );
              return FileParser::TT_ID; 
            }
 
{separator} {  printf( "A separator: %s\n", yytext );
              return FileParser::TT_OP;                      
            }

.           { printf( "Unrecognized character: %s\n", yytext );
              return FileParser::TT_UNDEF;
            }

%%

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

