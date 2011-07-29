
%option nounput

%{

/*
 * This is a simple Lex filter to put the Viterbi printing output
 * into a canonical form so that the different printing 
 * implementations can be diffed against each other or known-correct
 * output.
 *
 * -vitValsFile and -pVitVals file should generate the same 
 * value assignments for each segment, but they may be organized
 * differently. This filter just removes the organization...
 * Both printing implementations generate a line for each partition
 * consisting of comma separated variableName(frameNumber)=value
 * assignments. The filter just puts each of those on a separate
 * line, prefixed by its segment number. If that's piped through
 * sort, the results should be identical for -vitVals and -pVitVals
 */

#include <stdlib.h>

int segmentNumber;

%}

ws	[[:space:]]+
ident   [[:alpha:]]([[:alnum:]]|\_|\-)*
int	[[:digit:]]+
vitval	{ident}\({int}\)=[^,\n]+

%s segnum
%%

{ws}

Segment		{ BEGIN(segnum); }
<segnum>{int}	{ segmentNumber = atoi(yytext); BEGIN(INITIAL); }

{vitval}	{ printf("seg %d: %s\n", segmentNumber, yytext); }
.

%%

int
main(int argc, char *argv[]) {
  int rc = 0;
  if (argc > 1) {
    for (int i=1; i < argc; i+=1) {
      yyin = fopen( argv[i], "r" );
      rc = yylex();
      if (rc) goto end;
    }
  } else {
    yyin = stdin;
    rc = yylex();
  }
end:
  exit(rc);
}
