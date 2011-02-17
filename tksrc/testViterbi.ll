
%option nounput

%{

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
