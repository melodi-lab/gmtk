%option nounput

%{

/*
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 */

%}

%%

" \n"	printf("\n");
.       ECHO;

%%

int
main(int argc, char *argv[]) {
  yylex();
  exit(0);
}
