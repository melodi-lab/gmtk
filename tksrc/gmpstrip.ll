%option nounput noyywrap

%{

/*
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 * This program looks for lines that end with " \n"
 * and normalizes them to "\n". This is to implement
 * a test case for https://j.ee.washington.edu/trac/gmtk/ticket/161
 * which needs to compare the EM parameters learned by
 * the trunk and branch versions of gmtkEMtrain. Unfortunately,
 * the two versions differ occasionally in their line endings,
 * so this filters out the EOL "noise" from the diff.
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
