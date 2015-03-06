/*
 * Copyright (C) 2014 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

/*

This program parses the ASCII gmtkViterbi output and 
picks out the actual Viterbi value lines. It copies 
them to stdout with the variables in sorted order 
(frame # as major key, variable name as minor). This
is to support the GMTK test suite: within a gmtkViterbi
run, the variables should come out in the same order.
But between different gmtkViterbi invocations, the 
order can change. With the variables in sorted order,
the output from different gmtkViterbi invocations can
be compared.

*/

%option nounput
%option noinput
%option noyywrap

%{

#define MAX_PREFIX_LENGTH 20
char prefix[MAX_PREFIX_LENGTH];

#define MAX_VARNAME_LENGTH 1024
char varname[MAX_VARNAME_LENGTH];

int frame;

#define MAX_VALUE_LENGTH 1024
char value[MAX_VALUE_LENGTH];

#define VITVALS 260

#include <string>
#include <vector>
#include <algorithm>
using namespace std;


class VitValue {
 public:
  string prefix;
  string variableName;
  int    frameNum;
  string value;

  VitValue(char const *prefix, char const *varName, int frameNum, char const* val) 
   : prefix(prefix), variableName(varName), frameNum(frameNum), value(val)
  {}

  bool operator< (VitValue const &vv) const {
    if (frameNum < vv.frameNum) return true;
    if (frameNum > vv.frameNum) return false;
    if (variableName < vv.variableName) return true;
    return false;
  }
};

vector<VitValue> vitVals;

%}

%x VARNAME
%x FRAMENUM
%x VITVAL

%%

Ptn-[[:digit:]]+" "[PCE]'?:      strcpy(prefix, yytext); BEGIN(VARNAME);
<VARNAME>" "
<VARNAME>[^ (]+                  strcpy(varname, yytext); 
<VARNAME>"("                     BEGIN(FRAMENUM);

<FRAMENUM>[[:digit:]]+           frame = atoi(yytext);
<FRAMENUM>")="                   BEGIN(VITVAL);

<VITVAL>[^,\n]+                  strcpy(value, yytext);
<VITVAL>","                      BEGIN(VARNAME); {VitValue vv(prefix, varname, frame, value); vitVals.push_back(vv); }
<VITVAL>\n                       BEGIN(INITIAL); {VitValue vv(prefix, varname, frame, value); vitVals.push_back(vv); } return VITVALS;

.
\n

%%
   
int
main(int argc, char *argv[]) {
  int tok;
  int i;
  for (i=1; i < argc; i+=1) {
    if (strcmp(argv[i], "-") == 0) {
      yyin = stdin;
    } else {
      yyin = fopen(argv[i], "r");
      if (!yyin) {
        perror(argv[i]);
        continue;
      }
    }
    for (tok=yylex(); tok; tok=yylex()) {
      switch (tok) {
      case VITVALS:
        sort(vitVals.begin(), vitVals.end());
        vector<VitValue>::iterator it = vitVals.begin();
	printf("%s %s(%d)=%s", it->prefix.c_str(), it->variableName.c_str(), 
                               it->frameNum, it->value.c_str());
        for(++it; it != vitVals.end(); ++it) 
          printf(",%s(%d)=%s", it->variableName.c_str(), it->frameNum, it->value.c_str());
        printf("\n");
	vitVals.clear();
	break;
      }
    }
    if (fclose(yyin) != 0) {
      perror(argv[i]);
    }
  }
  exit(0);
}

