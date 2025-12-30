%{
#include "y.tab.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
%}

%x COMMENT

DIGIT   [0-9]
NUMBER  {DIGIT}+(\.{DIGIT}+)?

%%

"/*"            { BEGIN(COMMENT); }
<COMMENT>"*/"   { BEGIN(INITIAL); }
<COMMENT>\n     { }
<COMMENT>.      { }
"//".* { }

"matrix"        return MATRIX;
"transpose"     return TRANSPOSE;
"det"           return DET;
"trace"         return TRACE;
"eigenval"      return EIGENVAL;
"eigenvec"      return EIGENVEC;
"inverse"       return INVERSE;
"solve"         return SOLVE;
"cholesky"      return CHOLESKY;
"eye"           return EYE;
"zeros"         return ZEROS;
"symmetric"     return SYMMETRIC;
"sparse"        return SPARSE;
"pos_def"       return POS_DEF;
"dense"         return DENSE;
"norm"          return NORM;
"cond"          return COND;
"rank"          return RANK;

"="             return '=';
";"             return ';';
","             return ',';
"["             return '[';
"]"             return ']';
"+"             return '+';
"-"             return '-';
"*"             return '*';
"("             return '(';
")"             return ')';
".*"            return DOT_MUL;

{NUMBER} {
    yylval.dval = atof(yytext);
    return NUMBER;
}

[a-zA-Z_][a-zA-Z0-9_]* {
    yylval.sval = strdup(yytext);
    return ID;
}

[ \t\r\n]+      ;

.               { printf("Lexer Error: Unknown character '%s'\n", yytext); }

%%

int yywrap(void) {
    return 1;
}