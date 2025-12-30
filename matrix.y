%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CODE_BUFFER_SIZE 4096
#define TEMP_NAME_SIZE 32

typedef struct {
    char name[32];
    int rows;
    int cols;
} Matrix;

Matrix symtab[100];
int symcount = 0;

void add_matrix(char* name, int r, int c);
Matrix* find_matrix(char* name);
int yyerror(const char* s);
int yylex(void);
FILE* out;

int temp_count = 0;
char* new_temp() {
    char* s = malloc(TEMP_NAME_SIZE);
    if (!s) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    sprintf(s, "temp_%d", temp_count++);
    return s;
}

/* Helper to safely allocate code buffer */
char* alloc_code_buffer() {
    char* buf = malloc(CODE_BUFFER_SIZE);
    if (!buf) { fprintf(stderr, "Memory allocation failed\n"); exit(1); }
    return buf;
}

/* Helper to generate binary operation code */
char* gen_binary_op(char* code1, char* code3, char* temp, 
                    int rows, int cols, 
                    char* var1, char* var3, 
                    const char* op_str, const char* op_symbol) {
    char* code = alloc_code_buffer();
    snprintf(code, CODE_BUFFER_SIZE,
        "%s%s"
        "double %s[%d][%d];\n"
        "/* %s */\n"
        "for(i=0;i<%d;i++){\n"
        "  for(j=0;j<%d;j++){\n"
        "    %s[i][j] = %s[i][j] %s %s[i][j];\n"
        "  }\n"
        "}\n",
        code1, code3, temp, rows, cols, op_str,
        rows, cols, temp, var1, op_symbol, var3);
    return code;
}

double num_buffer[1024];
int num_count = 0;
%}

%union {
    double dval;
    char* sval;
    struct {
        char* code;
        char* result_var;
        int rows;
        int cols;
    } expr_val;
}

%token MATRIX TRANSPOSE DET TRACE EIGENVAL EIGENVEC INVERSE
%token SOLVE CHOLESKY EYE ZEROS SYMMETRIC SPARSE POS_DEF DENSE NORM COND RANK DOT_MUL
%token <sval> ID
%token <dval> NUMBER

%type <dval> number
%type <expr_val> expr

%left '+' '-'
%left '*'
%left TRANSPOSE

%%
program:
{
    out = fopen("output.c", "w");
    fprintf(out,
        "#include <stdio.h>\n"
        "#include <math.h>\n"
        "#include <stdlib.h>\n"
        "#include <string.h>\n\n"
        "#define MAX_N 100\n\n"
        "void lu_decompose(double A[][MAX_N], int n, double L[][MAX_N], double U[][MAX_N]) {\n"
        "    int i,j,k;\n"
        "    for(i=0;i<n;i++){\n"
        "        for(j=0;j<n;j++){\n"
        "            if(i<=j){\n"
        "                U[i][j]=A[i][j];\n"
        "                for(k=0;k<i;k++) U[i][j]-=L[i][k]*U[k][j];\n"
        "                L[i][j]=(i==j)?1:0;\n"
        "            }else{\n"
        "                L[i][j]=A[i][j];\n"
        "                for(k=0;k<j;k++) L[i][j]-=L[i][k]*U[k][j];\n"
        "                L[i][j]/=U[j][j];\n"
        "                U[i][j]=0;\n"
        "            }\n"
        "        }\n"
        "    }\n"
        "}\n\n"
        "double det_nxn(double A[][MAX_N], int n) {\n"
        "    double L[MAX_N][MAX_N], U[MAX_N][MAX_N];\n"
        "    double det = 1.0;\n"
        "    int i;\n"
        "    lu_decompose(A, n, L, U);\n"
        "    for(i=0;i<n;i++) det *= U[i][i];\n"
        "    return det;\n"
        "}\n\n"
        "int gauss_jordan(double A[][MAX_N], double inv[][MAX_N], int n) {\n"
        "    double aug[MAX_N][2*MAX_N];\n"
        "    int i,j,k,p;\n"
        "    double temp, pivot;\n"
        "    for(i=0;i<n;i++){\n"
        "        for(j=0;j<n;j++) aug[i][j]=A[i][j];\n"
        "        for(j=n;j<2*n;j++) aug[i][j]=(j-n==i)?1:0;\n"
        "    }\n"
        "    for(i=0;i<n;i++){\n"
        "        p=i;\n"
        "        for(j=i+1;j<n;j++){\n"
        "            if(fabs(aug[j][i])>fabs(aug[p][i])) p=j;\n"
        "        }\n"
        "        if(fabs(aug[p][i])<1e-10) return 0;\n"
        "        if(p!=i){\n"
        "            for(j=0;j<2*n;j++){\n"
        "                temp=aug[i][j]; aug[i][j]=aug[p][j]; aug[p][j]=temp;\n"
        "            }\n"
        "        }\n"
        "        pivot=aug[i][i];\n"
        "        for(j=0;j<2*n;j++) aug[i][j]/=pivot;\n"
        "        for(k=0;k<n;k++){\n"
        "            if(k!=i){\n"
        "                temp=aug[k][i];\n"
        "                for(j=0;j<2*n;j++) aug[k][j]-=temp*aug[i][j];\n"
        "            }\n"
        "        }\n"
        "    }\n"
        "    for(i=0;i<n;i++){\n"
        "        for(j=0;j<n;j++) inv[i][j]=aug[i][j+n];\n"
        "    }\n"
        "    return 1;\n"
        "}\n\n"
        "void power_iteration(double A[][MAX_N], int n, double *eigenval, double eigenvec[]) {\n"
        "    double v[MAX_N], v_new[MAX_N];\n"
        "    double norm, lambda_old=0, lambda_new;\n"
        "    int i,j,iter;\n"
        "    for(i=0;i<n;i++) v[i]=1.0;\n"
        "    for(iter=0;iter<1000;iter++){\n"
        "        for(i=0;i<n;i++){\n"
        "            v_new[i]=0;\n"
        "            for(j=0;j<n;j++) v_new[i]+=A[i][j]*v[j];\n"
        "        }\n"
        "        norm=0;\n"
        "        for(i=0;i<n;i++) norm+=v_new[i]*v_new[i];\n"
        "        norm=sqrt(norm);\n"
        "        lambda_new=0;\n"
        "        for(i=0;i<n;i++){\n"
        "            v_new[i]/=norm;\n"
        "            lambda_new+=v_new[i]*v[i];\n"
        "        }\n"
        "        if(fabs(lambda_new-lambda_old)<1e-8) break;\n"
        "        lambda_old=lambda_new;\n"
        "        for(i=0;i<n;i++) v[i]=v_new[i];\n"
        "    }\n"
        "    lambda_new=0;\n"
        "    for(i=0;i<n;i++){\n"
        "        double Av=0;\n"
        "        for(j=0;j<n;j++) Av+=A[i][j]*v[j];\n"  /* FIXED: was v[i], now v[j] */
        "        lambda_new+=Av*v[i];\n"
        "    }\n"
        "    *eigenval=lambda_new;\n"
        "    for(i=0;i<n;i++) eigenvec[i]=v[i];\n"
        "}\n\n"
        "int main(void){\n"
        "int i,j,k;\n"
        "double trace_val, det_val;\n"); 
}
stmts
{
    fprintf(out, "return 0;\n}\n");
    fclose(out);
    printf("[Compiler] Generating output.c\n");
    system("gcc output.c -lm -o output.exe");
    printf("[Compiler] Running output.exe\n\n");
    system("output.exe");
}
;

stmts:
      stmts stmt
    | stmt
;

stmt:
      matrix_decl
    | assign_stmt
    | expr ';' { 
        fprintf(out, "%s", $1.code);
        free($1.code);
        free($1.result_var);
    }
;

matrix_decl:
    MATRIX ID '(' number ',' number ')' '=' '[' num_list ']' ';'
    {
        int r = (int)$4;
        int c = (int)$6;
        int i, j, idx=0;

        if (num_count != r * c) {
            printf("Error: Matrix %s expects %d elements, got %d\n", $2, r*c, num_count);
            exit(1);
        }

        add_matrix($2, r, c);
        
        fprintf(out, "double %s[%d][%d] = {\n", $2, r, c);
        for(i=0; i<r; i++) {
            fprintf(out, "  {");
            for(j=0; j<c; j++) {
                fprintf(out, "%g", num_buffer[idx++]);
                if(j < c-1) fprintf(out, ",");
            }
            fprintf(out, "}");
            if(i < r-1) fprintf(out, ",");
            fprintf(out, "\n");
        }
        fprintf(out, "};\n");

        fprintf(out,
            "printf(\"%s = \\n\");\n"
            "for(i=0;i<%d;i++){\n"
            "  for(j=0;j<%d;j++) printf(\"%%g \", %s[i][j]);\n"
            "  printf(\"\\n\");\n"
            "}\n"
            "printf(\"\\n\");\n\n",
            $2, r, c, $2);
            
        free($2);
        num_count = 0;
    }
;

num_list:
      number { num_count = 0; num_buffer[num_count++] = $1; }
    | num_list ',' number { num_buffer[num_count++] = $3; }
;

assign_stmt:
    ID '=' expr ';'
    {
        add_matrix($1, $3.rows, $3.cols);
        
        fprintf(out, "%s", $3.code);
        fprintf(out, "double %s[%d][%d];\n", $1, $3.rows, $3.cols);
        
        fprintf(out,
            "for(i=0;i<%d;i++){\n"
            "  for(j=0;j<%d;j++) %s[i][j] = %s[i][j];\n"
            "}\n",
            $3.rows, $3.cols, $1, $3.result_var);
            
        fprintf(out,
            "printf(\"%s = \\n\");\n"
            "for(i=0;i<%d;i++){\n"
            "  for(j=0;j<%d;j++) printf(\"%%g \", %s[i][j]);\n"
            "  printf(\"\\n\");\n"
            "}\n"
            "printf(\"\\n\");\n\n",
            $1, $3.rows, $3.cols, $1);

        free($3.code);
        free($3.result_var);
        free($1);
    }
;

expr:
    ID 
    {
        Matrix* m = find_matrix($1);
        if(!m) { printf("Error: Undefined matrix %s\n", $1); exit(1); }
        $$.code = strdup(""); 
        $$.result_var = strdup($1);
        $$.rows = m->rows;
        $$.cols = m->cols;
        free($1);
    }
  | '(' expr ')'
    {
        $$.code = $2.code;
        $$.result_var = $2.result_var;
        $$.rows = $2.rows;
        $$.cols = $2.cols;
    }
  | expr '+' expr
    {
        if ($1.rows != $3.rows || $1.cols != $3.cols) {
            printf("Error: Dimension mismatch in addition (%dx%d + %dx%d)\n", 
                   $1.rows, $1.cols, $3.rows, $3.cols);
            exit(1);
        }
        char* temp = new_temp();
        char* code = gen_binary_op($1.code, $3.code, temp, $1.rows, $1.cols, 
                                   $1.result_var, $3.result_var, "Addition", "+");
        
        $$.code = code;
        $$.result_var = temp;
        $$.rows = $1.rows;
        $$.cols = $1.cols;
        
        /* Free intermediate code buffers */
        free($1.code);
        free($3.code);
        free($1.result_var);
        free($3.result_var);
    }
  | expr '-' expr
    {
        if ($1.rows != $3.rows || $1.cols != $3.cols) {
            printf("Error: Dimension mismatch in subtraction (%dx%d - %dx%d)\n",
                   $1.rows, $1.cols, $3.rows, $3.cols);
            exit(1);
        }
        char* temp = new_temp();
        char* code = gen_binary_op($1.code, $3.code, temp, $1.rows, $1.cols,
                                   $1.result_var, $3.result_var, "Subtraction", "-");
        
        $$.code = code;
        $$.result_var = temp;
        $$.rows = $1.rows;
        $$.cols = $1.cols;
        
        free($1.code);
        free($3.code);
        free($1.result_var);
        free($3.result_var);
    }
  | expr '*' expr
    {
        if ($1.cols != $3.rows) {
            printf("Error: Dimension mismatch in multiplication (%dx%d * %dx%d)\n",
                   $1.rows, $1.cols, $3.rows, $3.cols);
            exit(1);
        }
        char* temp = new_temp();
        char* code = alloc_code_buffer();
        int R = $1.rows;
        int C = $3.cols;
        int K = $1.cols;
        
        snprintf(code, CODE_BUFFER_SIZE,
            "%s%s"
            "double %s[%d][%d];\n"
            "/* Multiplication */\n"
            "for(i=0;i<%d;i++){\n"
            "  for(j=0;j<%d;j++){\n"
            "    %s[i][j]=0;\n"
            "    for(k=0;k<%d;k++){\n"
            "      %s[i][j]+=%s[i][k]*%s[k][j];\n"
            "    }\n"
            "  }\n"
            "}\n",
            $1.code, $3.code, temp, R, C, R, C, temp, K, temp, $1.result_var, $3.result_var);
            
        $$.code = code;
        $$.result_var = temp;
        $$.rows = R;
        $$.cols = C;
        
        free($1.code);
        free($3.code);
        free($1.result_var);
        free($3.result_var);
    }
  | DET '(' expr ')'
    {
        if ($3.rows != $3.cols) {
            printf("Error: Determinant requires square matrix.\n"); 
            exit(1);
        }
        char* code = alloc_code_buffer();
        char* temp = new_temp();
        int n = $3.rows;
        
        snprintf(code, CODE_BUFFER_SIZE,
            "%s"
            "{\n"
            "  double %s_copy[MAX_N][MAX_N];\n"
            "  for(i=0;i<%d;i++)\n"
            "    for(j=0;j<%d;j++) %s_copy[i][j]=%s[i][j];\n"
            "  det_val = det_nxn(%s_copy, %d);\n"
            "  printf(\"det(%s) = %%g\\n\\n\", det_val);\n"
            "}\n",
            $3.code, temp, n, n, temp, $3.result_var, temp, n, $3.result_var);
            
        $$.code = code;
        $$.result_var = strdup("0"); 
        $$.rows = 1; 
        $$.cols = 1;
        
        free($3.code);
        free($3.result_var);
        free(temp);
    }
  | TRANSPOSE '(' expr ')'
    {
        int R = $3.cols;
        int C = $3.rows;
        char* temp = new_temp();
        char* code = alloc_code_buffer();
        
        snprintf(code, CODE_BUFFER_SIZE,
            "%s"
            "double %s[%d][%d];\n"
            "/* Transpose */\n"
            "for(i=0;i<%d;i++){\n"
            "  for(j=0;j<%d;j++){\n"
            "    %s[i][j] = %s[j][i];\n"
            "  }\n"
            "}\n"
            "printf(\"transpose(%s) = \\n\");\n"
            "for(i=0;i<%d;i++){\n"
            "  for(j=0;j<%d;j++) printf(\"%%g \", %s[i][j]);\n"
            "  printf(\"\\n\");\n"
            "}\n"
            "printf(\"\\n\");\n\n",
            $3.code, temp, R, C, R, C, temp, $3.result_var, $3.result_var, R, C, temp);
            
        $$.code = code;
        $$.result_var = temp;
        $$.rows = R;
        $$.cols = C;
        
        free($3.code);
        free($3.result_var);
    }
  | TRACE '(' expr ')'
    {
        if ($3.rows != $3.cols) { 
            printf("Error: Trace requires square matrix.\n"); 
            exit(1); 
        }
        char* code = alloc_code_buffer();
        snprintf(code, CODE_BUFFER_SIZE,
            "%s"
            "trace_val = 0;\n"
            "for(i=0; i<%d; i++) trace_val += %s[i][i];\n"
            "printf(\"trace(%s) = %%g\\n\\n\", trace_val);\n",
            $3.code, $3.rows, $3.result_var, $3.result_var);
            
        $$.code = code;
        $$.result_var = strdup("0");
        $$.rows = 1; 
        $$.cols = 1;
        
        free($3.code);
        free($3.result_var);
    }
  | EIGENVAL '(' expr ')'
    {
        if ($3.rows != $3.cols) {
            printf("Error: Eigenvalues require square matrix.\n"); 
            exit(1);
        }
        char* code = alloc_code_buffer();
        char* temp = new_temp();
        int n = $3.rows;
        
        snprintf(code, CODE_BUFFER_SIZE,
            "%s"
            "{\n"
            "  double %s_copy[MAX_N][MAX_N], eigenval, eigenvec[MAX_N];\n"
            "  for(i=0;i<%d;i++)\n"
            "    for(j=0;j<%d;j++) %s_copy[i][j]=%s[i][j];\n"
            "  power_iteration(%s_copy, %d, &eigenval, eigenvec);\n"
            "  printf(\"dominant eigenvalue(%s) = %%g\\n\\n\", eigenval);\n"
            "}\n",
            $3.code, temp, n, n, temp, $3.result_var, temp, n, $3.result_var);
            
        $$.code = code;
        $$.result_var = strdup("0");
        $$.rows = 1; 
        $$.cols = 1;
        
        free($3.code);
        free($3.result_var);
        free(temp);
    }
  | EIGENVEC '(' expr ')'
    {
        if ($3.rows != $3.cols) {
            printf("Error: Eigenvectors require square matrix.\n"); 
            exit(1);
        }
        char* code = alloc_code_buffer();
        char* temp = new_temp();
        int n = $3.rows;
        
        snprintf(code, CODE_BUFFER_SIZE,
            "%s"
            "{\n"
            "  double %s_copy[MAX_N][MAX_N], eigenval, eigenvec[MAX_N];\n"
            "  for(i=0;i<%d;i++)\n"
            "    for(j=0;j<%d;j++) %s_copy[i][j]=%s[i][j];\n"
            "  power_iteration(%s_copy, %d, &eigenval, eigenvec);\n"
            "  printf(\"dominant eigenvector(%s):\\n  [\");\n"
            "  for(i=0;i<%d;i++) {\n"
            "    printf(\"%%g\", eigenvec[i]);\n"
            "    if(i<%d-1) printf(\", \");\n"
            "  }\n"
            "  printf(\"]\\n\\n\");\n"
            "}\n",
            $3.code, temp, n, n, temp, $3.result_var, temp, n, $3.result_var, n, n);
            
        $$.code = code;
        $$.result_var = strdup("0");
        $$.rows = 1; 
        $$.cols = 1;
        
        free($3.code);
        free($3.result_var);
        free(temp);
    }
  | INVERSE '(' expr ')'
    {
        if ($3.rows != $3.cols) {
            printf("Error: Inverse requires square matrix.\n"); 
            exit(1);
        }
        char* temp = new_temp();
        char* code = alloc_code_buffer();
        int n = $3.rows;
        
        snprintf(code, CODE_BUFFER_SIZE,
            "%s"
            "double %s[%d][%d];\n"
            "{\n"
            "  double %s_copy[MAX_N][MAX_N], %s_inv[MAX_N][MAX_N];\n"
            "  int success;\n"
            "  for(i=0;i<%d;i++)\n"
            "    for(j=0;j<%d;j++) %s_copy[i][j]=%s[i][j];\n"
            "  success = gauss_jordan(%s_copy, %s_inv, %d);\n"
            "  if(!success) {\n"
            "    printf(\"Warning: Matrix singular, inverse undefined.\\n\");\n"
            "    for(i=0;i<%d;i++)\n"
            "      for(j=0;j<%d;j++) %s[i][j]=0;\n"
            "  } else {\n"
            "    for(i=0;i<%d;i++)\n"
            "      for(j=0;j<%d;j++) %s[i][j]=%s_inv[i][j];\n"
            "  }\n"
            "}\n",
            $3.code, temp, n, n, temp, temp, n, n, temp, $3.result_var, 
            temp, temp, n, n, n, temp, n, n, temp, temp);
            
        $$.code = code;
        $$.result_var = temp;
        $$.rows = n;
        $$.cols = n;
        
        free($3.code);
        free($3.result_var);
    }
;

number:
    NUMBER { $$ = $1; }
;

%%

void add_matrix(char* name, int r, int c) {
    Matrix* m = find_matrix(name);
    if(m) {
        m->rows = r;
        m->cols = c;
        return;
    }
    strcpy(symtab[symcount].name, name);
    symtab[symcount].rows = r;
    symtab[symcount].cols = c;
    symcount++;
}

Matrix* find_matrix(char* name) {
    int i;
    for (i = 0; i < symcount; i++) {
        if (strcmp(symtab[i].name, name) == 0)
            return &symtab[i];
    }
    return NULL;
}

int yyerror(const char* s) {
    printf("Parse error: %s\n", s);
    return 0;
}

int main(void) {
    return yyparse();
}