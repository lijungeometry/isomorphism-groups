/*lib.d file input.h */
extern boolean find_keyword PARMS((char * label,  FILE * file));
extern read_char PARMS((FILE * rfile));
extern void find_char PARMS((char c,FILE * rfile));
extern boolean read_next_int PARMS((int * np, FILE * file));
extern boolean read_next_float PARMS((float * fp,FILE * file));
extern boolean read_next_letter PARMS((int * lp,FILE * file));
extern boolean read_next_string PARMS((char * cp,int n,FILE * file));
extern void bad_data PARMS((VOID));
extern void format_echocheck PARMS((char * number,FILE * rfile, FILE * wfile));
extern void format_check PARMS((char * number,FILE * rfile));

