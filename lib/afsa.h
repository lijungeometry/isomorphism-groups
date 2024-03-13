/* lib.d file afsa.h */
typedef struct {
	int states;
	int symbols;
	int base_symbols;
	boolean bfs;
	boolean min; 
	int variables; 
	boolean eos;
	int ** array;
}
afsa; 

typedef struct {
	int states;
	int symbols;
	int base_symbols;
	boolean bfs;
	boolean min;
	boolean eos;
	int *** array;
}
twoafsa; 


extern void afsa_init PARMS((afsa*));
extern void afsa_clear PARMS((afsa *));
extern afsa * afsa_read PARMS((word **,FILE *));
extern void afsa_print PARMS((FILE*,afsa *));
extern twoafsa * twoafsa_read PARMS((word **,FILE *));
extern void twoafsa_print PARMS((FILE *,twoafsa *));
extern afsa * afsa_eosdelete PARMS((afsa * fsap));
extern afsa * afsa_eosadd PARMS((afsa * fsap));
extern int afsa_num_negstates PARMS((afsa * fsap));
extern FILE * wfile;

#define ACCEPTSTATE 1
#define NONACCEPTSTATE 0
