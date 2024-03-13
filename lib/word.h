/* lib.d file word.h */
typedef int gen; /*name for generator*/

/* |INVALID_GEN| is used quite a lot in various files. 
Because space filled with zeroes can be allocated easily, we define
this constant to be zero.
This will allow for shorter coding, but code which relies on this
shoule be protected with a conditional compilation.
*/
#define IDENTITY 				-1
#define INVALID_GEN 	0
/*This assumes that the generators are numbered
from 1 to |num_gens|*/
/* The inverse of a generator is found by the function |inv()|. Both the
generator and its inverse are assumed to have program names.
*/
#define inv(g)		inv_of[(g)]

typedef struct word {
   gen *g;
   int first;	       /* the position of the first entry in the word */
   int last;	       /* position of last entry; if word is empty, 
                          define this to be position BEFORE first entry */
   int space;	       /* this should be a power of 2 */
   char type; /* initialised to 0, 'c' for commutator, or 's' for string */
   int n; /* used if the word is a power of that pointed to by the "g" field.
            0 means the same as 1, i.e. the word is not a proper power */
} word;

/* A |word_traverser| is used to run through a word generator by
generator.
*/
typedef struct word_traverser {
	word * wp;
	int posn;
} word_traverser;

extern int num_gens;
extern gen * inv_of;
extern boolean * default_involution;
extern word * user_gen_name;
extern int gen_array_size;

extern void word_print PARMS((FILE*, word *));
extern void user_word_print  PARMS((FILE*, word * )); /*print for humans*/
extern void word_init PARMS(( word * ));
extern void pc_word_init PARMS(( word * ));
extern void word_longinit PARMS(( word *, int ));
extern void word_clear PARMS(( word * ));
extern void word_reset PARMS(( word * ));
extern void pc_word_reset PARMS(( word * ));

extern void word_traverser_init PARMS((word_traverser * wtp, word
* wp));
 
/*
Nothing has to be done since no space is assigned to a |word_traverser|.
*/
#define word_traverser_clear(X) 

/* The signature of an ordering of two words.
Are the words pointed at by |wp1| and |wp2| in the order we would wish to
see them in in an ordered list (or priority queue), i.e. is |*wp1| the "smaller"
(shorter, or lexicographically less) of the two? The function returns 1 if
|*wp1| is earlier in the ordering than |*wp2|, -1 if the opposite is true
(|*wp2| is earlier than |*wp1|, 0 if the two words are equal).
We consider this function to be a kind of signature on the ordering of the
two words, hence the name.
This function expects pointers to words as its arguments. The same function,
but expecting datapointer arguments, is in the function table for lists of
words. 
*/
#define word_sgn(wp1,wp2)  word_sgn_dp((dp)(wp1), (dp)(wp2))

/* Copy the word pointed at by |oldp| into the address
specified by |newp|. Before the application of the function |newp| can point 
either to a freshly initialized (empty) word or to an actual word of any
length.
This function expects pointers to words as its arguments. The same function,
but expecting datapointer arguments, is in the function table for lists of
words.
*/
#define word_cpy(oldp,newp)  word_cpy_dp((dp) (oldp), (dp) (newp))
/*
Move, rather than copy, the word pointed at address oldwp to address newwp.
Afterwards old wp points to a nonitialised word. This function is much quicker> is equivalent in usage to the sequence word_cpy(old,new); word_clear(old);
but more efficient */
extern void word_mv PARMS((word * oldwp, word * newwp));


/* |word_next()| is used when traversing a word, generator by generator. It
returns |TRUE| when there is another generator, and sets |*gp| equal
to that generator.
The procedure does not change the word being traversed.
*/
extern boolean word_next PARMS(( word_traverser * wtp, gen *g));

/* Reduce a word.
Before the function is applied |*out| could be a freshly initialized word or
an actual word already in use. Space is assigned or reassigned as necessary
within this function.
*/
extern boolean word_reduce PARMS (( word * in, word * out ));

/* Cyclically reduce a word.
Before the function is applied |*out| could be a freshly initialized word or
an actual word already in use. Space is assigned or reassigned as necessary
within this function.
*/
extern boolean word_creduce PARMS (( word * in, word * out ));

/* Construct the inverse of a word.
Before the function is applied |*inverse| could be a freshly initialized word or
an actual word already in use. Space is assigned or reassigned as necessary
within this function.
*/

extern void word_inv PARMS(( word * given, word * inverse));
extern boolean word_eqinv PARMS((word * w1p, word * w2p));

/* The length of a word.
*/
#define word_length(wp)		(((wp)->last) + 1 - ((wp)->first))

/* Adding on a generator to the right hand end of a word (right
multiplication).
*/
extern void word_put_last PARMS ((word* wp,gen g)); 

/* Adding on a generator to the left hand end of a word (left
multiplication).
*/
extern void word_put_first PARMS ((word* wp,gen g)); 

/* Take the first generator in a word, and delete it from that word.
The function returns |FALSE| if the word is trivial, |TRUE|
otherwise.
*/
extern boolean word_delget_first PARMS ((word* wp, gen *gp));

/* Delete the first generator in a word.
The function returns |FALSE| if the word is trivial, |TRUE|
otherwise.
*/
extern boolean word_del_first PARMS ((word*  wp));

/* Get the first generator  in a word, without changing the word.
The function returns |TRUE| provided the word is non-trivial, |FALSE| 
otherwise.
The  pointer |gp| is set to point to the first generator in the word pointed
to by |wp|.
*/
extern boolean word_get_first PARMS ((word* wp, gen* gp)); 

/* Get the last generator  in a word, without changing the word.
The function returns |TRUE| provided the word is non-trivial, |FALSE| 
otherwise.
The  pointer |gp| is set to point to the last generator in the word pointed
to by |wp|.
*/
extern boolean word_get_last PARMS ((word* wp, gen* gp)); 


/* Delete the last generator in a word.
The function returns |FALSE| if the word is trivial, |TRUE|
otherwise.
*/
extern boolean word_del_last PARMS ((word*  wp)); 

/* Get the last generator  in a word. and delete it at the same time.
The function returns |TRUE| provided the word is non-trivial, |FALSE| 
otherwise.
The  pointer |gp| is set to point to the last generator in the word pointed
to by |wp|.
*/
extern boolean word_delget_last PARMS ((word* wp, gen* gp)); 


/* Concatenate two words. 
|wp1|, |wp2| and |wp3| are pointers to initialized words. |*wp1| and |*wp2|
are not changed by this procedure (provided that they are both distinct from
|*wp3|). At the end |*wp3| contains the concatenation of |*wp1| and |*wp2|.
*/
extern void word_concat PARMS ((word* wp1, word* wp2, word* wp3)); 
/* without changing |*wp2| put a copy of it onto the end of |*wp1| */
extern void word_append PARMS ((word * wp1, word * wp2));
/* without changing |*wp2| put a copy of its inverse onto the end of |*wp1| */
extern void word_invappend PARMS ((word * wp1, word * wp2));
extern void word_insert PARMS ((word * wp, word * wwp));
/* without changing |*wwp| put a copy of it onto the beginning of |*wp| */

extern boolean word_empty PARMS (( word * wp));


extern void gen_print PARMS((FILE*,gen g )); /*print in user notation*/ 
extern void gen_string PARMS((char * string, gen g));

/* Function table functions for lists of generators.
(See ``function tables'' in the index.)
*/
extern int gen_sgn_dp PARMS(( dp gp1, dp gp2 )); /*signature*/
/* Returns +1 if |*gp1|<|*gp2|, -1 if |*gp1|>|*gp2|, 0 if the two generators
are equal */ 
#define gen_sgn		gen_sgn_dp
extern void gen_cpy_dp PARMS (( dp oldp, dp newp)); /*copy*/
extern gen * gen_create(); /*constructor*/ 
extern void Free_dp PARMS(( dp )); /*destructor*/ 
extern void gen_print_dp PARMS((FILE*, dp )); /*print for humans*/ 

extern void gen_print PARMS((FILE* wfile,gen g )); /*print in user notation*/

/* Function table for manipulating words.
See ``function tables'' in the index.
Each generator has two names, a name used by the user (an alphabet
character) and a name used by the program ( a small integer). So there are
two slightly different kinds of lists of words depending on which generator
names are being used. The signature functions and the print functions in the
function tables vary according to the type of list being used.
*/
extern elt_fntab WORD_fntab;
#define WORD (&WORD_fntab)
extern elt_fntab GENWT_WORD_fntab;
#define GENWT_WORD (&GENWT_WORD_fntab)


/* Function table functions for lists of words. See ``function tables'' in the
index. 
*/
extern int word_sgn_dp   PARMS((dp w1p, dp w2p));/*signature*/
extern int genwt_word_sgn_dp   PARMS((dp w1p, dp w2p));/*signature*/
/*returns +1 if first word is smaller, -1 if bigger, 0 if equal*/ 
extern void word_cpy_dp  PARMS(( dp oldp, dp newp)); /*copy*/
extern word * word_create(); /* return initialized space for a |word| */
extern void word_kill_dp  PARMS (( dp ));
extern void word_print_dp   PARMS((FILE *, dp )); /*print for humans*/

extern boolean read_next_word PARMS((word * wp,FILE * file));
extern boolean read_next_gen PARMS((word * wp,FILE * file));
extern void read_gen_name_array PARMS ((FILE * file));
extern void define_next_gen PARMS((word * wp));
extern void delete_gen PARMS((gen g));
extern void read_inverse_array PARMS ((FILE * file));
extern void default_inverse_array PARMS ((VOID));
extern boolean read_next_rel PARMS((word * relp, FILE * file));
/* relator *relp is input (read from relator, relation or rewrite rule,
stored as a relator with program notation for generators */
extern void word2prog_gen PARMS((word * wp, gen * gp));
extern void word2prog_word PARMS((word * user_wp, word * prog_wp));
extern void gen2prog_gen PARMS((gen g, gen * gp));
extern void gen2user_name PARMS((gen g, word * wp));
extern void word2user_name PARMS ((word * prog_wp, word * user_wp));
extern void word_factor PARMS((word *wp,word *wwp,int *ep));
extern char * word2string PARMS((word * wp));
