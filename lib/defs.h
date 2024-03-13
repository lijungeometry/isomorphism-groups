/* lib.d file defs.h */
/*use upper case macros when there may be side-effects*/ 
#ifdef mips		/* MIPS cc rejects void * typedefs */
 typedef char * data_ptr;
 /* The previous line is what we have to use when the compiler does not
    recognize |void *|*/
#else
 typedef void * data_ptr;
#endif

typedef char boolean;
typedef data_ptr dp; /*abbreviation*/
typedef char metachar ;/* char, but can equal -1*/

#if (defined __STDC__ || defined c_plusplus)

#define PARMS(X)	X
#else
#define PARMS(X) () 
#endif

#ifdef __STDC__
#define VOID	void
#else
#define VOID
#endif

#define FALSE 0
#define TRUE 1
#ifndef NULL
#define NULL 0
#endif
#define MAXINT 2147483647 
#define INVALID_HEAP_INDEX -1

/* For historical reasons  some programs use some memory allocation
functions, some use others. They're all defined below. */
/* First those we inherited from the automata programs */

/* If space is assigned using |malloc()| or |calloc()|, we will not be
able to check, using |store_ptrs| that deallocation has been properly
carried out. All allocation should be done with one of |vzalloc1()|,
|vzalloc2()|, |valloc1()| and |valloc2()|. The zed means that all
space is initialized to zero before being handed over. The |v[z]alloc1|
procedures take a type as an argument. The |v[z]alloc2| procedure takes a
typename and an unsigned integer---the integer is the number of objects 
of that type for which one needs space.  */
#include <stdlib.h> /* Import malloc(), free(), calloc() */

extern char *malloc_value;
/*
#ifdef ultrix
extern char * calloc();
#endif
*/
extern int store_ptrs; /*used to check whether all deallocations have
				been carried out*/ 
extern void Free_dp PARMS(( dp ));
/*|Free_dp()| is the same as |free()|, except that |store_ptrs| is
updated and the argument is a |dp|*/
#define Malloc(X)	(store_ptrs++,\
			((malloc_value = malloc(X))?\
			malloc_value:\
			(char *)(void *)\
			(fprintf(stderr,"\nno more swap space\n"),exit(0),0)))
#define Calloc(X,Y)	(store_ptrs++,\
			(malloc_value = calloc((X),(Y)))?\
			malloc_value:\
			(char *)(void *)\
			(fprintf(stderr,"\nno more swap space\n"),exit(0),0))
/*|calloc| and |Calloc| set all assigned storage to zero*/
#define valloc1(A) (A *)Malloc(sizeof(A))
#define valloc2(A,N) (A *)Malloc((unsigned)(sizeof(A)*(N)))
#define vzalloc1(A) (A *)Calloc(1,(unsigned)sizeof(A))
#define vzalloc2(A,N) (A *)Calloc((unsigned)N,(unsigned)sizeof(A))

/* Next some we inherited from the Holt programs */
#define tcalloc(D,T,N) {D = (T *) calloc((unsigned)N,(unsigned)sizeof(T));\
 store_ptrs++;\
 if (D==0) { \
 fprintf(stderr,"Out of space.\n");\
 exit(0);}}
#define tmalloc(D,T,N) {D = (T *) malloc(sizeof(T)*(N)); \
 store_ptrs++;\
  if (D==0) { fprintf(stderr,"Out of space.\n"); exit(0);}}
#define tfree(D) {store_ptrs--; if (D) free( (char *) D); D=0;}

/*in these four definitions, |A| is a type name*/

/* We make |max| and |min| into routines rather than macros, to avoid
side-effects.
*/
int max PARMS ((int , int));
int min PARMS ((int, int));

/* If the compilation is done with the DEBUG flag turned on, then any
assert statement which is false leads to a message with filename, line
number and then a core dump.
*/
/*extern abort();*/
extern void error();
#ifdef DEBUG
#define assert(EX) {if ((EX)); else { \
	fprintf(stderr,"\nThe following assertion failed\n %s \nfilename: %s ",\
		"EX", __FILE__);\
	fprintf(stderr," line number %d\naborted\n", __LINE__) ;\
	fflush(stderr);\
	fflush(stdout);\
	abort();}}

#else
#define assert(EX)
#endif
