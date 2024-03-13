/* ========================== C MeatAxe =============================
   chop.c - Chop module (Calculate composition series)

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: chop.c,v 2.16 1994/03/26 06:34:04 mringe Exp $
 *
 * $Log: chop.c,v $
 * Revision 2.16  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.15  1994/03/02  11:06:14  mringe
 * Vermeide mehrfache Berechnung von Nullitaeten.
 *
 * Revision 2.14  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.13  1994/01/28  08:24:38  mringe
 * Schreibe die Erzeuger sofort raus, nicht erst am Ende.
 *
 * Revision 2.12  1993/12/14  22:38:12  mringe
 * Neu: Funktion gcd().
 *
 * Revision 2.11  1993/12/13  08:27:53  mringe
 * split(): Reihenfolge der Arguimente.
 *
 * Revision 2.10  1993/12/13  08:25:53  mringe
 * Reihenfolge der Fkt.-argumente vereinheitlicht.
 *
 * Revision 2.9  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.8  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.7  1993/12/07  16:44:24  mringe
 * Schreibe erst .cfinfo, dann die irreduziblen.
 * (Sonst funktioniert cfname() nicht richtig).
 *
 * Revision 2.6  1993/10/27  09:34:47  mringe
 * IdWords.
 *
 * Revision 2.5  1993/10/22  16:50:49  mringe
 * *** empty log message ***
 *
 * Revision 2.4  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.3  1993/10/21  08:37:53  mringe
 * Compiler warnings.
 *
 * Revision 2.2  1993/10/20  18:26:42  mringe
 * Bug bei matload() behoben.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.44  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.43  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.42  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.41  1993/09/20  19:47:59  mringe
 * *** empty log message ***
 *
 * Revision 1.40  1993/09/20  17:40:14  mringe
 * msg_prefix gesetzt.
 *
 * Revision 1.39  1993/08/27  15:27:26  mringe
 * Option -T
 *
 * Revision 1.38  1993/08/25  15:57:32  mringe
 * Verwaltung der Worte verbessert.
 *
 * Revision 1.37  1993/08/10  15:52:06  mringe
 * Test nullword(n) statt nullword(n->parent).
 *
 * Revision 1.36  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.35  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.34  1993/07/28  06:06:12  mringe
 * Probiere erst alle Woerter, die schon einmal zum Nachweis der
 * Irreduzibilitaet benutzt wurden (nicht nur bei passender Dimension).
 *
 * Revision 1.33  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.32  1993/07/19  15:00:29  mringe
 * Option -G (GAP output).
 *
 * Revision 1.31  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.30  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.29  1993/02/15  13:49:24  mringe
 * Funktionen aus cyclic.c -> yyy-Lib verschoben.
 *
 * Revision 1.28  1993/02/12  17:06:59  mringe
 * Woerter mit N Erzeugern.
 *
 * Revision 1.27  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.26  1993/02/10  19:23:03  mringe
 * !comp
 *
 * Revision 1.25  1993/01/28  07:35:51  mringe
 * Benutze "utils" Bibliothek
 *
 * Revision 1.24  93/01/19  07:21:10  07:21:10  mringe (  Michael Ringe)
 * Bug in checkspl() entfernt.
 * 
 * Revision 1.23  1993/01/15  14:50:01  mringe
 * ngen in nach files.h ausgelagert.
 *
 * Revision 1.22  1993/01/15  07:30:49  mringe
 * Erweiterung auf N Erzeuger.
 *
 * Revision 1.21  1993/01/09  14:50:26  mringe
 * `Save vectors' implementiert.
 *
 * Revision 1.20  1993/01/08  19:50:08  mringe
 * Benutze auch den Fingerprint bei der Suche nach einem
 * Wort mit minimaler Nullitaet.
 *
 * Revision 1.19  1993/01/08  10:53:51  mringe
 * Suchalgorithmus verbessert: (1) beruechsichtige die schon
 * gefundenen irreduziblen und (2) benutze die Worte, mit denen
 * vorher schon einal gesplittet wurde.
 *
 * Revision 1.18  1993/01/07  20:12:44  mringe
 * Teste zuerst, ob es einer der schon bekannten Irred. ist.
 *
 * Revision 1.17  1993/01/06  20:59:57  mringe
 * getopt in zgetopt() umbenannt.
 *
 * Revision 1.16  1992/10/13  18:00:59  mringe
 * No Changes
 *
 * Revision 1.15  1992/10/12  10:40:30  hiss
 * Neu: Option -s
 *
 * Revision 1.14  1992/10/02  16:40:00  mringe
 * Merke nullit"at=0 auch bei W"ortern -1..-6.
 *
 * Revision 1.13  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.12  1992/09/01  13:01:11  mringe
 * Bug in chop() beseitigt.
 *
 * Revision 1.11  1992/08/31  14:40:32  mringe
 * Benutze Worte -1..-6 zum Splitten.
 *
 * Revision 1.10  1992/07/29  08:27:07  mringe
 * Debug messages in extendbasis() entfernt.
 *
 * Revision 1.9  1992/07/28  08:33:55  mringe
 * Cleaned up checkspl() and findw1().
 *
 * Revision 1.8  1992/07/28  08:20:29  mringe
 * Option -m erweitert (Startwert.Maximum)
 *
 * Revision 1.7  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.6  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.5  1992/07/08  15:28:31  mringe
 * Erkennung des Zerf.k"orpers verbessert.
 *
 * Revision 1.4  1992/07/02  15:42:07  mringe
 * Sichere Erkennung des Zerf"allungsk"orpers.
 *
 * Revision 1.3  1992/07/01  14:43:09  mringe
 * Einige type casts eingef"ugt.
 *
 * Revision 1.2  1992/06/01  08:33:30  mringe
 * CL warnings entfernt.
 *
 * Revision 1.1  1992/05/23  16:46:00  mringe
 * Initial revision
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"
#include "lattice.h"
#include "files.h"
#include "split.h"
#include "words.h"


/* ------------------------------------------------------------------
   Some constants
   ------------------------------------------------------------------ */

#define MAXNULL_INIT 3
#define MAXNULL_MIN 1
#define MAXNULL_MAX 8
#define NN_LIM 15
#define MAXTRIES 50000	/* Number of tries before findw1() fails */


/* ------------------------------------------------------------------
   typedefs
   ------------------------------------------------------------------ */

typedef struct nodestruct
	{	struct nodestruct *sub, *quot, *parent;
		long dim, num;
		long *nullwords;
	}
	node_t;


typedef struct
	{	long dim, num;		/* Dimension and number */
		matrix_t *gen[MAXGEN];	/* Generators */
		long mult;		/* Multiplicity */
		long spl;		/* Degree of splitting field */
		long idword;		/* Word used for std basis */
		long fprint[MAXFP];	/* Fingerprint */
	}
	irred_t;


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void init _PL((void));
static void addnullword _PL((long w, node_t *n));
static int isanullword _PL((long w,node_t *n));
static void writetree _PL((node_t *n));
static void writeinfo _PL((node_t *root));
static void splitnode _PL((node_t *n, split_t *spl, int tr));
static node_t *chop _PL((matrix_t *gen[], node_t *parent,
    matrix_t *proj));
static matrix_t *extendbasis _PL((matrix_t *basis, matrix_t *space));
static int checkspl _PL((basis_t *basis, matrix_t *gen[], long w));
static long findidword _PL((basis_t *basis, matrix_t *gen[],
    long ggt,node_t *parent));
static int equiv _PL((basis_t *basis, matrix_t *gen[], int cf));
static void stdbasis _PL((basis_t *basis, matrix_t *gen[], int cf));
static int newirred _PL((basis_t *basis, matrix_t *gen[],
	long idword, long ggt, node_t *parent));


/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

matrix_t *generator[MAXGEN];		/* Generators */
long maxnull_init = MAXNULL_INIT;	/* Initial nullity threshold */
long maxnull_max = MAXNULL_MAX;		/* Maximal threshold */
irred_t *irred[MAXCF];			/* List of irreducibles */
long firstword = -1;
int opt_G = 0;	/* GAP output */
set_t *goodwords;			/* List of `good' words */


static char *helptext[] = {
"SYNTAX",
"    chop [<Options>] <Name>",
"",
"OPTIONS",
"    -G            GAP output (implies -Q).",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -g <NGen>     Set number of generators (default is 2).",
"    -m <N>[.<M>]  Set start value (<N>) and upper limit for the",
"                  `maxnull' parameter.",
"    -T <MaxTime>  Set CPU time limit",
"",
"FILES",
"    The input must be <NGen> square matrices. They are read from",
"    <Name>.1, <Name>.2, ... <Name>.<NGen>.",
"",
"    Output is written to various files:",
"      <Name>.cfinfo        Information about the composition factors.",
"      <Name><Dim><a>.<N>   Generators on the composition factors.",
"                           <Dim> is the dimension, <a> is a letter",
"                           discriminating different factors, and <N>",
"                           runs from 1 to <NGen>.",
NULL};

static proginfo_t pinfo =
   { "chop", "Composition Series", "$Revision: 2.16 $", helptext };




/* ------------------------------------------------------------------
   init () - Read generators
   ------------------------------------------------------------------ */

static void init()

{
    char fn[100];
    char fn0[200];
    int i;
    matrix_t *x;

    for (i = 0; i < ngen; ++i)
    {
	sprintf(fn,"%s.%d",cfbasename,i+1);
	if (i == 0) strcpy(fn0,fn);
	x = generator[i] = matload(fn);
	if (x->nor != x->noc)
	    errexit(ERR_NOTSQUARE,fn);
	if (i > 0 && x->nor != generator[0]->nor)
	    errexit(ERR_INCOMPAT,strcat(strcat(fn0," and "),fn));
    }
}


/* ------------------------------------------------------------------
   addnullword() - Add a word to the 'nullwords' list
   isanullword() - Look up a word in the 'nullwords' list
   ------------------------------------------------------------------ */

#define ALLOCBLKSIZE 100

static void addnullword(w,n)
long w;
node_t *n;

{
    int i;

    if (n->nullwords == NULL)
    {
	n->nullwords =(long *)malloc(ALLOCBLKSIZE*sizeof(long));
	n->nullwords[0] = ALLOCBLKSIZE;
	n->nullwords[1] = 0;
    }
    else if (n->nullwords[1] >= n->nullwords[0] - 2)
    {
	n->nullwords[0] += ALLOCBLKSIZE;
	n->nullwords = (long *) realloc(n->nullwords,
		(size_t)(n->nullwords[0]*sizeof(long)));
    }
    if (n->nullwords == NULL)
	FATAL("OUT OF MEMORY");
    i = (int) ++n->nullwords[1];
    n->nullwords[i+1] = w;
}



static int isanullword(w,n)
long w;
node_t *n;

{	int i;

	if (n == NULL) return 0;
	if (n->nullwords != NULL)
	{	for (i = 2; i <= (int) n->nullwords[1] + 1; ++i)
			if (n->nullwords[i] == w) return 1;
	}
	return isanullword(w,n->parent);
}



/* ------------------------------------------------------------------
   writetree() - Write the composition series to stdout
   ------------------------------------------------------------------ */

static void writetree(n)
node_t *n;

{
    static int count = 0;	/* Characters written in one line */
    static int first = 1;
    int i;

    if (n->sub == NULL)		/* Irreducible */
    {
	if (count >= 20)
	{   printf("\n");
	    count = 0;
	}
	for (i = 0; cfinfo[i].dim != n->dim || cfinfo[i].num != n->num;
	    ++i);
	if (opt_G)
	{
	    if (!first) printf(","); else { printf(" "); first = 0; }
	    printf("%d",i+1);
	}
	else
	    printf("%s ",cfname(i));
	++count;
    }
    else
    {
	writetree(n->sub);
	writetree(n->quot);
    }
}


/* ------------------------------------------------------------------
   writeinfo() - Write some information to stdout and create the
	.cfinfo file.
   ------------------------------------------------------------------ */

static void writeinfo(root)
node_t *root;

{
    int i, k;

    MESSAGE(0,(
    "\n\nChopping completed: %d different composition factors\n",ncf));

    /* Write the cfinfo file
       --------------------- */
    MESSAGE(0,("Writing %s.cfinfo\n",cfbasename));
    for (i = 0; i < ncf; ++i)
    {
    	cfinfo[i].dim = irred[i]->dim;
	cfinfo[i].num = irred[i]->num;
	cfinfo[i].mult = irred[i]->mult;
	cfinfo[i].spl = irred[i]->spl;
	cfinfo[i].idword = irred[i]->idword;
    }
    writecfinfo();


    /* Write composition factors
       ------------------------- */
    if (opt_G)
    {
	printf("MeatAxe.CompositionFactors := [\n");
    	for (i = 0; i < ncf; ++i)
    	{
	    printf("  [ \"%s\", %ld, %ld ]",cfname(i),
	      irred[i]->mult,irred[i]->spl);
	    if (i < ncf-1) printf(",");
	    printf("\n");
	}
	printf("\n];\n");
    }
    else
    {
    	printf("\nName   Mult  SF  Fingerprint\n");
    	for (i = 0; i < ncf; ++i)
    	{
	    printf("%-6s %4ld  %2ld  ",cfname(i),
	      irred[i]->mult,irred[i]->spl);
	    for (k = 0; k < MAXFP; ++k)
		printf("%ld%s",irred[i]->fprint[k],k==MAXFP-1?"\n":",");
	}
    }


	/* Write the composition series
	   ---------------------------- */
    if (opt_G)
	printf("MeatAxe.CompositionSeries := [\n");
    else
	printf("\nAscending composition series:\n");
    writetree(root);
    if (opt_G)
	printf("];\n");
    else
	printf("\n");
}


/* ------------------------------------------------------------------
   splitnode() - Split a node
   ------------------------------------------------------------------ */

static void splitnode(n,spl,tr)
node_t *n;		/* Node to split */
split_t *spl;		/* Result of previous split */
int tr;			/* Indicates that it was a `dual split' */

{	


	/* If it was a dual split, subspace and quotient
	   have been calculated in the dual module. To split
	   the original module, transpose again and exchange
	   sub and quot.
	   --------------------------------------------------- */
	if (tr)
	{	int i;
		matrix_t *x, *y;

		for (i = 0; i < ngen; ++i)
		{	x = mattr(spl->sub[i]);
			matfree(spl->sub[i]);
			y = mattr(spl->quot[i]);
			matfree(spl->quot[i]);
			spl->sub[i] = y;
			spl->quot[i] = x;
		}
	}

	/* Chop the subspace and quotient
	   ------------------------------ */
	n->sub = chop(spl->sub,n,NULL);
	n->quot = chop(spl->quot,n, tr ? NULL : spl->proj);

	/* Clean up
	   -------- */
	splfree(spl);
}




/* ------------------------------------------------------------------
   chop() - Chop a module
   ------------------------------------------------------------------ */

static node_t *chop(gen, parent, proj)
matrix_t *gen[];		/* Generators */
node_t *parent;			/* Parent node */
matrix_t *proj;			/* Vectors to try first */

{
    node_t *n;
    basis_t basis;		/* Used by the word generator */
    int nn = 0;
    long maxnull;
    enum {GOOD,FP,STD} status;
    int iset;
    long w = 0, nul;
    long ggt = 0;
    long w1 = -1;		/* Word with nullity 1 */
    matrix_t *nsp;
    split_t *spl;
    int tflag = 0;		/* Generators have been transposed */
    matrix_t *wtr;		/* Transposed word */
    matrix_t *gentr[MAXGEN];	/* Transposed generators */

    MESSAGE(0,("Chop: Dim=%ld\n",gen[0]->nor));

    /* Allocate a new node
       ------------------- */
    n = (node_t *)malloc(sizeof(node_t));
    n->parent = parent;
    n->dim = gen[0]->nor;
    n->sub = n->quot = NULL;
    n->nullwords = NULL;

    /* Initialize the word generator
       ----------------------------- */
    initbasis(gen,ngen,&basis);

    /* Dimension 1 is simple...
      ------------------------- */
    if (n->dim == 1)
    {
	int cf = newirred(&basis,gen,(long)-1,(long)1,parent);
	n->num = irred[cf]->num;
	freebasis(&basis);
	return n;
    }

    /* Try saved vectors.
       ------------------ */
    if (proj != NULL && proj->nor > 0)
    {
	spl = split(proj,ngen,gen,0);
	if (spl->result != 0) 		/* Module was split */
	{
	    freebasis(&basis);
	    MESSAGE(0,(
		"Split (saved vectors):Subspace=%ld,Quotient=%ld\n",
		spl->sub[0]->nor,spl->quot[0]->nor));
	    splitnode(n,spl,0);
	    return n;
	}
    }


    /* Try words
       --------- */
    maxnull = maxnull_init;
    for (status = GOOD, iset = 0; ; )
    {	
	/* Make the next word, transpose and
	   calculate the null-space.
	   ------------------------------------- */
	do
	{
	    if (status == GOOD)
	    {
		if (iset >= goodwords->len)
		{    w = 0;
		    status = FP;
		}
		else
		    w = goodwords->buf[iset++];
	    }
	    if (status == FP)
	    {
		if (--w < -6)
		{    w = 0;
		    status = STD;
		}
	    }
	    if (status == STD)
	    {
		w = nextword(w);
	    }
	}
	while (isanullword(w,n) ||
	       (status != GOOD && set_contains(goodwords,w))
	     );

	mkword(&basis,w);
	wtr = mattr(basis.w);
	nsp = nullsp(basis.w);
	nul = nsp->nor;
	ggt = gcd(ggt,nul);

	/* Adjust maxnull, if too many words with
	   'large' (or 'small') nulities appear.
	   -------------------------------------- */
	if (nn >= NN_LIM)
	{   nn = 0;
	    if (maxnull < maxnull_max) ++maxnull;
	}
	else if (nn <= -NN_LIM)
	{   nn = 0;
	    if (maxnull > MAXNULL_MIN) --maxnull;
	}
	MESSAGE(0,("  Word %5ld, nul=%3ld, maxnull=%3ld, gcd=%ld\n",
	    w,nul,maxnull,ggt));


	/* Nullity 0: Discard, but remember
	   -------------------------------- */
	if (nul == 0)
	{
  	    matfree(nsp);
	    matfree(wtr);
	    addnullword(w,n);
	    continue;
	}

	/* Remember if a word with nullity 1 is found.
	   It will be used for the standard basis later.
	   --------------------------------------------- */
	if (w > 0 && nul == 1 && w1 == -1)
	{
	    w1 = w;
	    set_insert(goodwords,w1);
	}

	/* Adjust nn, depending on the difference
	   between nullity and the current maxnull value.
	   ---------------------------------------------- */
	if (nul > maxnull)
	{   nsp->nor = (long)1;
	    ++nn;
	}
	else if (nul < maxnull) --nn;

	/* Try to split
	   ------------ */
	spl = split(nsp,ngen,gen,0);
	matfree(nsp);
	if (spl->result != 0) 		/* Module was split */
	{
	    freebasis(&basis);
	    matfree(wtr);
	    MESSAGE(0,("Split: Subspace=%ld, Quotient=%ld\n",
		spl->sub[0]->nor,spl->quot[0]->nor));
	    set_insert(goodwords,w);
	    splitnode(n,spl,0);
	    break;
	}

	/* Not split, try dual split
	   ------------------------- */
	splfree(spl);
	if (nul > maxnull)	/* Sorry, nullity too large */
	{   matfree(wtr);
	    continue;
	}
	if (!tflag)		/* Transpose generators */
	{   int i;
	    for (i = 0; i < ngen; ++i)
		gentr[i] = mattr(gen[i]);
	    tflag = 1;
	}
	nsp = nullsp(wtr);
	matfree(wtr);
	if (nsp->nor != nul)	/* Should be equal... */
	    FATAL("CHOP ERROR");
	spl = split(nsp,ngen,gentr,1);
	matfree(nsp);
	if (spl->result == 0) 	/* Irreducible! */
	{
	    int cf = newirred(&basis,gen,w1,ggt,parent);
	    n->num = irred[cf]->num;
	    splfree(spl);
	    fflush(stdout);
	    freebasis(&basis);
	    set_insert(goodwords,w);
	    break;
	}
	else	/* The dual was split */
	{
	    freebasis(&basis);
	    MESSAGE(0,("Dual split: Subspace=%ld, ",spl->sub[0]->nor));
	    MESSAGE(0,("Quotient=%ld\n",spl->quot[0]->nor));
	    set_insert(goodwords,w);
	    splitnode(n,spl,1);
	    break;
	}
    }

    if (tflag)	/* Clean up */
    {
	int i;
	for (i = 0; i < ngen; ++i)
	    matfree(gentr[i]);
    }
    return n;
}



/* ------------------------------------------------------------------
   extendbasis() - Find a vector in 'space' which is not in the
	(linear) span of 'basis'. 'basis' must be linearly
	independent.
   ------------------------------------------------------------------ */

static matrix_t *extendbasis(basis, space)
matrix_t *basis, *space;

{
    long i, j, piv;
    long dimb = basis->nor;
    long dims = space->nor;
    PTR tmp, x, y;
    FEL f;


    /* Concatenate basis and space
       --------------------------- */
    zsetlen(zfl,basis->noc);
    tmp = zalloc(dimb+dims);
    memcpy(tmp,basis->d,zsize(dimb));
    x = tmp;
    zadvance(&x,dimb);
    memcpy(x,space->d,zsize(dims));

    /* Clean with basis
       ---------------- */
    for (i = 1, x = tmp; i <= dimb; zadvance(&x,(long)1), ++i)
    {
	piv = zfindpiv(x,&f);
	if (piv == 0) FATAL("extendbasis(): zero vector in basis");
	y = x;
	for (j = i+1; j <= dimb+dims; ++j)
	{
	    zadvance(&y,(long)1);
	    zaddmulrow(y,x,zsub(F_ZERO,zdiv(zextract(y,piv),f)));
	}
    }


    /* Find the first non-zero row
       --------------------------- */
    x = tmp;
    zadvance(&x,dimb);
    for (j = 1; zfindpiv(x,&f) == 0; ++j, zadvance(&x,(long)1));
    free(tmp);
    if (j > dims)
    {
	FATAL("extendbasis() failed");
    }
    return matextract(space,j,(long)1);
}



/* ------------------------------------------------------------------
   checkspl() - Test for [E:F] > 1. Returns 1 if w is a word with
	minimal nullity, i.e., nullity(w) = [E:F], where E=splitting
	field. Otherwise, the return value is 0.

	The algorithm uses that
	(i) [E:F] = dim End_A(V)    (where A is the algebra)
	(ii) End_A(V) operates regularly on ker(w)
   ------------------------------------------------------------------ */

#define MAXENDO 10	/* Max. dimension of endomorphism ring */


static int checkspl(basis,gen,w)
basis_t *basis;
matrix_t *gen[];
long w;

{
    matrix_t *nsp, *v1;		/* Null space & first seed vector */
    matrix_t *sb1, *sb2,	/* Standard bases */
	     *sbi1,*sbi2;	/* Inverse of sb1, sb2 */
    matrix_t *span = NULL;	/* Span of v1 under endo[0..nendo-1] */
    matrix_t *g1[MAXGEN];	/* Generators in standard basis sb1 */
    matrix_t *g2[MAXGEN];	/* Generators in standard basis sb2 */
    matrix_t *endo[MAXENDO];	/* Endomorphisms */
    int nendo = 0;		/* # of endomorphisms obtained so far */
    long dim = basis->b[1]->nor;
    int i, result;

    /* Calculate the word and its null-space. Take
       the first vector and change to standard basis.
       ---------------------------------------------- */
    mkword(basis,w);
    nsp = nullsp(basis->w);
    v1 = matextract(nsp,(long)1,(long)1);
    sb1 = sbasis(v1,ngen,gen);
    sbi1 = matinv(sb1);
    for (i = 0; i < ngen; ++i)
    {
	g1[i] = matdup(sb1);
	matmul(g1[i],gen[i]);
	matmul(g1[i],sbi1);
    }

    sb2 = sbi2 = NULL;	/* Mark as unused */
    while (1)
    {  
	matrix_t *v2;

	/* Spin up v1 under all endomorphisms found so far. If this
	   yields the whole null-space, we know that the endomorphism
	   ring has at least dimension dim ker (w).
	   ---------------------------------------------------------- */
	if (span != NULL) matfree(span);
	span = matspin(v1,nendo,endo);
	if (span->nor == nsp->nor)
	{
	    result = 1;	/* Successfull! */
	    break;
	}

	/* Take a vector which is not in span(v1)
	   and make the standard basis as with v1.
	   --------------------------------------- */
	if (sb2 != NULL)
	{
	    int j;
	    matfree(sb2); matfree(sbi2);
	    for (j = 0; j < ngen; ++j) matfree(g2[j]);
	}
	v2 = extendbasis(span,nsp);

	sb2 = sbasis(v2,ngen,gen);
	matfree(v2);
	sbi2 = matinv(sb2);
	for (i = 0; i < ngen; ++i)
	{
	    g2[i] = matdup(sb2);
	    matmul(g2[i],gen[i]);
	    matmul(g2[i],sbi2);
	}

	/* Compare the two representations. If they are identical,
	   the matrix which transforms one into the other is an
	   endomorphism. If they are different, we know that the
	   splitting field degree must be smaller than dim ker (w).
	   -------------------------------------------------------- */
	result = 1;
	for (i = 0; result && i < ngen; ++i)
	    if (memcmp(g1[i]->d,g2[i]->d,zsize(dim)))
	        result = 0;
	if (result == 0) break;	/* Not successfull */

	/* We have found an endomorphism
	   ----------------------------- */
	if (nendo >= MAXENDO) FATAL("Too many endomorphisms");
	endo[nendo] = matdup(sbi2);
	matmul(endo[nendo],sb1);
	++nendo;

    }

    /* Clean up
       -------- */
    matfree(nsp);
    matfree(span);
    matfree(sb1); matfree(sbi1);
    for (i = 0; i < ngen; ++i) matfree(g1[i]);
    if (sb2 != NULL)
    {
	matfree(sb2); matfree(sbi2);
        for (i = 0; i < ngen; ++i) matfree(g2[i]);
    }
    while (nendo > 0)
	matfree(endo[--nendo]);

    return result;
}



/* ------------------------------------------------------------------
   findidword() - Try to find a word with nullity = degree of splitting
	field = [E:F]. Use the g.c.d. of all nullities as an upper
	bound for [E:F].
   ------------------------------------------------------------------ */

static long findidword(basis,gen,ggt,parent)
basis_t *basis;
matrix_t *gen[];
long ggt;
node_t *parent;

{
    long i, nul;
    long gcdiv = ggt;	/* G.c.d. of all nullities */
    long count = 0;	/* Number of words with nullity != 0 */


    /* Main loop: Try all words
       ------------------------ */
    for (i = -1; count <= MAXTRIES; i = nextword(i))
    {
	if (isanullword(i,parent)) continue;
	mkword(basis,i);
	nul = nullity(basis->w);
	if (nul == 0) continue;	/* Use only nullities > 0 */
	if (nul == 1)		/* Nullity 1: OK */
	    return i;
	++count;
	gcdiv = gcd(gcdiv,nul);
	if (count > 3 && gcdiv > 1 && nul == gcdiv)
	{
	    if (checkspl(basis,gen,i))
		return i;
	}
    }

    /* Failed...
       --------- */
    FATAL("findidword() failed");
    return 0;
}



/* ------------------------------------------------------------------
   equiv() - Check if a given module is isomorphic to one of the
	composition factors in the list.
   ------------------------------------------------------------------ */

static int equiv(basis,gen,cf)
basis_t *basis;
matrix_t *gen[];
int cf;

{
    matrix_t *b, *bi, *seed, *seeds, *sv, *g;
    int result = 0, j;
    long i;
    size_t datasize;

    mkword(basis,irred[cf]->idword);
    seed = nullsp(basis->w);
    if (seed->nor != irred[cf]->spl)	/* Nullities don't match */
    {	matfree(seed);
	return 0;
    }

    /* Try each seed vector. This is because we don't know
       which seed vector has been used to bring module #cf
       into standard form.
       --------------------------------------------------- */
    seeds = mkseed(seed);
    datasize = zsize(irred[cf]->dim);
    for (i = 1; result == 0 && i <= seeds->nor; ++i)
    {	
	sv = matextract(seeds,i,(long)1);
	b = sbasis(sv,ngen,gen);
	bi = matinv(b);
	for (j = 0; result == 0 && j < ngen; ++j)
	{
	    g = matdup(b);
	    matmul(g,gen[j]);
	    matmul(g,bi);
	    if (memcmp(g->d,irred[cf]->gen[j]->d,datasize))
		result = 1;
	    matfree(g);
	}
	matfree(b);
	matfree(bi);
	matfree(sv);
    }
    matfree(seed);
    matfree(seeds);
    return (result == 0);
}


/* ------------------------------------------------------------------
   stdbasis() - Change to standard basis. It is assumed that basis.w
	holds an appropriate word (i.e., a word with minimal
	nullity). The generators in the new basis are stored at
	position 'cf' in the list of irreducibles.
   ------------------------------------------------------------------ */

static void stdbasis(basis,gen,cf)
basis_t *basis;
matrix_t *gen[];
int cf;

{	matrix_t *nsp, *b, *bi;
	int i;

	/* Calculate null space and spin up one
	   (randomly chosen) null vector
	   ------------------------------------ */
	nsp = nullsp(basis->w);
	irred[cf]->spl = nsp->nor;
	b = sbasis(nsp,ngen,gen);
	bi = matinv(b);

	/* Transform all generators into standard basis
	   -------------------------------------------- */
	for (i = 0; i < ngen; ++i)
	{	irred[cf]->gen[i] = matdup(b);
		matmul(irred[cf]->gen[i],gen[i]);
		matmul(irred[cf]->gen[i],bi);
	}

	/* Clean up
	   -------- */
	matfree(b);
	matfree(bi);
	matfree(nsp);
}


/* ------------------------------------------------------------------
   newirred() - Check if a given irreducible module is already
	contained in the list of composition factors. If yes, return
	the index. If not, insert the new irreducible and return its
	index.
   ------------------------------------------------------------------ */

static int newirred(basis, gen, idword, ggt, parent)
basis_t *basis;
matrix_t *gen[];
long idword;		/* A word with nullity 1, or -1 */
long ggt;		/* g.c.d. of nullities found in chop() */
node_t *parent;		/* We use the nullword list */

{
    int i, k;
    long dim, fp[MAXFP];

    dim = gen[0]->nor;
    zsetlen(zfl,dim);
    makefp(basis,fp);	/* Calculate fingerprint */

    /* Check if the module is already in the list
	   ------------------------------------------ */
    for (i = 0; i < ncf && dim >= irred[i]->dim; ++i)
    {
	/* Compare dimensions and fingerprints
	   ----------------------------------- */
	if (dim < irred[i]->dim)
	    continue;
	if (memcmp(fp,irred[i]->fprint,sizeof(fp)))
	    continue;

	/* Check for equivalence
	   --------------------- */
	if (equiv(basis,gen,i))
	{
	    ++irred[i]->mult;
    	    cfinfo[i].dim = irred[i]->dim;  /* cfname() needs this */
    	    cfinfo[i].num = irred[i]->num;
	    MESSAGE(0,("Irreducible (%s)\n",cfname(i)));
	    return i;
	}
    }

    /* Allocate a new irred_t structure
       -------------------------------- */
    if (ncf >= MAXCF)
	FATAL("TOO MANY COMPOSITION FACTORS");
    for (k = ncf-1; k >= i; --k)
	irred[k+1] = irred[k];
    irred[i] = (irred_t *) malloc(sizeof(irred_t));
    ++ncf;

    /* Set the fields
       -------------- */
    irred[i]->dim = dim;
    irred[i]->num = (i == 0 || irred[i-1]->dim != dim) ?
	0 : irred[i-1]->num + 1;
    irred[i]->mult = 1;
    if (idword == -1)
    	irred[i]->idword = findidword(basis,gen,ggt,parent);
    else
	irred[i]->idword = idword;
    mkword(basis,irred[i]->idword);
    stdbasis(basis,gen,i);	/* gen and spl fields are set here */
    memcpy(irred[i]->fprint,fp,sizeof(fp));
    cfinfo[i].dim = irred[i]->dim;  /* cfname() needs this */
    cfinfo[i].num = irred[i]->num;
    MESSAGE(0,("Irreducible (%s",cfname(i)));
    if (irred[i]->spl > 1)
	MESSAGE(0,(", Splitting field has degree %ld",irred[i]->spl));
    MESSAGE(0,(")\n"));

    /* Write out the generators
       ------------------------ */
    for (k = 0; k < ngen; ++k)
    {
	char fn[200];
	sprintf(fn,"%s%s.%d",cfbasename,cfname(i),k+1);
   	matsave(irred[i]->gen[k],fn);
    }
    return i;
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{   node_t *root;
    int i;

    mtxinit();
    goodwords = set_alloc();

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gs:g:m:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 's':
		firstword = getint();
		if (firstword < 1 || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-s");
		break;
	    case 'g':
		ngen = getint();
		if (ngen < 2 || ngen > MAXGEN || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-g");
		break;
	    case 'm':
		maxnull_init = getint();
		if (*opt_text_ptr == '.' || *opt_text_ptr == ',')
		{	++opt_text_ptr;
			maxnull_max = getint();
		}
		if (maxnull_init < 1 || maxnull_max < maxnull_init ||
		    *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-m");
		break;
	}
    }

    if (opt_ind != argc-1) errexit(ERR_NARGS,"chop");
    setbasename(argv[argc-1]);
    MESSAGE(0,("*** CHOP MODULE ***\n\n"));
    init();
    root = chop(generator,NULL,NULL);
    writeinfo(root);
    if (msg_level >= 0)
    {
        printf("\n");
        prtimes();
    }
    return 0;
}



