#include        <stdio.h>
/* Next 2 lines added by dfh 23/6/92 */
#include <sys/types.h>
#include <sys/times.h>



/****************************************************************************
**
*T  TypGen  . . . . . . . . . . . . . . . . . . . . . . .  type of generators
**
**  The group generators are represented using unsigned short  integers.  The
**  first generators 'a' is represented by 0, the second 'b' by 1 and so  on.
*/
typedef unsigned short  TypGen;


/****************************************************************************
**
*F  INV(<G>)  . . . . . . . . . . . . . . . . . . . . . .  invert a generator
**
**  To invert a generator we invert  the  6th  bit,  by  xor-ing  with  0x20.
*/
#define INV(G)          ((G) ^ 0x20)


/****************************************************************************
**
*V  NrGen . . . . . . . . . . . . . . . . . number of generators of the group
**
**  NrGen is the number of generators of the group,  not  counting  inverses.
*/
unsigned long   NrGen;


/****************************************************************************
**
*V  Rel[] . . . . . . . . . . . . .  array storing the relations of the group
*V  Sgr[] . . . . . . . . . . . . . . . array storing the subgroup generators
**
**  The following variables are used to store information  about  the  group.
**  Rel[]  is used to store the  relations,  Sgr[]  the  subgroup generators.
**
**  Both are a list of words, where each word starts with a size  field,  and
**  is followed by the 'characters' of the word. A word of size 0 terminates.
**  Thus  5,1,1,1,1,1,4,1,0,1,0,0,... corresponds to  b^5,(ba)^2.
*/
#define MAX_REL_LEN     4096

TypGen          Rel [MAX_REL_LEN];

TypGen          Sgr [MAX_REL_LEN];


/****************************************************************************
**
*T  TypCos  . . . . . . . . . . . . . . . . . . type used to represent cosets
*T  PtrCos  . . . . . . . . . . . . type used to represent pointers to cosets
**
**  All the cosets are simply represented  as  unsigned  integer  quantities.
**
**  If less than (<nr. of gens.> + <nr.  of nonivol.  gens.>  + 2)  / 4 MByte
**  memory is available you can't enumerate more than 2^16 coset and can save
**  memory by changing the definition of PtrCos to either (unsigned int *) or
**  (unsigned short *), depending whether int or  short is 16 bit  wide.  The
**  definition if  TypCos should   in this case be  changed  to unsigned int.
**  This can be easily achieved by defining the preprocessor symbol SMALL.
**
**  You should  never  change  TypCos  to  unsigned short  for  two  reasons:
**
**  1)  most compilers generate better code for  (unsigned int)  because  all
**      all arithmetic in C is done with integers, and using short's makes it
**      necessary to do lots of conversions between short's and int's.
**  2)  It would exhibit a bug in the  Lattice C-Compiler  on  the  Atari ST,
**      which treats unsigned short integers as signed if  indexing  into  an
**      array. E.g. if US is an unsigned short with value 65535 and  A  is an
**      array,  A[US] will be equal to  A[-1]  instead of  A[65535].  At some
**      places in the program a comment "Lattice-C Bug"  will warn you not to
**      change the code to an obvious easier expression, which also result in
**      faulty behavior.
*/
#ifndef SMALL
  typedef unsigned long   TypCos;
  typedef unsigned long   * PtrCos;
#else
  typedef unsigned int    TypCos;
  typedef unsigned short  * PtrCos;
#endif


/****************************************************************************
**
*V  Coset[][] . . . . . . . . . . . . . . . . . . . . . . . . . . coset table
**
**  Coset[g][cos] is the result of applying the generator g to the coset cos.
**  If an entry  in  this  table  is  zero  it  has  not  yet  been  defined.
**
**  The following condition is always (at least outside  Coinc())  preserved:
**  If  Coset[g][cos1] = cos2  then  Coset[INV(g)][cos2] = cos1.
*/
TypCos          Coset [256][256];


/****************************************************************************
**
** The following variables are used to store the cyclic permutations  of  the
** relations needed during the felsch strategy enumeration.
** CycGen[g] stores the index i where the list of cyclic permutations of  the
**           relations starting with  the  generator  g  start  in  CycInd[].
** CycInd[i] stores the index j where the i-th relation starts in  CycRels[].
** CycLen[i] stores the length of the i-th relation.
** CycRels[] stores for every relation  rel strings  rel rel,  rel^-1 rel^-1.
*/
unsigned short  CycGen [128];

unsigned short  CycInd [2*MAX_REL_LEN];

unsigned short  CycLen [2*MAX_REL_LEN];

TypGen          CycRels [4 * MAX_REL_LEN];


/****************************************************************************
**
*V  DedFirst  . . . . . . index of the first deduction on the deduction queue
*V  DedLast . . . . . . .  index of the last deduction on the deduction queue
*V  DedGen[i] . . . .  generator of the i-th deduction on the deduction queue
*V  DedCos[i] . . . . . .  coset of the i-th deduction on the deduction queue
**
**  The variables describe the deduction queue used in a  felsch enumeration.
*/
unsigned short  DedFirst;

unsigned short  DedLast;

TypGen          DedGen [65536];

TypCos          DedCos [65536];

short           Earlier [ 256 ];
short           EarlierLev [ 256 ];

unsigned long   NrCosets = 1;           /*     number of alive cosets.     */
unsigned long   MaxCosets = 10;         /* maximal index                   */
unsigned long   NoCoinc;
unsigned long   Perm;                   /* -p: 1 if representation wanted. */
int		MaxTime = 0;		/* time limit. added by dfh 23/6/92 */
struct 	tms 	buffer;			/* for timing. added by dfh 23/6/92 */


main ( argc, argv )
        int     argc;
        char    * argv [];
{
        char            name [64];

        while ( argc > 1 && argv[1][0] == '-' && argv[1][1] != '\0' ) {
                switch ( argv[1][1] ) {
                case 'm': ++argv; --argc; MaxCosets = atoi(argv[1]); break;
/* next line added by dfh 23/6/92 */
                case 't': ++argv; --argc; MaxTime = atoi(argv[1]); break;
                case 'p': Perm = 1; break;
                default : argc = 0; break;
                }
                ++argv; --argc;
        }

        if ( argc != 2 ) {
                fprintf(stderr,"usage:  lowindex [-m <Max>]  <Grp>\n");
                exit( 1 );
        }

        ReadGrp( Rel, argv[1] );
        InitCyc( );
        TryChoicesAt( 1, 0, 1 );
        exit( 0 );
}


/****************************************************************************
*/
#define DBG( x )
/*DBG( unsigned long           lev,  j; )*/
TryChoicesAt ( cos, gen, lev )
    TypCos              cos;
    TypGen              gen;
    long                lev;
{
    TypCos              choice;
    unsigned long       DedCurr;
    unsigned long       i;
    TypCos              c;
    TypGen              g;
    unsigned long       new;
    TypCos              lc,  rc,  hc;
    TypGen              * l,  * r;
    long                isRep, nrConj;
    TypCos              d;

DBG( lev++; )

/* Time check. Added by dfh 23/6/92 */
    if (MaxTime>0)
      { times(&buffer);
        if (buffer.tms_utime/60 > MaxTime)
        { fprintf(stderr,"Timeout.\n"); exit(0);}
      }
    DedCurr = DedLast;

    /* run over all possible choices for this entry                        */
    for ( choice=1; choice <= NrCosets+1 && choice <= MaxCosets; choice++ ) {
        if ( Coset[INV(gen)][choice] == 0 ) {

DBG( for(j=0;j<lev;j++)putchar(' ');printf("choose %d*%d = %d\n",cos,gen,choice); )

            /* enter the entry into the coset table an on the entries list */
            Coset[gen][cos] = choice;
            Coset[INV(gen)][choice] = cos;
            DedGen[DedLast] = gen;  DedCos[DedLast] = cos;
            DedFirst = DedLast;  DedLast++;
            if ( choice == NrCosets+1 ) {
                new = 1;
                NrCosets++;
            }
            else {
                new = 0;
            }

            /* make all deductions that are possible                       */
            while ( DedFirst < DedLast ) {
                for ( i = CycGen[DedGen[DedFirst]]; CycLen[i] != 0; i++ ) {

                    /* set up pointers                                     */
                    lc = DedCos[DedFirst];  l = CycRels+CycInd[i];
                    rc = lc;  r = l + CycLen[i] - 1;

                    /* scan as long as possible from the right to the left */
                    while ( l < r && (hc = Coset[INV(*r)][rc]) )  {
                        rc = hc;  r--;
                    }

                    /* scan as long as possible from the left to the right */
                    while ( l < r && (hc = Coset[*l][lc])      )  {
                        lc = hc;  l++;
                    }

                    /* if relation does close up nontrivially              */
                    if ( r <= l && Coset[*l][lc] != rc ) {

                        /* look for a deduction or a coincedence           */
                        if (      Coset[*l][lc] != 0 ) {
DBG( for(j=0;j<lev;j++)putchar(' ');printf("coinc %d = %d\n",Coset[*l][lc],rc); )
                            goto coinc;
                        }
                        else if ( Coset[INV(*r)][rc] != 0 ) {
DBG( for(j=0;j<lev;j++)putchar(' ');printf("coinc %d = %d\n",Coset[INV(*r)][rc],lc); )
                            goto coinc;
                        }
                        else {
DBG( for(j=0;j<lev;j++)putchar(' ');printf("deduce %d*%d = %d\n",lc,*l,rc); )
                            Coset[*l][lc] = rc;  Coset[INV(*r)][rc] = lc;
                            DedGen[DedLast] = *l;   DedCos[DedLast] = lc;
                            DedLast++;
                        }

                    }

                }
                DedFirst++;
            }

            /* find an empty position                                      */
            c = cos;  g = gen;
            while ( c <= NrCosets && Coset[g][c] != 0 ) {
                g = g+1; /* (g&0x20) ? INV(g)+1 : INV(g); */
                if ( g == NrGen ) {
                    g = 0;  c++;
                }
            }

            /* test if there is a lexicographically earlier renumbering    */
            isRep = 1;
            for ( d = 2; d <= NrCosets; d++ ) {
                if ( Earlier[d] == 0 ) {
                    Earlier[d] = Renumbered( d );
                    if ( Earlier[d] != 0 )
                        EarlierLev[d] = lev;
                }
                if ( Earlier[d] == -1 ) {
                    isRep = 0;
                    break;
                }
            }

            /* if there is no empty position we have found a subgroup      */
            if ( isRep && NrCosets < c ) {
                nrConj = 1;
                for ( d = 2; d <= NrCosets; d++ )
                    if ( Earlier[d] == 2 )
                        nrConj++;
                if( nrConj == NrCosets )
                    printf( "Class of index %d, normal\n",
                           NrCosets );
                else
                    printf( "Class of index %d, length %d\n",
                           NrCosets, NrCosets/nrConj );
                if ( Perm )  PrintPres();
            }

            /* otherwise we have to try all choices at this position       */
            else if ( isRep ) {
                TryChoicesAt( c, g, lev+1 );
            }

        coinc:
            for ( d = 2; d <= NrCosets; d++ ) {
                if ( EarlierLev[d] == lev ) {
                    EarlierLev[d] = 0;
                    Earlier[d] = 0;
                }
            }
            /* undo all entries again                                      */
            while ( DedCurr < DedLast ) {
                Coset[ INV(DedGen[DedLast-1]) ]
                     [ Coset[ DedGen[DedLast-1] ][ DedCos[DedLast-1] ] ] = 0;
                Coset[ DedGen[DedLast-1] ][ DedCos[DedLast-1] ] = 0;
                DedLast--;
            }
            if ( new ) {
                NrCosets--;
            }

        }

    }
DBG( lev--; )
}


/****************************************************************************
*/
Renumbered ( cos )
    TypCos              cos;
{
    TypGen              g;
    TypCos              c,  max;
    TypCos              R [256];
    TypCos              S [256];

    for ( c = 1; c <= NrCosets; c++ )
        R[c] = S[c] = 0;
    R[cos] = 1;  S[1] = cos;

    /* loop over all cosets                                                */
    c = 1;  max = 1;
    while ( c <= max ) {

        /* compare <row c> with R[ <row S[c]> ]                            */
        for ( g = 0; g != NrGen; g++ ) {

            /* if one entry is unbound, we cant say                        */
            if ( Coset[g][c] == 0 || Coset[g][S[c]] == 0 )
                return 0;

            /* if the coset in the second row has no name then             */
            if ( R[ Coset[g][S[c]] ] == 0 ) {

                if ( Coset[g][c] == max+1 ) {
                    max++;
                    R[ Coset[g][S[c]] ] = max;
                    S[max] = Coset[g][S[c]];
                }
                else if ( Coset[g][c] <= max ) {
                    return 1;
                }
                else {
                    printf("this should not happen\n");
                    exit( 1 );
                }

            }

            /* otherwise compare the entries                               */
            if ( Coset[g][c] < R[ Coset[g][S[c]] ] )
                return 1;
            else if ( R[ Coset[g][S[c]] ] < Coset[g][c] )
                return -1;
        }

        c++;
    }

    return 2;
}


/****************************************************************************
**
*F InitCyc()  . . . . . . . . initialize cyclic permutations of the relations
** Initializes the tables with the cyclic permutations of the relations  that
** are needed for the felsch strategy.
*/
InitCyc ()
{
    TypGen              * cyc,  * rel,  * r;
    TypGen              g;
    unsigned long       ind,  len;
    unsigned long       i,  l;

    /* for every relation rel store the strings rel rel, rel^-1 rel^-1     */
    cyc = CycRels;
    for ( rel = Rel; *rel != 0; rel += *rel + 1 ) {
        for ( r = rel + 1; r <= rel + *rel; ++r, ++cyc )
            *cyc = *r;
        for ( r = rel + 1; r <= rel + *rel; ++r, ++cyc )
            *cyc = *r;
        *cyc++ = 0x80;

        for ( r = rel + *rel; r >= rel + 1; --r, ++cyc )
            *cyc = INV(*r);
        for ( r = rel + *rel; r >= rel + 1; --r, ++cyc )
            *cyc = INV(*r);
        *cyc++ = 0x80;
    }
    *cyc++ = 0x80;

    /* for all generators set up the CycGen[], CycInd[] and CycLen[]       */
    ind = 0;
    for ( g = 0; g != NrGen; g = (g & 0x20) ? INV(g)+1 : INV(g) ) {

        /* the entry in CycInd[] and CycLen[] for g starts at ind          */
        CycGen[g] = ind;

        /* loop through all relations                                      */
        for ( cyc = CycRels; *cyc != 0x80; ++cyc ) {

            /* find the length of this relation                            */
            for ( len = 0, r = cyc; *r != 0x80; r += 2 )
                ++len;

            /* loop through entries in this relation                       */
            while ( *(cyc+len) != 0x80 ) {
                if ( *cyc == g ) {

                    /* make sure we haven't included this relation yet     */
                    for ( i = CycGen[g]; i < ind; i++ ) {
                        for ( l = 0;  l < CycLen[i];  l++ ) {
                            if ( (CycRels+CycInd[i])[l] != cyc[l] )
                                break;
                        }
                        if ( l == CycLen[i] )
                            break;
                    }
                    if ( i == ind ) {
                        CycInd[ ind ] = cyc - CycRels;
                        CycLen[ ind ] = len;
                        ind = ind + 1;
                    }
                }
                ++cyc;
            }

            /* goto the end of the relation                                */
            cyc = cyc + len;
        }

        /* the entry in CycInd[] and CycLen[] for g is done                */
        CycLen[ ind ] = 0;
        ind = ind + 1;
    }
}


/****************************************************************************
**
*F ReadGrp( <rel>, <name> ) . read the relations from file <name> into <rel>.
** reads the presentation from the file named <name> into the array  <rel>[].
*/
ReadGrp ( rel, name )
        TypGen  * rel;
        char    * name;
{
        char    fname [64];
        FILE *  fptr;
        char    line [1024];
        char    ch;

        strcpy( fname, name );  strcat( fname, ".grp" );
        if ( (fptr=fopen(fname,"r")) == NULL ) {
                fprintf(stderr,"index can't find group file %s.\n",fname);
                exit(1);
        }

        do {
                if ( ! fgets(line,1024,fptr) ) {
                        fprintf(stderr,"index: can't find a PRESENTATION.\n");
                        exit(1);
                }
        } while ( strcmp(line,"PRESENTATION\n") );

        while ( ( ch = getc(fptr) ) != ';' ) {
                if ( 'a' <= ch && ch <= 'z' ) {
                        rel[++rel[0]] = ch-'a';
                        if ( NrGen < ch-'a'+1 ) NrGen = ch-'a'+1;
                }
                else if ( 'A' <= ch && ch <= 'Z' ) {
                        rel[++rel[0]] = ch-'A';
                        if ( NrGen < ch-'A'+1 ) NrGen = ch-'A'+1;
                }
                else if ( ch == '\'' ) {
                        rel[rel[0]] = INV(rel[rel[0]]);
                }
                else if ( ch == ',' ) {
                        rel += rel[0]+1;
                }
                else if ( ch != ' ' && ch != '\t' && ch != '\n' ) {
                        fprintf(stderr,
                            "index: found invalid char %c in relation.\n",ch);
                        exit(1);
                }
        }
}


PrintPres ()
{
    TypCos  i, j;
    TypGen  g;
    int     ident;
    char    t[ 128 ],  s[ 32 ];
    int     new[ 256 ];

    for ( g = 0; g < NrGen; ++g ) {
        for ( i = 1; i <= NrCosets; i++ )  new[i] = 1;
        ident = 1;
        sprintf(t," %c := ",g+'A');
        for ( i = 1; i <= NrCosets; i++ ) {
            if ( Coset[g][i] != i && new[i] != 0 ) {
                ident = 0;
                strcat(t,"(");
                for ( j=i; Coset[g][j]!=i; j=Coset[g][j] ) {
                    sprintf(s,"%lu,",(long)j);
                    strcat(t,s);
                    if ( strlen(t) >= 75 ) {
                        printf("%s\n",t);
                        strcpy(t,"  ");
                    }
                    new[j] = 0;
                }
                sprintf(s,"%lu)",(long)j); strcat(t,s);
                if ( strlen(t) >= 75 ) {
                    printf("%s\n",t);
                    strcpy(t,"  ");
                }
                new[j] = 0;
            }
        }
        if ( ident )  strcat(t,"(1)");
        printf("%s;\n",t);
    }
}
