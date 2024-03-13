/* lib.d file reduce.c */
#include <stdio.h>
#include "defs.h"
#include "list.h"
#include "word.h"
#include "afsa.h"
#include "reduce.h"
int ndiff;
int ngens;
int dollar;
int *** difftable;
int ** fsa;
int ** lhsb, ** lhse, ** rhsb, ** rhse;
int callct = 0;
#define identity 1
#define dt(a,b,c) difftable[a][b][c]
void diff_reduction(DIFF,wp)
	twoafsa * DIFF;
	word * wp;
{
	int i;
	int  *wsaddress, *weaddress, wlm1;
	gen * wpg = wp->g;
	difftable = DIFF->array;
	ndiff = DIFF->states;
	dollar=DIFF->base_symbols;
	ngens = dollar-1;
	wsaddress = wpg + wp->first;
	weaddress = wpg + wp->last;
	if (reducewd(&wsaddress,&weaddress)!=0){
		fprintf(stderr,"\t# Error in reduction procedure.\n");
		exit(2);
	}
	wlm1 = weaddress - wsaddress;
	wp->last=wp->first+wlm1;
}
	
int ***
reduction_rules(rulewords,nrules)
	word ** rulewords;
	int nrules;
{
	int *** rules;
	int i,j;
	rules = vzalloc2(int**,4);
	for (i=0;i<=3;i++)
		rules[i] = vzalloc2(int*,nrules+1); 
	for (j=1;j<=nrules;j++){
		rules[0][j] = (rulewords[0][j]).g + (rulewords[0][j]).first;
		rules[1][j] = (rulewords[0][j]).g + (rulewords[0][j]).last;
		rules[2][j] = (rulewords[1][j]).g + (rulewords[1][j]).first;
		rules[3][j] = (rulewords[1][j]).g + (rulewords[1][j]).last;
	}
	return rules;
}

void
clear_reduction_rules(rules)
	int *** rules;
{
	int i;
	for (i=0;i<=3;i++){
		Free_dp((dp)rules[i]); 
		rules[i]=0;
	}
}

void wa_reduction(WA,rules,wp)
	afsa * WA;
	int *** rules;
	word * wp;
{
	int  *wsaddress, *weaddress, wlm1;
	gen * wpg = wp->g;
	fsa=WA->array;
	lhsb=rules[0]; lhse=rules[1]; rhsb=rules[2]; rhse=rules[3];
	wsaddress=wpg + wp->first;
	weaddress=wpg + wp->last;
	if (reducewd2(&wsaddress,&weaddress)!=0){
		fprintf(stderr,"\t# Error in reduction procedure.\n");
		exit(2);
	}
	wlm1 = weaddress-wsaddress;
	wp->last = wp->first + wlm1;
}
	
/* The next bit of this file is a copy of Derek's file reduceproc */
/* 
INFORMATION on use of procedure reducewd.

This procedure assumes that the following symbols are defined externally.
1. int ngens; The number of generators. It is assumed that generators are
   integers numbered  1,2,3,...,ngens.
2. int ndiff; The number of word differences. It is assumed that these are
   positive integers.
3. int dollar; The number of the string padding generator 'dollar'. This is
   likely to be 0 or ngens+1. Clearly must not lie in range [1,ngens].
4. int identity; This is the number of the (or a) word difference that is
   known to map onto the identity in the group. Probably identity=1.
5. dt(a,b,c) int a,b,c; This function defines the word difference
   transition mapping  A x A U {$} x D -> D U {0}, where  A = set of generators,
   D = set of word differences. Undefined images must map to zero.
   dt must either be defined as an external function or by a macro.

It is assumed that the word to be reduced is stored as a string of integers.
*wordstart and *wordend are pointers to the first and last entries of the
word. (In the case of the empty word, *wordend = *wordstart-1.) The
procedure works by adjusting the input word. *wordstart is not changed, but
if the word gets shorter, then *wordend will decrease. For this reason, the
parameters are double pointers.

reducewd returns 0 on successful completion, and -1 otherwise. The latter can
only happen if a call to malloc fails. No checks are made that the generators,
word differences, and image of  dt  lie in the correct ranges.

It also uses the system call char *malloc() for assigning space, and stderr
for error messages.
*/
   
reducewd(wordstart,wordend) int **wordstart,**wordend;
{ int gmax,gct,*gpref,wordlen,*ws,level,gen1,gen2,diff,diffct,newdiff,olen,
      nlen,i,j;
  char deqi,donesub,*cf,*space1,*space2,*space3,*space4;
  struct vertexd
     { int genno; int diffno; int sublen;
       struct vertexd *backptr;
     } *gptr,*ngptr,*substruc;

/* vertexd is the structure used to store a vertex in the graph of strings
   for possible substitution. The components are as follows.
   backptr - points back to another vertexd, or to zero.
   genno - the generator number,
   diffno - the word difference number of the string defined by following
            backptr back to zero (using genno), relative to the corresponding
            part of the word being reduced.
   sublen - plus or minus the length of this string. sublen is positive if and
            only if the string lexicographically precedes the
            corresponding part of the word being reduced.
   (sublen and subgen1 are put in to save time. Another essential component
    of a vertexd is its level, but we always know this already, because of
    the integers defined by gpref. (See below))
*/
  gmax=65536;
  ws= *wordstart-1; wordlen= *wordend-ws;
  if (wordlen<=0) return(0);

  if ((space1=malloc((unsigned)ndiff))==0)
  { fprintf(stderr,"malloc failed.\n"); return(-1);}
  cf=space1-1;
/* cf is used as a characteristic function, when constructing a subset of the
  set  D  of word differences.
*/

  if ((space2=malloc((unsigned) (wordlen+1)*sizeof(int)))==0)
  { fprintf(stderr,"malloc failed.\n"); return(-1);}
  gpref= (int *) space2; gct= -1; gpref[0]= -1;
/* gpref[n]+1 is the number of vertices that have been defined after reading the
  first n elements of the word. These vertices are gptr[0],...,gptr[gpref[n]].
  We start by allocating space for gmax vertices.
*/

  if ((space3=malloc((unsigned) gmax*sizeof(struct vertexd)))==0)
  { fprintf(stderr,"malloc failed.\n"); return(-1);}
  gptr= (struct vertexd *) space3;

/* Now we start reading the word. */
  level=0;
  while (++level<=wordlen)
  { for (i=1;i<=ndiff;i++) cf[i]=0;
/* Read the element of the word at position level. */
    gen1= ws[level];
/* The next loop is over the identity and the subset of D defined at the
   previous level, level-1.
*/
    diff=identity;
    while (1)
    { deqi= diff==identity;
/* First look for a possible substitution of a shorter string */
      newdiff=dt(gen1,dollar,diff);
      if (newdiff==identity)
/* Make substitution and reduce length of word by 1. */
      { i=level-1;
        if (deqi==0)
        { substruc=gptr+diffct;
          do
          { ws[i]= substruc->genno;
            substruc=substruc->backptr; i--;
          } while (substruc);
        }
        for (j=level;j<wordlen;j++) ws[j]=ws[j+1];
        (*wordend)--; wordlen--;
/* Whenever we make a substitution, we have to go back one level more than
   expected, because of our policy of looking ahead for substitutions
   that reduce the length by 2.
*/
        level= i>0 ? i-1 : i;
        gct=gpref[level]; break;
      }
      else if (newdiff && level<wordlen)
      { j=dt(ws[level+1],dollar,newdiff);
        if (j==identity)
/* Make substitution and reduce length of word by 2. */
        { i=level-1;
          if (deqi==0)
          { substruc=gptr+diffct;
            do
            { ws[i]= substruc->genno;
              substruc=substruc->backptr; i--;
            } while (substruc);
          }
          for (j=level;j<wordlen-1;j++) ws[j]=ws[j+2];
          (*wordend)-=2; wordlen-=2;
          level= i>0 ? i-1 : i;
          gct=gpref[level]; break;
        }
      }

      donesub=0;
/* Now we loop over the generator that is a candidate for substitution
   at this point.
*/
      for (gen2=1;gen2<=ngens;gen2++) if (newdiff=dt(gen1,gen2,diff))
      { if (newdiff==identity)
        { if (deqi)
          { if (gen2<gen1) gen1=ws[level]=gen2;}
          else if (gptr[diffct].sublen > 0)
/* Make a substitution (by a string of equal length). */
          { ws[level]=gen2;
            i=level-1;
            substruc=gptr+diffct;
            do
            { ws[i]= substruc->genno;
              substruc=substruc->backptr; i--;
            } while (substruc);
            level= i>0 ? i-1 : i;
            gct=gpref[level]; donesub=1; break;
          }
        }
        else
        { if (cf[newdiff]) for (i=gpref[level-1]+1;;i++)
/* We have this word difference stored already, but we will check to see if
   the current string precedes the exisiting one.
*/
          { substruc=gptr+i;
            if (substruc->diffno == newdiff)
            { olen=substruc->sublen;
              nlen= deqi ? (gen2<gen1 ? 1 : -1) : 
                       (j=(gptr[diffct].sublen))>0 ? j+1 : j-1;
              if (nlen>olen)
/* The new string is better than the existing one */
              { substruc->genno=gen2;
                substruc->sublen=nlen;
                substruc->backptr= deqi ? 0 : gptr+diffct;
              }
              break;
            }
          }
          else
/* This is a new word difference at this level, so we define a new vertexd in
   graph.
*/
          { gct++;
            if (gct>=gmax)
/* We need more space for vertices. Allocate twice the preceding space and
   copy existing data.
*/
            { if ((space4=malloc((unsigned) 2*gmax*sizeof(struct vertexd)))==0)
              { fprintf(stderr,"malloc failed.\n"); return(-1);}
              printf("Allocating more space.\n");
              ngptr= (struct vertexd *) space4;
              for (i=0;i<gmax;i++)
              { ngptr[i].genno=gptr[i].genno;
                ngptr[i].diffno=gptr[i].diffno;
                ngptr[i].sublen=gptr[i].sublen;
                substruc=gptr[i].backptr;
                if (substruc==0) ngptr[i].backptr=0;
                else for (j=i-1;;j--) if (substruc==gptr+j)
                { ngptr[i].backptr=ngptr+j; break;}
              }
              free(space3); space3=space4; gmax*=2; gptr=ngptr;
            }
/* Define the new vertexd. */
            substruc=gptr+gct;
            nlen= deqi ? (gen2<gen1 ? 1 : -1) : 
                       (j=(gptr[diffct].sublen))>0 ? j+1 : j-1;
            substruc->genno=gen2; substruc->diffno=newdiff;
            substruc->sublen=nlen;
            substruc->backptr= deqi ? 0 : gptr+diffct;
            cf[newdiff]=1;
          }
        }
      } /*End of loop over gen2 */
      if (donesub) break;

/* Go on to next word difference from previous level. */
      if (diff==identity)
      { if (level==1) break;
        diffct=gpref[level-2]+1;
      }
      else diffct++;
      if (diffct>gpref[level-1]) break;
      diff=gptr[diffct].diffno;
    } /* end of loop over word differences at previous level */
    gpref[level]=gct;
  }

  free(space2); 
  free(space3);
  free(space1); 
  return(0);
}


/* The last bit is essentially a copy of Derek's file reduceproc2 */

/* Information on procedure reducewd2.
   Calling it is the same as for reducewd.
   The externally defined symbols are as follows.

   fsa  The word acceptor automaton. fsa[i][j]  (i,j positive) is the
        image of state  j  under generator  i.   It is assumed that
        generators are numbered  1,2,3,...  and that there is no eos-symbol.
        The accept states are numbered  1,2,3,...  All states with negative
        numbers are failure states. fsa[i][j]= -k  means that the left hand
        side of relation number k  has just been read, and so the procedure
        will replace it by its right hand side.

  lhsb, lhse:  lhsb[k]  and  lhse[k]  are pointers to the first and last
               addresses of the word representing the left hand side of
               relation number  k.

  rhsb, rhse:  rhsb[k]  and  rhse[k]  are pointers to the first and last
               addresses of the word representing the right hand side of
               relation number  k.
*/

reducewd2(wordstart,wordend) int **wordstart,**wordend;
/* Reduce "word", by replacing any occurrences of the LHS of the list of
   relations by their RHS. The FSA is used for this. The complete sequence
   of states that it goes through on reading "word" is remembered in the
   array "history".
*/
{ int l,len,st,longer,*midwd,*ptr1,*ptr2,*ptr2e,*ptr,*ws,*we,wordlen,*history;
  char *histspace;
  ws= *wordstart; we= *wordend;
  wordlen= we - ws +1;
  if (wordlen<=0) return(0);
  if ((histspace=malloc((wordlen+1)*sizeof(int)))==0)
  { fprintf(stderr,"malloc failed.\n"); return(-1);}
  history= (int *) histspace;
  midwd=ws; len=0;
  history[0]=1; st=1;

  while (midwd<=we)
  { st=history[++len]=fsa[*midwd][st];
    if (st<0)
/* st= -k means that we have just read the LHS or relation  k. Replace it by
   RHS, and go back to the beginning of that subword.
*/
    { st= -st; ptr1=midwd;
      ptr2=rhsb[st]-1; ptr2e=rhse[st];
      len-=(lhse[st]-lhsb[st]+1);
      midwd=ws+len-1; ptr=midwd;
      while (++ptr2<=ptr2e) *(++ptr)= *ptr2;
      if (ptr!=ptr1)
      { while (++ptr1<=we) *(++ptr)= *ptr1;
        we=ptr; *wordend=we;
      }
      st=history[len];
    }
    midwd++;
  }

  free(histspace);
  return(0);
}
