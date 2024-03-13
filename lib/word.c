/* lib.d file word.c */
#include <stdio.h> 
#include <ctype.h> 
#include "defs.h"
#include "list.h"
#include "word.h"
#include "input.h"
int num_gens=0;
gen * inv_of=0;
word * user_gen_name=0;
word * user_gen_prefix=0;
int * user2prog=0;
int * prog2user=0;
int prefix_len=0;
int max_suffix=0;
int gen_array_size=32;
boolean no_inverses = TRUE;
static void word_expand PARMS((word * wp));
static void insert_gen PARMS((gen h, word * wp));
static void word_stretch_reset PARMS((word * wp, int N, int a));
static void word_double PARMS((word * wp));


elt_fntab
WORD_fntab = {  
  word_sgn_dp,  
  word_cpy_dp,  
  (V2DP)word_create,  
  word_kill_dp,  
  word_print_dp, 
}; 

elt_fntab
GENWT_WORD_fntab = {  
  genwt_word_sgn_dp,  
  word_cpy_dp,  
  (V2DP)word_create,  
  word_kill_dp,  
  word_print_dp, 
}; 

void
pc_word_init (wp) 
word *wp; 
{ 
   wp->g =valloc2(gen,16); 
   wp->first = 0;
   wp->last = -1;
   wp->space = 16;
   wp->type = 0;		/* later set as 's' or 'c' */
   wp->n = 1;			/* only one for trivial or unexpanded words. After initialisation
				   this is set to 0 if the word is not a proper power, and to the power otherwise.
				   */
} 

void
word_init (wp)
word * wp;
{
   wp->g =valloc2(gen,16); 
   wp->first = 0;
   wp->last = -1;
   wp->space = 16;
   wp->type = 's'; 
   wp->n = 0; 
} 

void
word_longinit(wp,n) 
  word *wp; 
  int n;
{ 
  assert(wp);
       wp->g =valloc2(gen,n); 
    wp->first = 0;
    wp->last = -1;
    wp->space = n;
   wp->type = 's'; 
   wp->n = 0; 
} 

void
pc_word_reset (wp) 
  word *wp; 
{ 
   wp->last = wp->first - 1;
   wp->type = 0;
   wp->n = 1;
} 

void
word_reset (wp) 
  word *wp; 
{ 
   wp->last = wp->first - 1;
   if (wp->type != 0){
     wp->type = 's'; 
     wp->n = 0; 
   }
} 

void word_clear(wp) 
    word * wp; 
{ 
    assert(wp); 
    Free_dp((dp)wp->g); wp->g = 0; 
} 

/* 
The word |*wp| has certainly been initialized, and may already contain a
non-trivial word. But the information it contains, if any, is redundant, and
can be thrown away. We just require the space for new information but there
may not be enough of it. If the word currently contains space for less than
n generators we increase the amount of space to this. We set the position
for the first word entry at a, and the last at a-1 (since at this stage the
word is empty).  
*/
static void 
word_stretch_reset(wp,N,a) 
  word *wp;
  int N;
  int a;
{
  if (wp->space < N) {
    Free_dp((dp) wp->g); wp->g = 0; 
    wp->g = valloc2(gen,N); 
    wp->space = N;
  }
    wp->first = a;
    wp->last = a - 1;
}


/* Initialize a |word_traverser| to make it traverse a particular
word, whose address is given in the second argument.
*/
void 
word_traverser_init(wtp,wp)
  word_traverser *wtp; 
  word *wp; 
{ 
  wtp->wp = wp; 
  wtp->posn = wp->first;
} 

/* |word_get_first()| gets hold of the first generator in a word without
changing it. It returns |TRUE| provided the word is non-trivial, |FALSE| 
otherwise.
The  pointer |gp| is set to point to the first generator in the word pointed
to by |wp|.
*/
boolean
word_get_first(wp,gp) 
  word *wp; 
  gen *gp; 
{ 
  assert(wp); 
  assert(wp->first >=0);
  if (wp->first <= wp->last) {
    *gp = (wp->g)[wp->first];
    assert(*gp != INVALID_GEN);
    return TRUE;
  }
  else {
    *gp = INVALID_GEN;
    return FALSE; 
  }
} 

/* |word_next()| is used when traversing a word, generator by generator. It
returns |TRUE| while the |list_traverser| points to a generator, |FALSE| 
once the last generator in the word has been passed.
The  pointer |gp| is set to point to the generator at which the
|word_traverser| is currently positioned..
The procedure does not change the word being traversed.
*/
boolean
word_next(wtp,gp) 
  word_traverser *wtp; 
  gen *gp; 
{ 
  assert(wtp->posn >= wtp->wp->first);
  if (wtp->posn <= wtp->wp->last) {
    *gp = (wtp->wp->g)[wtp->posn]; 
    wtp->posn++;
    return TRUE; 
  } 
  else {
    *gp = INVALID_GEN;
    return FALSE; 
  }
} 


/* Reduce a word.
The word pointed to by |in| is reduced by the cancellation of adjacent
inverse generators to give the word pointed to by |out|. The function 
returns a boolean value, |TRUE| if a reduction has been made so that the
word pointed to by |out| is actually shorter than the word pointed to by
|in|, |FALSE| if no reduction is possible.
Before the function is applied |*out| could be a freshly initialized word or
an actual word already in use. Space is assigned or reassigned as necessary
within this function (via |word_cpy|).
*/
boolean
word_reduce(in,out) 
  word *in, *out; 
{ 
  int i=0; 
  boolean ans = FALSE; 
  gen * outgp=0;

  assert(in); 
  assert(out); 

  outgp = valloc2(gen,in->space); /* we have to do this in case in and
                            out are identical */

  out->space = in->space;
  out->first = in->first;
  out->last = in->last;
  for (i=out->first;i<=out->last;i++)
    outgp[i] = in->g[i];
  
  i = out->first;
  while (i<out->last ) { 
    if (outgp[i] == inv(outgp[i+1])) {
    /* the generators  in positions
|i| and |i+1| are mutually inverse */
      int j;
      ans = TRUE; 
      out->last -=2;
      for ( j = i; j <= out->last; j++) 
        outgp[j] = outgp[j+2]; 
      if (i>out->first)  
        i--; 
    } 
    else 
      i++; 
  }
  Free_dp((dp)out->g); out->g = 0;
  out->g = outgp;
  return ans; 
} 


/* Cyclically reduce a word.
The boolean valued function |word_creduce()| operates just like
|word_reduce()| except that it also cancels inverse pairs of generators at
the two ends of the word. The function operates in two stages. First the
word is reduced using |word_reduce()|, then the inverse generators at the
two ends are cancelled.
*/
boolean word_creduce(in, out) 
 
    word *in, * out; 
{ 
  boolean ans; 

  ans = word_reduce(in,out); 
  while (out->first < out->last
      && out->g[out->first]==inv(out->g[out->last])) {
    out->first++;
    out->last--;      
    ans = TRUE; 
  } 
  return ans;
}


/* Construct the inverse of a word.
If a word is read as a sequence of generators the inverse word consists of
the string of the inverses of these generators, but in the reverse order.
Before the function is applied |*inverse| could be a freshly initialized word or
an actual word already in use. It might be the same word as |*given|.
 Space is assigned or reassigned as necessary
within this function.
*/
void word_inv(given, inverse)
  word *given, *inverse; 
{
  int i;
  int j;
  gen * gp=0;
  gp = valloc2(gen,given->space);
  for (i=given->first,j=given->last;i<=given->last;i++,j--)
    gp[i] = inv(given->g[j]);  
  inverse->space=given->space;
  inverse->first=given->first;
  inverse->last = given->last;
  Free_dp((dp)inverse->g);inverse->g=0;
  inverse->g = gp;
}

boolean word_eqinv(w1p,w2p)
  word * w1p;
  word * w2p;
{
  int m;
  if ((m=((w1p->last) - (w1p->first))) != ((w2p->last) - (w2p->first)))
    return FALSE;
  else if (m>=0){
    int i;
    gen * g1p = w1p->g + w1p->first;
    gen * g2p = w2p->g + w2p->last;
    for (i=0;i<=m;i++)
      if (g1p[i]!=inv(g2p[-i]))
        return FALSE;
  }
  return TRUE;
}

/* Find the length of a word.
The length of a word is the number of non-trivial generators in the string
that makes up the word. The trivial word (the identity) has length 0.
*/
#if 0
int
word_length(wp) 
  word *wp; 
{ 
  int n;
  n = (wp->last) + 1 - (wp->first);
  return n;
} 
#endif

    
/* Adding on a generator to the right hand end of a word (right
multiplication).
*/
void
word_put_last(wp,g) 
  word *wp; 
  gen g; 
{ 
  int i;
  int n;
  n=wp->space;
  if (wp->last == n - 1){ /* there's an entry in the rightmost
piece of space */
    if (wp->first>n/2) {
      int k;
      k=n/4;
      wp->first -=k;
      wp->last -=k;
      for (i=wp->first;i<=wp->last;++i)
        wp->g[i]=wp->g[i+k];
    }
    else
      word_double(wp);
  }
  (wp->last)++;
  wp->g[wp->last] = g;
}   

/*Adding on a generator to the left hand end of a word (left
multiplication).
*/
void
word_put_first(wp,g) 
  word *wp; 
  gen g; 
{ 
  if (wp->first == 0) {
    int i;
    int n;  
    n = wp->space;
    if (wp->last<=n/2) {
      int k;
      k=n/4;
      wp->first +=k;
      wp->last +=k;
      for (i=wp->last;i>=wp->first;--i)
        wp->g[i]=wp->g[i-k];
    }
    else  
      word_double(wp);
  }
  (wp->first)--;
  wp->g[wp->first] = g;
}   

/* Double the amount of space available to a word, while preserving the
information stored in the word.
*/
static void
word_double(wp)
  word * wp;
{
  int n;
  int k;
  int i;
  gen * gp;
  n = wp->space;
  k = n/2;
  gp = valloc2(gen,2*n);
  for (i=k;i<3*k;++i) 
    gp[i] = (wp->g)[i-k];
  Free_dp((dp)wp->g); wp->g=0;
  wp->g = gp;
  wp->space *= 2 ;
  wp->first += k;
  wp->last += k; 
}


/*Get the first generator out of a word, and at the same time delete it from 
that word. 
The function returns |FALSE| if the word is trivial, |TRUE| otherwise.
*/
boolean
word_delget_first(wp,gp) 
  word *wp; 
  gen *gp; 
{ 
  assert(wp); 
  if (wp->first > wp->last) {
    *gp = INVALID_GEN;
    return FALSE;
  }
  else {
    *gp =  wp->g[wp->first++];
    return TRUE; 
  } 
} 
    
/*  Delete the first generator in a word.
The function returns |FALSE| if the word is trivial, |TRUE| otherwise.
*/
boolean
word_del_first(wp) 
  word *wp; 
{ 
  assert(wp); 
  if (wp->first > wp->last) 
    return FALSE;
  else  {
    wp->first++;
    return TRUE; 
  }
} 
    
/* Delete the last generator in a word.
The function returns |FALSE| if the word is trivial, |TRUE| otherwise.
*/
boolean
word_del_last(wp) 
  word *wp; 
{
  assert(wp); 
  if (wp->first > wp->last) 
    return FALSE;
  else  {
    wp->last--;
    return TRUE; 
  }
} 

/* Concatenate two words.
|wp1|, |wp2| and |wp3| are pointers to initialized words. |*wp1| and |*wp2| are
not changed by this procedure (provided that they are both distinct from
|*wp3|). At the end |*wp3| contains the concatenation
of |*wp1| and |*wp2|.
*/
void
word_concat(wp1,wp2,wp3) 
  word *wp1, *wp2, *wp3; 
{ 
  int n1,n2,n; 
  int N;
  int i;
  gen *factor,*concat; 

  n1 = word_length(wp1);
  n2 = word_length(wp2);
  n = n1 + n2;
  if (n<=0){
    wp3->last = wp3->first -1;
    return;
  } /* avoid malloc()ing 0 bytes */
  concat = valloc2(gen,n); 
  factor = wp1->g + wp1->first; 
  for (i=0;i<n1;i++)
    concat[i] = factor[i];
  factor = wp2->g + wp2->first; 
  for (i=n1;i<n;i++)
    concat[i] = factor[i-n1]; 
  N = max(wp1->space,wp2->space);
  while (N<2*n)
    N *=2;
  word_stretch_reset(wp3,N,n/2);
  for (i=0;i<n;++i)
    wp3->g[wp3->first + i] = concat[i];
  wp3->last = wp3->first + n - 1;
  Free_dp((dp)concat); concat = 0;
}   

void
word_insert(wp,wwp)
        word * wp, * wwp;
{
        int N = wp->space;
        int first = wp->first; /* the free space at the beginning of wp */
        int last = wp->last;
        int n = (wwp->last) - (wwp->first) + 1; /* the length of wwp */
        int i;
        gen * genp1, * genp2, * genp3;
        if (first  <  n){
                gen * genp;
/* there isn't enough free space at the beginning of
 wp for wwp, but it may be possible to create space by shifting wp towards the
 end of its space. If this isn't the case then we increase the amount
 of space available. */
                while (first + N - last <n)
                        N *= 2;
                genp=vzalloc2(gen,N);
                for (i=wp->first;i<=last;i++)
                        genp[i+n]=(wp->g)[i];
                Free_dp((dp)(wp->g));
                wp->g = genp;
                wp->space = N;
                wp->first += n;
                wp->last += n;
        }
        genp1 = (wp->g) + first;
        genp3 = genp1 - n;
        genp2 = (wwp->g) + (wwp->last);
        while (--genp1>=genp3) *genp1 = *(genp2--);
        wp->first -= n;

}

void
word_append(wp1,wp2)
  word * wp1, * wp2;
{
  int N = wp1->space;
  int last = wp1->last;
  int n2 = (wp2->last) - (wp2->first) + 1;
  int i;
  gen * genp1, * genp2, * genp3;
  while (N - 1 - last <= n2)
    N *= 2;
  if (N > wp1->space){
    gen * genp=vzalloc2(gen,N);
    for (i=wp1->first;i<=last;i++)
      genp[i]=(wp1->g)[i];
    Free_dp((dp)(wp1->g));
    wp1->g = genp;
    wp1->space = N;
  }
  genp1 = (wp1->g) + last;
  genp3 = genp1 + n2;
  genp2 = (wp2->g) + (wp2->first);
  while (++genp1<=genp3) *genp1 = *(genp2++);
  wp1->last += n2;
}

void
word_invappend(wp1,wp2)
  word * wp1, * wp2;
{
  int N= wp1->space;
  int last = wp1->last;
  int n2 = (wp2->last) - (wp2->first) + 1;
  int i;
  gen * genp1, * genp2, * genp3;
  while (N - 1 - last <= n2)
    N *= 2;
  if (N > wp1->space){
    gen * genp=vzalloc2(gen,N);
    for (i=wp1->first;i<=wp1->last;i++)
      genp[i]=(wp1->g)[i];
    Free_dp((dp)(wp1->g));
    wp1->g = genp;
    wp1->space = N;
  }
  genp1 = (wp1->g) + last;
  genp3 = genp1 + n2;
  genp2 = (wp2->g) + (wp2->last);
  while (++genp1<=genp3) *genp1 = inv(*(genp2--));
  wp1->last += n2;
}

/* The signature of an ordering of two words. 
Words are ordered by length and lexicographically (that is, w1 is
earlier than w2 in the ordering if it is shorter than w2 or if the two words
have the same length but in the first position where they differ the
generator in w1 is lower valued than the generator in w2). This function returns
1 if the two words are given with the earlier ordered one first, -1 if the
earlier ordered comes second, 0 if the two words are equal.
The arguments of the function must be datapointers, since this is a function
in the function table for lists of words; another version of the function
exists which uses pointers to words as its arguments. 
 We compare the two words lexicographically first, calculating the
difference between them in the first position where they differ. Then we
compare their lengths. Lexicographic comparison is only relevant for two words
of the same length.
*/
int
word_sgn_dp(dtp1,dtp2) 
  dp dtp1, dtp2; 
{ 
  int diff=0; 
  int ans = 0; 
  int n1, n2;
  int i=0;
  word * wp1, *wp2;
  gen * gp1, * gp2; 

  wp1 = (word *)dtp1;
  wp2 = (word *)dtp2;

  n1 = word_length(wp1);
  n2 = word_length(wp2);
  gp1 = wp1->g + wp1->first; 
  gp2 = wp2->g + wp2->first; 
  if (n1 < n2)
    ans = 1;
  else if (n1 > n2)
    ans = -1;
  else  { /* the two words have the same length */
    while(diff == 0 && i < n1) { 
      diff = (int)gp1[i] - (int)gp2[i]; 
      i++;
    } 
    if (diff > 0) 
      ans = -1;  
    else if (diff < 0) 
      ans = 1; 
  }
  return ans; 
} 

int 
genwt_word_sgn_dp(dp1,dp2)
  dp dp1, dp2; 
{ 
  int ans = 0; 
  word * wp1, *wp2;
  int m1=0;
  int m2=0;
  word_traverser wt;
  gen g;

  wp1 = (word *)dp1;
  wp2 = (word *)dp2;

  word_traverser_init(&wt,wp1);
  while (word_next(&wt,&g))
    if (g>m1)
      m1=g;  
  word_traverser_clear(&wt);
  word_traverser_init(&wt,wp2);
  while (word_next(&wt,&g))
    if (g>m2)
      m2=g;  
  word_traverser_clear(&wt);
  if (m1 < m2)
    ans = 1;
  else if (m1 > m2)
    ans = -1;
  else  
    ans=word_sgn(wp1,wp2);
  return ans; 
} 

/* Copy the word pointed at by |oldp| (a data pointer) into the address
specified by |newp| (also a data pointer).
Before the application of the function |newp| can point either to a
freshly initialized (empty) word or to an actual word of any length.
We need to interpret the data pointers as pointers to  words,  |oldwp| and
|newwp|. Once we have reallocated space (according to the length of the word
|*oldwp|) to the pointer |newwp| we copy the word across one generator at a 
time. 
*/
void word_cpy_dp(oldp,newp) 
  dp oldp, newp; 
{ 
    int m;
    int n;
    int i;
    word *newwp, *oldwp; 
    gen *newgp, *oldgp; 
    
    assert(newp); 
    assert(oldp);
    if (newp == oldp) return; /* nothing to do */ 
    newwp = (word*)newp; 
    oldwp = (word*)oldp;
    m = word_length(oldwp); 
    if ((n=newwp->space) < 2*m) {
      while (n<2*m)
        n *= 2;
      Free_dp((dp)newwp->g); newwp->g = 0; 
      newwp->g = valloc2(gen,n); 
      newwp->space = n;
    }
    newwp->first = m/2;
    newwp->last = newwp->first+m-1;
    newgp = newwp->g + newwp->first; 
    oldgp = oldwp->g + oldwp->first; 
    for (i=0;i<m;++i)
      newgp[i]=oldgp[i];
} 

void word_mv(oldwp,newwp)
word * oldwp, * newwp;
{
  if (newwp->g)
   Free_dp((dp)newwp->g);
  newwp->g = oldwp->g;
  newwp->first = oldwp->first;
  newwp->last = oldwp->last;
  newwp->space = oldwp->space;
  oldwp->g = 0;
}

/* Create initialized space for a word.
The function returns a pointer of type |word *| assigned initialized space equal
in size to the size of a word.
*/
word *  word_create()

{ 
    word * wp; 
    wp = vzalloc1(word); 
    word_init(wp); 
    return (wp); 
} 

/* Retrieve the space occupied by a word.
Both the space the word itself occupies and the space assigned to it as a
pointer to generators must be reclaimed. 
*/
void word_kill_dp(dtp)
  dp dtp; 
{ 
  word_clear((word *) dtp); 
  Free_dp(dtp); 
} 

/* Print a word.
This prints out a word  currently stored as a string of generators in the
program's notation for generators as a string of generators 
in the notation used by the user.
This is a function table function, so accepts a datapointer argument. It
also exists coerced to accept a |word *| argument.
*/
void word_print_dp(wfile,dtp) 
  FILE * wfile;
  dp dtp; 
{ 
    word * wp; 
    word user_word;
    word_init(&user_word);
    wp = (word *)dtp; /* coerce to type |word *| */ 
    word2user_name(wp,&user_word);
    user_word_print(wfile,&user_word);
    word_clear(&user_word);
} 

/* Print out a word. This is the same as the above except that it takes a
|word *| rather than a |dp| argument.
*/
void word_print(wfile,wp) 
  FILE * wfile;
  word * wp; 
{ 
    word user_word;
    word_init(&user_word);
    word2user_name(wp,&user_word);
    user_word_print(wfile,&user_word);
    word_clear(&user_word);
} 

/* Print a word.
This prints out a word as a string of generators interpreted as alphabet
characters.
*/

void user_word_print(wfile,wp)
  FILE * wfile;
  word * wp;
{
    gen * gp;
    int i;
    int ct=0;
    gp = wp->g;
    assert(gp);
    for (i=wp->first;i<=wp->last;i++) {
      assert(isascii(gp[i]) && isprint(gp[i]));
      if (gp[i]=='$')
        fprintf(wfile,"$");
      else
        fprintf(wfile,"%c",gp[i]);
      ct++;
      if (gp[i]=='*'){
        if (ct>=72){
           fprintf(wfile,"\n  ");
           ct=1;
        }
      }
    }
    if (word_length(wp)==0)
      fprintf(wfile,"1");
}





boolean read_next_gen(wp,file)
  word * wp;
  FILE * file;
{
  boolean ans=TRUE;
  int c;
  word_reset(wp);
  while ((c=read_char(file))!=EOF && !(isalpha(c) || c=='_')){
    if (c=='$') 
      return TRUE;
    if (c=='}') {
  /* a '}' is used to terminate the input of gens */
      ungetc('}',file);
      ans = FALSE;
      break;
    }
    else if (c!=' ' && c!=','){
      fprintf(stderr,"Generators must start with a letter or underscore.\n");
      bad_data();
    }
  }
  if (c==EOF){
    fprintf(stderr,"Unexpected end of file.\n");
    bad_data();
  }
  if (ans==TRUE){
    do {
      if (c==EOF)
        break;
      if (c=='^'){
        int d,e;
        d=read_char(file); e=read_char(file);
        if (d!='-' || e!='1'){
          fprintf(stderr,
"Invalid generator name.\n");
          bad_data();
        }
        word_put_last(wp,c);
        word_put_last(wp,d);
        word_put_last(wp,e);
        c=read_char(file);
        break;
      }
      
      if (!isalpha(c) && !isdigit(c) && c!='_' && c!='.'){
        fprintf(stderr,
"Only letters, digits, .'s and underscores are allowed in generator names.\n");
        bad_data();
      }
      word_put_last(wp,c);
      c = read_char(file);
    } while (isalpha(c) || isdigit(c)
      ||c=='^'||c=='-'||c=='_'||c=='.');
    if (c!=EOF && c!= ' ')
      ungetc(c,file);
  }
  return ans;
}

char read_next_word (wp,file)
  word * wp;
  FILE * file;
{
/*
A word is recognisable as a power of a single commutator iff
a) the first symbol after any number of ('s is a [ and
b) it contains no *'s, no ['s or ^'s inside the []'s, and only ('s, )'s ^'s -'s
and digits outside.
wp->type is initialised as 0, then we reset it
when we encounter the first symbol that is not a (. If set to
c it will be reset to s once an inappropriate symbol is encountered.
is encountered .
A word is recognised as a power if it contains no *'s
unless enclosed by ()'s. wp->n is reset to 0 as soon
as such a '*' is encountered.
*/
  char ans=1;
  int c;
  word_reset (wp);
  while ((c=read_char (file))!=EOF && !(isalpha(c)) &&!(isdigit(c))&&
    c!= '_' && c!='(' && c!='['){
    if (c=='$') 
      return 1;
    if (c=='}') {
  /* a '}' is used to terminate the input of words */
      ungetc('}',file);
      ans = 0;
      break;
    }
  }
  if (c==EOF){
    fprintf(stderr,"Unexpected end of file.\n");
    bad_data ();
  }
  if (ans==1 && c=='1') return ans;
/* 1 is used to represent the identity. It will only be the first symbol in a 
word when it's the whole word. */
  if (ans==1){
    gen separator = '*'; 
    int bracket_level=0;
    int comm_level=0;
    do {
      if (c==EOF)
        break;
      if (c==' '||(c==','&&comm_level==0)){
/* if this happens, the previous character must have been an
alphanumeric character or a ), and we are not within commutators (since 
spaces after ^,( and * and inside commutators are dealt with
elsewhere). If the next non-space character is not *,^ or ) we must have
reached the end of the word. */
        do {
          c=read_char (file);
          if (c==EOF) break;
        } while (c==' ');
        if (c==EOF||(c!='^'&&c!='*'&&c!=']'&&c!=')')) 
          break; /* we have got to the end of the word */
      }
      if (c=='['){
        comm_level++;
        if (wp->type!=0) wp->type = 's';
        if (wp->type==0) wp->type = 'c';
      }
      else if (c==']')
        comm_level--;  
      else if (c=='(')
        bracket_level++;
      else if (c==')')
        bracket_level--;  
      if ((c=='*')|| (wp->type==0 && c!='(')|| (comm_level==0 && c!=']' 
                     && c!='('&& c!=')' && c!='^' && c!='-' && (!isdigit(c))))
        wp->type = 's';
      if (bracket_level==0 && comm_level==0 && c=='*') wp->n=0;
      if (bracket_level>=0 && comm_level>=0)
        word_put_last (wp,c);
      else break;
/* word reading is terminated by an unmatched right bracket of either type */
      if (c=='*'||c=='^'||c=='('|| comm_level!=0){
/* '*','^' and '(' are non-terminal characters, so we can automatically 
skip over any
white spaces that appear after them, without further investigation */
         do {
           c=read_char (file);
           if (c==EOF) break;
         } while (c==' ');
       }
       else { 
         c = read_char (file);
         if (c==EOF) break;
       }
    } while (isalpha(c) ||isdigit(c) || c=='_' ||c=='.'
      ||c=='^'||c=='-'||c=='('||c==')'||c=='['||
             c==']'||c==','||c=='*'||c==' ');
    if (c!=EOF) ungetc(c,file);
  }
  if (c==EOF){
    fprintf(stderr,"Unexpected end of file.\n");
    bad_data ();
  }
  else if (word_length(wp)==0) 
    wp->type = 's'; /* then we'd have the trivial word */
  else word_expand (wp);
  return ans;
}

static void
word_expand (wp) word *wp;
{
  word expansion;
  word buf;
  word temp;
  word temp2;
  gen g;
  char type = wp->type;
  int n = wp->n;
  int separator= '*';
  int bracket_level=0;
  int commutator_level=0;
  int exponent;
  int sign=1;
  if (type=='c'){
    if (n){
      int i = wp->first; /* i is the position we're currently at in gp */
      gen * gp = wp->g + wp->first;
      gen * ggp = wp->g + wp->last;
      while (gp<=ggp){
        if ((*gp)=='^'){
          int j=i;
          exponent = 0;
          gp++; i++; 
          if ((*gp)=='-') sign = -1; else exponent = (*gp) - '0';
          gp++; i++;
          while (gp<=ggp && isdigit(*gp)){
            exponent = 10*exponent + (*gp) - '0'; gp++; i++;
          }
          exponent *= sign;
          if (commutator_level==0){ 
/* in the case the whole commutator is raised to a power */
            wp->last = j-1; /* this is to delete the exponent from the end */
            if (exponent==0) word_reset(wp);
            else n *= exponent;
            break;
          }
          else if (exponent != -1){
 /* if we see any exponent other than a -1 inside a commutator we expand it */
             wp->type = type = 's';
             break;
           }
        }
        else {
          if ((*gp)=='[') commutator_level=1;
          else if ((*gp)==']') commutator_level=0;
          gp++; i++; 
        }
      }
      wp->n = n;
    }
    if (wp->type == 'c') return;
  }
  word_init (&expansion);
  word_init (&buf);
  word_init (&temp);
  word_init (&temp2);
  while (word_delget_first (wp,&g)){
    if (isalpha(g) || isdigit(g) || g=='_' ||g=='.')
      word_put_last (&buf,g);
    else if (g==separator){
      if (word_length(&buf)>0){
        word_append (&expansion,&buf);
        word_put_last (&expansion,g);
        word_reset (&buf);
      }
    }
    else if (g=='^'){
/*  read the exponent that follows and replace whatever's currently in
the buf by the appropriate number of copies of it or its inverse */ 
      int i;
      sign=1;
      (void)word_delget_first (wp,&g);
      if (g=='-'){
        sign = -1;
        word_delget_first (wp,&g);
      }
      if (isdigit(g)){
        exponent=0;
      /*Now g must be a digit */
        exponent= g - '0' + 10*exponent;
        while (word_get_first (wp,&g) && isdigit(g)){
            exponent= g - '0' + 10*exponent;
            word_delget_first (wp,&g);
        }
        if ((n) && word_length(wp)==0){
/* wp->n non-zero means we believe we're looking for a power at the end of 
the word. We know we are at the end of the word if we've deleted every character */
          if (exponent==0) word_reset(wp);
          else n = sign*exponent;
        }
        else {
          if (sign == -1) {
            gen h;
            word2prog_word (&buf,&temp);
            if (word_length(&temp)==1 && word_get_last (&temp,&h) &&
                 (inv_of==0 || inv_of[h]==0)){
  /* in this case we're reading the name of a new generator */
              word_put_last (&buf,'^');
              word_put_last (&buf,'-');
              word_put_last (&buf,'1');
            }
            else {
              word_inv (&temp,&temp);
              word2user_name (&temp,&buf);
            }
          }
          word_cpy (&buf,&temp);
          word_reset (&buf);
          for (i=1;i<=exponent;i++){
            word_append (&buf,&temp);
            if (i<exponent && word_length(&buf)>0)
              word_put_last (&buf,separator);
          }
          word_reset (&temp);
        }
      }
      else {
        word exp;
        bracket_level=0;
        commutator_level=0;
        word_init (&exp);
        do{
          if (commutator_level==0 && bracket_level==0 &&
             (g=='*'||g=='^')){
            word_put_first (wp,g);
            break;
          }
          if (g=='(')
            bracket_level++;
          else if (g==')')
            bracket_level--;
          if (g=='[')
            commutator_level++;
          else if (g==']')
            commutator_level--;
          word_put_last (&exp,g);
        } while (word_delget_first (wp,&g));
        word_expand (&exp);
        if (word_length(&exp)>0){
          word2prog_word (&exp,&temp);
          word_inv (&temp,&temp);
          word2user_name (&temp,&temp2);
          word_put_last (&temp2,separator);
          word_concat (&temp2,&buf,&buf);
          word_put_last (&buf,separator);
          word_append (&buf,&exp);
          word_reset (&temp);
          word_reset (&temp2);
        }
        word_clear (&exp);
      }
    }
    else if (g=='(') {
      if (word_length(&buf)>0) bad_data();
/* the word buf should be empty at this stage. If it isn't, it could
be that the user has forgotten to be a * between two subwords */
      bracket_level = 1;
      while (word_delget_first (wp,&g)){
        if (g=='(')
          bracket_level++;
        else if (g==')')
          bracket_level--;
        if (bracket_level==0)
          break;
        word_put_last (&buf,g);
      }
      word_expand (&buf);
    }
    else if (g=='[') {
      word w1, w2;
      if (word_length(&buf)>0) bad_data();
/* the word buf should be empty at this stage. If it isn't, it could
be that the user has forgotten to be a * between two subwords */
      commutator_level = 1;
      word_init (&w1);
      while (word_delget_first (wp,&g)){
        if (g==',' && commutator_level==1)
          break;
        else if (g=='[')
          commutator_level++;
        else if (g==']')
          commutator_level--;
        word_put_last (&w1,g);
      }
      word_expand (&w1);
      word_init (&w2);
      while (word_delget_first (wp,&g)){
        if (g=='[')
          commutator_level++;
        else if (g==',' && commutator_level==1){
/* convert to nested binary commutators, to make expansion easier */
          word_put_first(&w1,'[');
          word_put_last(&w1,',');
          word_append(&w1,&w2);
          word_put_last(&w1,']');
          word_expand(&w1);
          word_reset(&w2);
          continue; /* we don't want to put the ',' onto the end of w2 */
        }
        else if (g==']'){
          commutator_level--;
          if (commutator_level==0)
            break;
        }
        word_put_last (&w2,g);
      }
      word_expand (&w2);
      if (word_length(&w1)>0 && word_length(&w2)>0){
        word2prog_word (&w1,&temp);
        word_inv (&temp,&temp);
        word2user_name (&temp,&buf);
        word_put_last (&buf,separator);
        word2prog_word (&w2,&temp);
        word_inv (&temp,&temp);
        word2user_name (&temp,&temp2);
        word_append (&buf,&temp2);
        word_put_last (&buf,separator);
        word_append (&buf,&w1);
        word_put_last (&buf,separator);
        word_append (&buf,&w2);
        word_reset (&temp);
        word_reset (&temp2);  
      }
      word_clear (&w1); word_clear (&w2); 
    }
  }
  if (bracket_level!=0 || commutator_level!=0){
    fprintf(stderr,"Unmatched bracket in relation.\n");
    bad_data ();
  }
  word_append (&expansion,&buf);
/* there's no * at the end */
  if (word_get_last (&expansion,&g) && g==separator) word_del_last (&expansion);
  word_cpy (&expansion,wp);
  word_clear (&expansion);
  word_clear (&temp);
  word_clear (&temp2);
  word_clear (&buf);
  if (n==1) wp->n = 0;
  else wp->n = n;
  wp->type = type;
}

/*
\Pre |wp| points to an initialized word.
\Post |gp| points to the last letter of |*wp|, |*wp| is unchanged.
word.
\Returns |TRUE| if the word pointed to by |wp| is non-empty, and |FALSE|
if it is empty.
*/
boolean
word_get_last(wp,gp)
  word *wp;
  gen *gp;
{
  boolean ans = TRUE;
  if (wp->first > wp->last) {
    *gp = INVALID_GEN;
    ans = FALSE;
  }
  else
    *gp = wp->g[wp->last];
  return ans;
}
  
/*
\Pre |wp| points to an initialized word.
\Post |gp| points to the last letter of |*wp|, which is deleted from that
word.
\Returns |TRUE| if the word pointed to by |wp| is non-empty, and |FALSE|
if it is empty.
*/
boolean
word_delget_last(wp,gp)
  word *wp;
  gen *gp;
{
  boolean ans = TRUE;
  if (wp->first > wp->last) {
    *gp = INVALID_GEN;
    ans = FALSE;
  }
  else {
    *gp = wp->g[wp->last];
    wp->last--;
  }
  return ans;
}
  
void
gen_print(wfile,g)
  FILE * wfile;
  gen g;
{
  word w;
  word_init(&w);
  if (g==IDENTITY)
    fprintf(wfile,"1");
  else if (g==num_gens+1)
    fprintf(wfile,"$");
  else {
    gen2user_name(g,&w);
    user_word_print(wfile,&w);
  }
  word_clear(&w);
}
    
int
gen_sgn_dp(gp1,gp2)
  dp gp1, gp2; /*really these are pointers to |gen|'s*/
{
  boolean ans = 0;
  gen *g1, *g2;
  g1 = (gen *)gp1;
  g2 = (gen *)gp2;
  if (*g1 < * g2)
    ans = 1;
  else if (*g1 > *g2)
    ans = -1;
  return ans;
}

void
gen_cpy_dp(oldp,newp)
  dp oldp, newp; /*these are really |gen *|'s*/
{
  gen * newg = (gen *)newp;
  gen * oldg = (gen *)oldp;
  *newg = *oldg;
}
gen *
gen_create()
{
  return valloc1(gen);
}
void
gen_print_dp(wfile,gp)
  FILE * wfile;
  dp gp; /*really a |gen *|*/
{
  word w;
  word_init(&w);
  gen2user_name(*((gen *)gp),&w);
  user_word_print(wfile,&w);
  word_clear(&w);
}

void
genstring(string,g)
  char * string;
  gen g;
{ 
  if (g==0){
    string[0]='0';
    string[1]='\0';
  }
  else {
    word w;
    gen h;
    int length;
    int i=0;
    word_init(&w);
    gen2user_name(g,&w);
    length = word_length(&w);  
    while(word_delget_first(&w,&h))
    string[i++]=(char)h;
    string[length]= '\0';
    word_clear(&w);
  }
}

void
word2prog_gen(user_namep,prog_genp) 
  word* user_namep;  
  gen * prog_genp;
{ 

  gen g; 
  gen * g1p, * g2p;
  int j;

  *prog_genp=INVALID_GEN;
  if (user2prog){
      g1p = user_namep->g + user_namep->first + prefix_len;
      g2p = user_namep->g + user_namep->last;
      j = 0;
      while (g1p <= g2p && isdigit(*g1p)){
        j = 10*j + (*g1p - '0');
        g1p ++;
      }
      if (j<=max_suffix){
        if (g1p < g2p) /* we had an inverse */
          *prog_genp = inv_of[user2prog[j]];
        else *prog_genp = user2prog[j];
      }
  }
  else {
    for (g=1;g<=num_gens;++g) 
       if (word_sgn(user_namep,user_gen_name+g)==0) { 
             *prog_genp = g; 
             break; 
       } 
  }
}  

void
word2prog_word(user_wordp,prog_wordp)
  word * user_wordp;
  word * prog_wordp;
{
  word  user_gen;
  word_traverser wt;
  gen g;
  gen h;
  int i;
  int n=word_length(user_wordp);
  char* epsilon ="epsilon";
  word identity;
  word_init(&identity);
  for (i=0;i<=6;i++)
    word_put_last(&identity,(gen)epsilon[i]);
  i=0;
  word_reset(prog_wordp);
  prog_wordp->type = user_wordp->type;
  prog_wordp->n = user_wordp->n;
  if (word_sgn(user_wordp,&identity)!=0){
    word_init(&user_gen);
    word_traverser_init(&wt,user_wordp);
    while (word_next(&wt,&g)){
      i++;
/* the  symbols '[' ']' '(,')' and ',' might appear if we weren't expanding 
commutators - at the moment 
we're not leaving commutators unexpanding, but we may want to do so later */
      if (g!='*' && g!=',' && g!='[' && g!=']' && g!='(' && g!=')')
        word_put_last(&user_gen,g);
      if ((g=='*'||g==','||g==']'||g==')'||i==n) && word_length(&user_gen)!=0){
        word2prog_gen(&user_gen,&h);
        word_reset(&user_gen);
        word_put_last(prog_wordp,h);
      }
    }
    word_traverser_clear(&wt);
    word_clear(&user_gen);
  }
  word_clear(&identity);
}

void
gen2prog_gen(user_gen,prog_genp)
  gen user_gen;
  gen * prog_genp;
{
  word w;
  word_init(&w);
  word_put_last(&w,user_gen);
  word2prog_gen(&w,prog_genp);
  word_clear(&w);
}


void
gen2user_name(program_gen,user_wordp) 
  gen program_gen;  
  word * user_wordp; 
{ 

  word_reset(user_wordp);
  if (program_gen==0||(num_gens!=0 && program_gen>num_gens)){
    gen g='$';
    word_put_last(user_wordp,g);
  }
  else 
     word_cpy(user_gen_name+program_gen,user_wordp); 

}  

void
word2user_name(prog_wordp,user_wordp)
  word * prog_wordp, * user_wordp;
{
  word w;
  gen separator='*';
  gen g=INVALID_GEN;
  word_traverser wt;
  word_traverser_init(&wt,prog_wordp);
  word_reset(user_wordp);
  word_init(&w);
  while (word_next(&wt,&g)){
    gen2user_name(g,&w);
    word_append(user_wordp,&w);
    word_put_last(user_wordp,separator);
    word_reset(&w);
  }
  word_traverser_clear(&wt);
  word_del_last(user_wordp);
  word_clear(&w);
}

/* if user_gen_prefix is non-zero at the end of this function, then every input 
generator is of the form prefix%d or prefix%d^-1 for a common string prefix which 
is the content of the word pointed at by user_gen_prefix */
void
read_gen_name_array(file)
  FILE * file;  
{
  word w;
  int i,j,k;
  gen * g1p, * g2p, * g3p, * g4p;
  num_gens=0;
  user_gen_name = vzalloc2(word,gen_array_size);
  user_gen_prefix = vzalloc1(word);
  for (i=0;i<gen_array_size;++i)
    word_init(user_gen_name+i);
  find_char('{',file);
  word_init(&w);
  word_init(user_gen_prefix);
  while (read_next_gen(&w,file)){
    define_next_gen(&w);
    if ((no_inverses) && word_length(&w)>=3 && (w.g)[(w.last)-2]=='^')
        no_inverses = FALSE;
    if (num_gens==1 || user_gen_prefix){
      j = w.last;
      if (word_length(&w)>=3 && w.g[j]=='1' && w.g[j-1]=='-' && w.g[j-2]=='^') j -= 3;
      k = j;
      while (isdigit((w.g)[j])) j--;
      if (j<k){
        if (num_gens==1){
          word_cpy(&w,user_gen_prefix);
          user_gen_prefix->last = user_gen_prefix->first + j - w.first;
/* the word pointed at by user_gen_prefix is set to prefix, where w = prefix%d or
prefix %d^-1 */
        }
        else {
          g1p = w.g + w.first; g2p = w.g+j;
          g3p = user_gen_prefix->g + user_gen_prefix->first;
          g4p = user_gen_prefix->g + user_gen_prefix->last;
          while (g3p <= g4p){
            if (g1p > g2p || *g1p != *g3p) break;
            else { g1p++; g3p++; }
          }  
          if (g3p<=g4p) { 
/* the word w is not of the form prefix%d or prefix%d^-1 */
            word_clear(user_gen_prefix);
            Free_dp((dp) user_gen_prefix);
            user_gen_prefix = 0;
          }
        }
      }
      else if (user_gen_prefix) {
/* the word w just read in has no terminating digits, or none before a terminating "^-1"
*/
        word_clear(user_gen_prefix); 
        Free_dp((dp) user_gen_prefix);
        user_gen_prefix = 0;
      }
    }
    word_reset(&w);
  }
  word_clear(&w);
  find_char('}',file);
  if (user_gen_prefix) prefix_len = word_length(user_gen_prefix);
}

void
define_next_gen(wp)
  word * wp;
{
  word * ucpy;
  int size = gen_array_size;
  num_gens++;
  if (num_gens>=gen_array_size){ 
  /* need more space */
    int i;
    gen_array_size *= 2;
    ucpy = vzalloc2(word,gen_array_size);
    for (i=0;i<gen_array_size;++i) 
      word_init(ucpy+i);
    for (i=0;i<num_gens;++i){ 
      word_cpy(user_gen_name+i,ucpy+i);
      word_clear(user_gen_name+i);
    }
    Free_dp((dp)user_gen_name); user_gen_name=0;
    user_gen_name = ucpy;
    ucpy=0;
  }
  word_cpy(wp,user_gen_name+num_gens);
}

void
read_inverse_array(file)
  FILE * file;
{
  int c;
  gen g,h,k;
  int i,j;
  gen * g1p, * g2p;
  word w;
  if (inv_of){
    find_char('{',file);
    find_char('}',file);
    return;
  }
  inv_of=vzalloc2(gen,gen_array_size);
  word_init(&w);
  find_char('{',file);
  while ((c=read_char(file))==' ');
  if (c=='c'){ /* first letter of case_change */
    gen l=INVALID_GEN;
    for (g = 1; g <= num_gens; ++g) { 
      if (inv_of[g]!=0)
        continue; /* entries are defined in pairs, so this entry may
have already been set */
      word_cpy(user_gen_name+g,&w);
      if (word_length(&w)!=1){
        fprintf(stderr,
          "Case change inverses are only defined for generators\n");
        fprintf(stderr,
              "which are single alphabet characters.\n"); 
        bad_data();
      }
      (void)word_delget_first(&w,&h);
      if (islower(h))
        k=toupper(h);
      else
        k=tolower(h);
      gen2prog_gen(k,&l);
      if (l==INVALID_GEN){
/* There's no generator yet defined with user name k, so we have to slot one
in, just after g.  */
        word_put_first(&w,k);
        l=g+1;
        insert_gen(l,&w);
      }
      inv_of[g] = l; inv_of[l] = g;
      word_reset(&w);
    }
  }
  else if (c=='f'){ /* first letter of formal */
    default_inverse_array();
  }
  else if (c=='i'){ /*first letter of inv */
    char * label;
    label=vzalloc2(char,4);
    ungetc(c,file);
    while (read_next_string(label,3,file)){
      if (strcmp(label,"inv")!=0){
        fprintf(stderr,
"Syntax error in inverse table.\n");
        bad_data();
      }
      find_char('\(',file);
      read_next_gen(&w,file);
      word2prog_gen(&w,&g);
      if (g==INVALID_GEN){
        fprintf(stderr,
"Undeclared generator in inverse table.\n");
        bad_data();
      }
      find_char('=',file);
      read_next_gen(&w,file);
      word2prog_gen(&w,&k);
      if (k==INVALID_GEN){
/* There's no generator yet defined with user name w, so we have to slot one
in, just after g.  */
        insert_gen(g+1,&w);
        k = g+1;
      }
      inv_of[g]=k; inv_of[k]=g;
    }
    Free_dp((dp)label);
    no_inverses = FALSE;
    default_inverse_array();
  }
  else  {
    if (c!=' ') ungetc(c,file);
    g=1;
    while (read_next_int(inv_of+g,file))
      g++;
    if (g==1) /* i.e. no entries */
      default_inverse_array();
  }
  word_clear(&w);
  find_char('}',file);
  if (user_gen_prefix && prog2user==0 ){
    prog2user = vzalloc2(int,gen_array_size);
    word_init(&w);
    for (i=1;i<=num_gens;i++){ 
      g1p = user_gen_name[i].g + user_gen_name[i].first + prefix_len;
      g2p = user_gen_name[i].g + user_gen_name[i].last;
      j = 0;
      while (g1p <= g2p && isdigit(*g1p)){
        j = 10*j + *g1p;
        g1p ++;
      }
      if (j > max_suffix) max_suffix = j;
      if (g1p<g2p) /* must be a ^-1 at the end */
        j *= -1;
      prog2user[i] = j;
    }
    user2prog = vzalloc2(int,max_suffix+1);
    for (i=1;i<=num_gens;i++){
      j = prog2user[i];
      if (j>0) user2prog[j] = i;
    }
  }
}

static void
insert_gen(h,wp)
  gen h;
  word * wp;
{
  word * ucpy;
  gen * icpy;
  int i;
  int size = gen_array_size;
  num_gens++;
  if (num_gens >= gen_array_size) 
  /* need more space */
    gen_array_size *= 2;
  ucpy = vzalloc2(word,gen_array_size);
  icpy = vzalloc2(gen,gen_array_size);
  for (i=0;i<gen_array_size;++i) 
    word_init(ucpy+i);
  for (i=0;i<h;++i){ 
    word_cpy(user_gen_name+i,ucpy+i);
  }
  word_cpy(wp,ucpy+h);
  for (i=h;i<=num_gens-1;++i) 
    word_cpy(user_gen_name+i,ucpy+i+1);
  for (i=0;i<size;++i)
    word_clear(user_gen_name+i);
  Free_dp((dp)user_gen_name); user_gen_name=0;
  user_gen_name = ucpy;
  ucpy=0;
  for (i=0;i<h;++i){
    if (inv_of[i]>=h) icpy[i] = inv_of[i]+1;
    else icpy[i]=inv_of[i];
  }
  for (i=h;i<num_gens;i++){
    if (inv_of[i]>=h) icpy[i+1] = inv_of[i]+1;
    else icpy[i+1] = inv_of[i];
  }
    Free_dp((dp)inv_of);
  inv_of = icpy;
  icpy = 0;
}

void
delete_gen(h)
  gen h;
{
  word * ucpy;
  gen * icpy;
  int i;
  num_gens--;
  ucpy = vzalloc2(word,gen_array_size);
  icpy = vzalloc2(gen,gen_array_size);
  for (i=0;i<gen_array_size;++i) 
    word_init(ucpy+i);
  for (i=0;i<h;++i){ 
    word_cpy(user_gen_name+i,ucpy+i);
  }
  for (i=h+1;i<=num_gens+1;++i) 
    word_cpy(user_gen_name+i,ucpy+i-1);
  for (i=0;i<gen_array_size;++i)
    word_clear(user_gen_name+i);
  Free_dp((dp)user_gen_name); user_gen_name=0;
  user_gen_name = ucpy;
  ucpy=0;
  for (i=0;i<h;++i){
    if (inv_of[i]>h) icpy[i] = inv_of[i]-1;
    else icpy[i]=inv_of[i];
  }
  for (i=h+1;i<num_gens;i++){
    if (inv_of[i]>h) icpy[i-1] = inv_of[i]-1;
    else icpy[i-1] = inv_of[i];
  }
    Free_dp((dp)inv_of);
  inv_of = icpy;
  icpy = 0;
}

void
default_inverse_array()
{
  word w;
  word_traverser wt;
  gen g,h,k;
  int i,j;
  gen * g1p, * g2p;
  gen l=INVALID_GEN;
  word * ucpy;
  word inv_suffix;
  if (no_inverses){
/* no generators exist yet with "^-1" at the end, so we can simply add
all of them, interlacing them with the originally declared generators */
    if (2*num_gens >= gen_array_size){ 
      gen_array_size *= 2; 
      if (inv_of) Free_dp((dp)inv_of); 
    }
    if (inv_of==0)
      inv_of=vzalloc2(gen,gen_array_size);
    ucpy = vzalloc2(word,gen_array_size);
    word_init(&inv_suffix);
    word_put_last(&inv_suffix,'^');
    word_put_last(&inv_suffix,'-');
    word_put_last(&inv_suffix,'1');
    for (g=1;g<=num_gens;g++){
      word_init(ucpy+2*g-1);
      word_cpy(user_gen_name+g,ucpy+2*g-1);
      word_init(ucpy+2*g);
      word_cpy(user_gen_name+g,ucpy+2*g);
      word_clear(user_gen_name+g);
      word_append(ucpy + 2*g,&inv_suffix);
      inv_of[2*g-1] = 2*g;
      inv_of[2*g] = 2*g-1;
    }
    word_clear(&inv_suffix);
    num_gens *= 2;
    Free_dp((dp)user_gen_name);
    user_gen_name = ucpy;
  }
  else {
    if (inv_of==0)
      inv_of=vzalloc2(gen,gen_array_size);
    word_init(&w);
    for (g = 1; g <= num_gens; ++g) { 
      if (inv_of[g]!=0)
        continue;
    /* we may have set some inverses already explicitly */
      if (word_length(user_gen_name+g)==1){
  /* we're just checking here to see that the user hasn't assumed the old
  default of case change. Therefore if the lower and upper case versions of
  the same character appear, the program exits with a bad data message. */
        (void)word_get_first(user_gen_name+g,&h);
        if (islower(h))
          k=toupper(h);
        else
          k=tolower(h);
        gen2prog_gen(k,&l);
        if (l!=INVALID_GEN){
  /* h and k are lower and upper case versions of the same alphabet character.
  We don't allow this unless the generators are inverse to each other. This
  isn't the case here. */
          fprintf(stderr,
    "%c and %c aren't allowed as generator names\n",h,k);
          fprintf(stderr,
    "unless they're defined to be inverse to each other.\n");
          bad_data();
        }
      }
  /* now we'll set the inverse of g appropriately */
      word_traverser_init(&wt,user_gen_name+g);
      while (word_next(&wt,&h)){
        if (h=='^')
          break;
        else
          word_put_last(&w,h);
      }
      word_traverser_clear(&wt);
      if (h!='^'){
        word_put_last(&w,'^');
        word_put_last(&w,'-');
        word_put_last(&w,'1');
      }
      word2prog_gen(&w,&k);
      if (k==INVALID_GEN){
  /* There's no generator yet defined with user name w, so we have to slot one
  in, just after g.  */
        insert_gen(g+1,&w);
        k = g+1;
      }
      inv_of[g]=k;
      inv_of[k]=g;
      word_reset(&w);
    }
    word_clear(&w);
  }
  if (user_gen_prefix && prog2user==0 ){
    prog2user = vzalloc2(int,gen_array_size);
    word_init(&w);
    for (i=1;i<=num_gens;i++){ 
      g1p = user_gen_name[i].g + user_gen_name[i].first + prefix_len;
      g2p = user_gen_name[i].g + user_gen_name[i].last;
      j = 0;
      while (g1p <= g2p && isdigit(*g1p)){
        j = 10*j + (*g1p - '0');
        g1p ++;
      }
      if (j > max_suffix) max_suffix = j;
      if (g1p<g2p) /* must be a ^-1 at the end */
        j *= -1;
      prog2user[i] = j;
    }
/*
    for (i=1;i<=num_gens;i++)
*/
    user2prog = vzalloc2(int,max_suffix+1);
    for (i=1;i<=num_gens;i++){
      j = prog2user[i];
      if (j>0) user2prog[j] = i;
    }
  }
} 

boolean
read_next_rel(relp,rfile)
  word * relp;
  FILE * rfile;
{
  boolean ans=FALSE;
  word uw1;
  word_traverser wt;
  gen g;
  int posn; /* position in the file. If we get a string of the
form a = b = c .... we shall need to read b twice, once as the right hand side 
of a=b, once as the left hand side of b=c, since relations are stored as 
relators. */
  word_init(&uw1);
  ans = read_next_word(&uw1,rfile);
  if (ans){
    int c;
    word2prog_word(&uw1,relp);
    while ((c=read_char(rfile))==' ');
    if (c=='=' || c=='>'){
      int count=0;
      word uw2, pw2;
      posn = ftell(rfile); /* mark the position of the beginning of the second
in the chain, in case we need to recover it */
      word_init(&uw2);
      do {
        read_next_word(&uw2,rfile);
        count++;
        while ((c=read_char(rfile))==' ');
      }
      while (c=='=' || c=='>'); 
      ungetc(c,rfile); 
/* uw2 is now the word at the end of the chain of '='s or '>s */
      word_init(&pw2);
      word2prog_word(&uw2,&pw2);
      word_inv(&pw2,&pw2);
      word_append(relp,&pw2);
      word_clear(&uw2);
      word_clear(&pw2);
      if (count>1) fseek(rfile,posn,0); 
/* if there's more than one = or > then go back to the beginning of the 
second word in the chain it'll be the lhs of another relation */
    }
    else
      ungetc(c,rfile);
  }
  word_clear(&uw1);
  word_traverser_init(&wt,relp);
  while (word_next(&wt,&g)){
    if (g==INVALID_GEN){
      fprintf(stderr,"Undeclared generator or syntax error in relator.\n");
      user_word_print(stderr,relp);
      fprintf(stderr,"\n");
      bad_data();
    }
  }
  word_traverser_clear(&wt);
  return ans;
}


void
word_factor(wp,wwp,ep)
  word * wp,* wwp;
  int  * ep;
{
  int length = word_length(wp);
  int baselength = 0;
  word power;
  word_init(&power);
  word_reset(wwp);
  if (length==0) return;
  while (baselength <= length){
    word_traverser wt;
    int count = 0;
    gen g;
    int i;
    word_traverser_init(&wt,wp);
    while (word_next(&wt,&g)){
      count++;
      if (count<=baselength) continue;
      word_put_last(wwp,g);
      baselength++;
      if (baselength>length/2 ||length%baselength==0)
        break;
    }  
    word_traverser_clear(&wt);
    if (baselength > length/2){
      baselength = length;
      *ep = 1;
      break;
    }  
    else {
      *ep = length/baselength;
      for (i=1;i<=*ep;i++)
        word_append(&power,wwp);
      if (word_sgn(wp,&power)==0)
        break;
    }  
    word_reset(&power);
  }
  word_clear(&power);
  if (baselength == length)
    word_cpy(wp,wwp);
}

char * word2string(wp)
  word * wp;
{
  char *  string;
  int len = word_length(wp);
  int j,k;
  string = vzalloc2(char,len+1);
  j = 0; k= wp->first;
  while (j<len) string[j++] = (wp->g)[k++];
  string[j] = '\0';
  return string;
}

  

