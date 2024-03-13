/* isom.d/lib file permfns.c */
#include <stdio.h>
#include <ctype.h>
#include "defs.h"
#define MAXINT 2147483647

extern  int npt,**cp,**lcp,**ucp,**pptr,*pno,stop;


image(pt) int pt;
/* This assumes that we have a permutation p stored as a product of permutations
   pointed at by lcp, lcp+1,...,ucp.
   The image of pt under p is calculated.
*/
{ int **p;
  p=lcp-1;
  while (++p<=ucp) pt=(*p)[pt];
  return(pt);
}

addsvb(pt,sv) int pt,**sv;
/* Again there is already a list of permutations lcp,lcp+1,...,ucp.
   sv is a Schreier vector and pt is a point.
   We use sv to get a product of permutations which takes the base point to
   pt, and insert this product at the beginning of the list.
*/
{ int *p;
  while ((p=sv[pt])!= &stop) {pt=p[pt]; lcp--; *lcp =p-(npt+1); }
  return;
}

addsvf(pt,sv) int pt,**sv;
/* Again there is already a list of permutations lcp,lcp+1,...,ucp.
   sv is a Schreier vector and pt is a point.
   We use sv to get a product of permutations which takes pt to the base point
   and insert this product at the end of the list.
*/
{ int *p;
  while ((p=sv[pt])!= &stop) {pt=p[pt]; ucp++; *ucp=p; }
  return;
}

exprep(pt,ptr,sv) int pt,*ptr,**sv;
/* The word for pt is computed using Schreier vector sv, and the corresponding
   permutation that takes the base point to pt is  stored in ptr.
   The inverse of the perm is stored in ptr+npt
   Note that it is the inverse that is actually calculated by addsvf.
*/
{ int i;
  ucp=lcp-1; addsvf(pt,sv);
  for (i=1;i<=npt;i++) ptr[npt+i]=image(i);
  invert(npt,ptr+npt,ptr);
  return(0);
}

sgs(base,anb,lorb,svptr,anp2,actgen,baseknown,ord_stab)
   int *base,*anb,*lorb,***svptr,*anp2,*actgen,ord_stab;
   char baseknown;
/* 
   This is the Schreier-Sims algorithm for computing a base and strong
   generating set.
   The initial generators are assumed to be stored in
   pptr[0],pptr[1],...,pptr[*anp2],
   where pptr[2n+1] is the inverse of pptr[2n].
   Any initial base points are in base[1],...,base[nb] (where *anb=nb).
   The basic orbit lengths will be calculated as lorb[1],...,lorb[nb],
   and the Schreier vectors as svptr[1],...svptr[nb].
   (Each Schreier vector is a list of npt pointers to permutations.)
   actgen[i] should normally be initialized to 1 for 1<=i<= *anp2 before
   calling this procedure for the first time on a group.
   (If further generators are added to the group later, it may be
    called again, and then actgen[i] should be put 1 on the new generators.)
   baseknown=1 means base[1],...,base[nb] is already known to be a base.
   baseknown=2 means in addition that generators pptr[0],..., are already
   known to be a strong generating set
   ord_stab>0 means that the order of the stabiliser of the first base
   point can be assumed to be ord_stab.
*/
/* Externals pptr, pno, cp, lcp, ucp, npt */
{ int i,j,k,l,m,bno,u,v,**w,**x,y,z,**lsv,np2,nb,*p,npt1,*orb;
  char trivrel,id,got_ord_stab;
  npt1=npt+1;
  tmalloc(orb,int,npt1)
  nb= *anb;
  np2= *anp2;
/* First set pptr[i][npt+1] for each even value of i. This is set equal to
  j, where the permutation lies in P(j)-P(j+1) i.e. it fixes the first
  j-1 base points but not the j-th.
  If any permutation is found fixing all base points, then add a new base point.
*/
  for (i=0;i<np2;i+=2)
  { pptr[i][npt1]=0;
    for (j=1;j<=nb;j++) if (pptr[i][base[j]]!=base[j])
    {pptr[i][npt1]=j; break;}
    if ((baseknown==0 && pptr[i][npt1]==0) || nb==0)
    { for (k=1;k<=npt;k++)
        if (pptr[i][k]!=k)
        { nb++; base[nb]=k;
          pptr[i][npt1]=nb;
          tmalloc(svptr[nb],int *,npt1)
          break;
        }
    }
  }

/* Just in case all generators are trivial: */
  if (nb==0)
  { nb=1; base[1]=1; tmalloc(svptr[nb],int *,npt1)}
  for (i=1;i<=nb;i++) lorb[i]=1;
  got_ord_stab=0;

/* Now we start the main loop. We assign space for a new candidate strong
   generator and its inverse.
   bno is the number for which we are currently calculating the
   stabilizer of base[bno] in P(bno). This stabilizer should be equal to
   P(bno+1). If not, we add it as a new generator of P(bno+1).
*/
  bno=nb;
  tmalloc(pptr[np2],int,2*npt1) pptr[np2+1]=pptr[np2]+npt1;

loop:
/* We start by making a list pno[1],pno[2],...,pno[n] (where n=pno[0])
   of those even numbers 2n such that pptr[2n] lies in P(bno).
  (i.e. it fixes base[1],...,base[bno-1])
  We only take those for which actgen[i]<=bno. This is because those
  generators with actgen[i]>bno are known to be redundant as generators
  of P(bno).
*/
  *pno=0;
  for (i=0;i<np2;i+=2)
  { if (pptr[i][npt1]>=bno && actgen[i]<=bno)
    { (*pno)++; pno[*pno]=i; }
  }
/* Calculate orbit and Schreier vector */
  lorb[bno]=orbitsv(base[bno],svptr[bno],orb);
  if (ord_stab>0 && !got_ord_stab) {
  /* Check if our current order of point stabiliser is ord_stab */
    char overflow=0; int cord_stab=1;
    for (i=2;i<=nb;i++) {
      if (cord_stab>=MAXINT/lorb[i]) overflow=1;
      else cord_stab *= lorb[i];
    }
    if (overflow || cord_stab>ord_stab) {
      fprintf(stderr,
         "WARNING: Order given as stabiliser order is too small.\n");
      ord_stab=0;
    }
    else if (cord_stab==ord_stab)
      got_ord_stab=1;
      /*printf("ord_stab=%d, cord_stab=%d, nb=%d\n",ord_stab,cord_stab,nb);*/
  }
  if (*pno!=0 && baseknown<2 && !got_ord_stab)
  { y=np2+1;
/* Now start calculating the Schreier generators of the stabilizer of base[bno]
   These Schreier generators will be stored as a list of pointers to
   permutations, lcp,lcp+1,...,ucp.
*/
    i= bno==1 && ord_stab>0 ? abs((rand()*rand())%lorb[bno]) + 1 : 1;
    while (i<=lorb[bno])
    { ucp=lcp-1; addsvf(orb[i],svptr[bno]);
      for (w=lcp,x=ucp;w<=x;w++,x--)
      { if (w==x) (*w)-=npt1; else {p= *w; *w= *x-npt1; *x=p-npt1;}}
      lsv= ucp;
      for (j=1;j<=*pno;j++)
      { ucp=lsv;
        trivrel = (ucp>=lcp) ? *ucp==pptr[pno[j]] : 0;
/* if trivrel is true, then the Schreier generator must be trivial so we
   skip over it
*/
        if (trivrel==0)
        { *(++ucp)=pptr[pno[j]+1]; id=1;
/* Now we check to see if this Schreier generator is in our current P(bno+1)
*/
          for (l=bno;l<=nb;l++)
          { v=base[l]; u=image(v); if (svptr[l][u]==0) {id=0; break;}
            addsvf(u,svptr[l]);
          }
          if (id==0)
/* It isn't so we have to add on the stripped part as a new generator
   (of P(l)).
*/
          { pptr[np2][npt1]=l; actgen[np2]=bno+1;
/* This records the fact that the new generator is redundant as a generator
   of P(j) for any j>bno
*/
            for (m=1;m<=npt;m++) {u=image(m); pptr[np2][m]=u; pptr[y][u]=m;}
            np2+=2;
            tmalloc(pptr[np2],int,2*npt1) pptr[np2+1]=pptr[np2]+npt1;
            bno=l; goto loop;
          }
          if (baseknown==0)
/*We still haven't finished the check. We still have to check that the
 stripped Schreier genrator is equal to the identity */
          { for (l=1;l<=npt;l++) if (image(l)!=l)
/* It isn't, so we have a new base point as well as a new strong generator. */
            { nb++; base[nb]=l; tmalloc(svptr[nb],int *,npt1);
              pptr[np2][npt1]=nb; actgen[np2]=bno+1;
              for (m=1;m<=npt;m++) {u=image(m); pptr[np2][m]=u; pptr[y][u]=m;}
              np2+=2;
              tmalloc(pptr[np2],int,2*npt1) pptr[np2+1]=pptr[np2]+npt1;
              bno=nb;
              goto loop;
            }
          }
        }
      }  
      i= bno==1 && ord_stab>0 ? abs((rand()*rand())%lorb[bno])+1 : i+1;
    }  
  }
  bno--;

  if (bno==0)
  { *anp2=np2;
    *anb=nb;
    tfree(pptr[np2])
    tfree(orb)
    return(0);
  }
  else goto loop;
}

ingp(base,nb,svptr,p) int *base,nb,***svptr, *p;
/*Check whether the permutation p is in the permutation group with
  Schreier vectors svgptr[1],...,svgptr[nb]
  Externals: cp, lcp, ucp. 
*/
{ int i,im;
  ucp=lcp; *ucp=p;
  for (i=1;i<=nb;i++)
  { im=image(base[i]);
    if (svptr[i][im]==0) return(0);
    addsvf(im,svptr[i]);
  }
  return(1);
}

orbitsv(pt,svec,orb) int pt,**svec,*orb;
/*  Calculate the orbit of the point pt under the permutations
  pptr[pno[1]], pptr[pno[2]], ..., pptr[pno[n]] (where n=pno[0]),
  and calculate its Schreier vector as sv.
  Externals: npt,stop,pno,pptr.
*/
{ int u,v,w,x,y,z,lo;
  for (u=1;u<=npt;u++) svec[u]=0;
  orb[1]=pt; lo=1; svec[pt]= &stop;
 
  for (x=1;x<=lo;x++)
  { z=orb[x];
    for (y=1;y<= *pno;y++)
    { w=pno[y]; v=pptr[w][z];
      if (svec[v]==0) { lo++; orb[lo]=v; svec[v]=pptr[w+1]; }
    }
  }  
  return(lo);
}
