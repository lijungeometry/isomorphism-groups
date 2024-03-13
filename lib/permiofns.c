/* isom.d/lib file permiofns.c */
#include <stdio.h>
#include <ctype.h>
#include "defs.h"

readperm(npt,ptr,rfile)  int npt,*ptr; FILE *rfile;
/* The next permutation on npt points from rfile is read into array ptr.
   If the data is invalid (if it is not a permutation for example), we
   exit(2).
   The permutation may either be in cyclic form or as a list of images.
   It may optionally be preceded by one or more labels, which will be
   ignored. However, the word "identity" is interpreted as the
   identity permutation, so is not valid as a label.
*/
{ int i,c,l; char *had,string[10];
  tmalloc(had,char,npt+1)
  for (i=1;i<=npt;i++)  had[i]=0;
  c=read_char(rfile);
  while (c!='(' && !isdigit(c) && !isalpha(c)) c=read_char(rfile);
  while (isalpha(c))
/* If the permutation has a label, then we want to ignore it */
  { ungetc(c,rfile);
    read_next_string(string,9,rfile);
    if (strcmp(string,"identity ")==0)
    { for (i=1;i<=npt;i++) ptr[i]=i;
      return(0);
    }
    c=read_char(rfile);
    while (c!='(' && !isdigit(c) && !isalpha(c)) c=read_char(rfile);
  }

  if (isdigit(c))
/* The permutation will be given as a list of images */
  { ungetc(c,rfile);
    for (i=1;i<=npt;i++)
    { read_next_int(&l,rfile);
      if (l<=0 || l>npt || had[l])
      { fprintf(stderr,"Invalid or repeated point in permutation %d\n",l);
        exit(2);
      }
      had[l]=1;
      ptr[i]=l;
    }
    return(0);
  }
  else
/* The permutation is in cyclic form. Any character other than whitespace
   (spaces, newlines, tabs) or "(" will terminate the permutation.
 */
  { char cyc; int j,k;
    for (i=1;i<=npt;i++) {ptr[i]=i; had[i]=0;}
    c=read_char(rfile);
/* An immediate closing bracket means the identity permutation, as in GAP */
    if (c!=')')
    { ungetc(c,rfile);
      while(1)
      { cyc=1; j=0;
        { cyc=1; j=0;
          while (cyc)
          { read_next_int (&l,rfile);
            if (l<=0 || l>npt || had[l])
            { fprintf(stderr,"Invalid or repeated point in permutation%d\n",l);
              exit(2);
            }
            had[l]=1;
            if (j==0) {j=l;k=l;} else {ptr[k]=l;k=l;}
            while ((c=read_char(rfile))==' ');
            if (c==')') {cyc=0; ptr[k]=j;}
          }  
        }
        while ((c=read_char(rfile))==' ');
        if (c!='(')
        { ungetc(c,rfile);
          break;
        }
      }
    }
    return(0);
  }
}

readbaselo(anb,abase,alorb,rfile) int *anb,**abase,**alorb; FILE *rfile;
/* The nb base points are read into base, and the nb orbit lengths into lorb,
   from rfile.
   Here base = *abase, lorb = *alorb and  nb = *nb.
*/
{ int i,nb,*base,*lorb;
  char * string,ans;
  tmalloc(string,char,20);
  *string=0;
  while (strcmp(string,"base_pts"))
  { ans=read_next_string(string,8,rfile);
    if (ans==0)
    { fprintf(stderr,"Could not find keyword 'base_pts'.\n"); return(-1);}
  }
  read_next_int(anb,rfile); nb = *anb;
  tmalloc(*abase,int,nb+1) base= *abase;
  tmalloc(*alorb,int,nb+1) lorb= *alorb;
  find_char('{',rfile);
  for (i=1;i<=nb;i++) read_next_int(base+i,rfile);
  find_char('}',rfile);
  *string=0;
  while (strcmp(string,"basic_orbit_lengths"))
  { ans =read_next_string(string,19,rfile);
    if (ans==0)
    { fprintf(stderr,"Could not find keyword 'basic_orbit_lengths'.\n");
      return(-1);
    }
  }
  find_char('{',rfile);
  for (i=1;i<=nb;i++) read_next_int(lorb+i,rfile);
  find_char('}',rfile);
  return(0);
}

readallperms(npt,nb,anperms,aperm,pptr,rfile)
  int npt,nb,*anperms,**aperm,**pptr;
  FILE *rfile;
/* nperms permutations  are read into perm  and inverted.
   nperms = *anperms, perm = *aperm.
   The permutations are packed into the space pointed at by perm. The 
   permutations themselves are pointed at by pptr[0],pptr[1],...,
   where pptr[2n] is the (n+1)-st permutation read in and pptr[2n+1] is
   its inverse.
   pptr[i][j] is the action of pptr[i] on point j, with 1<=j<=npt.
   Each permutation is in fact allotted npt+1 points.
   pptr[i][npt+1] will be used later to record its positing in the
   stabilizer chain (pptr[i][npt+1]=j will mean it lies in P(j)-P(j+1)).
   I'VE JUST NOTICED THAT nb DOESN'T NEED TO BE A PARAMETER!
*/
{ int i,j,k,*ptr,*temp,npt1,nperms,*perm;
  char * string,ans;
  npt1=npt+1;
  tmalloc(string,char,12);
  *string=0;
  while (strcmp(string,"gens       ") && strcmp(string,"strong_gens"))
  { ans=read_next_string(string,11,rfile);
    if (ans==0)
    { fprintf(stderr,"Could not find keyword 'gens' or 'strong_gens'.\n");
      return(-1);
    }
  }
  read_next_int(anperms,rfile); nperms= *anperms;
  tmalloc(*aperm,int,npt1*2*nperms); perm= *aperm;
  find_char('{',rfile);
  ptr=perm-1;
  for (i=1;i<=nperms;i++)
  { readperm(npt,ptr,rfile); invert(npt,ptr,ptr+npt1);
    pptr[2*i-2]=ptr; pptr[2*i-1]=ptr+npt1;
    ptr+= 2*npt1;
  }
  find_char('}',rfile);
  return(0);
}

printbaselo(nb,base,lorb,wfile) int nb,*base,*lorb; FILE *wfile;
/* Write the base and basic orbit lengths to wfile in the right format */
{ int i;
  fprintf(wfile," base_pts %d {",nb);
  for (i=1;i<=nb;i++) fprintf(wfile," %d",base[i]);
  fprintf(wfile,"}\n");
  fprintf(wfile," basic_orbit_lengths {");
  for (i=1;i<=nb;i++) fprintf(wfile," %d",lorb[i]);
  fprintf(wfile,"}\n");
  return(0);
}

printperm(npt,p,cycle,label,wfile) int npt,*p; char cycle, *label; FILE *wfile;
/* The permutation p is printed with label to wfile. If cycle is true it is
   printed in cycles - otherwise as point images.
  Note that label should not contain the following colon. That is added on.
  It's a headache working out where to put in carriage returns.
*/
{ int i,l;
  if ((l=strlen(label))>0)
  { fprintf(wfile,"  %s:",label);
    if (l<5) for (i=1;i<=5-l;i++) putc(' ',wfile);
  }
  if (cycle)
  { int *temp,m,im,ct; char id;
    ct=0;
    tmalloc(temp,int,npt+1)
    for (i=1;i<=npt;i++) temp[i]=1;
    id=1;
    for (m=1;m<=npt;m++) if (temp[m])
    { if ((im=p[m])!=m)
      { if ((ct=ct+2+numdigits(m))>72)
        { fprintf(wfile,"\n        "); ct=2+numdigits(m);}
        fprintf(wfile,"(%d,",m);
        if ((ct=ct+1+numdigits(im))>72)
        { fprintf(wfile,"\n        "); ct=1+numdigits(im);}
        fprintf(wfile,"%d",im);
        id=0; temp[im]=0;
        while ((im=p[im])!=m)
        { temp[im]=0;
          fprintf(wfile,",");
          if ((ct=ct+1+numdigits(im))>72)
          { fprintf(wfile,"\n        "); ct=1+numdigits(im);}
          fprintf(wfile,"%d",im);
        }
        fprintf(wfile,")");
      }
    }
    if (id) fprintf(wfile," identity\n"); else fprintf(wfile,"\n");
    tfree(temp)
  }
  else
  { int n=numdigits(npt);
    if (n==1) for (i=1;i<=npt;i++) fprintf(wfile,"%2d",p[i]);
    else if (n==2) for (i=1;i<=npt;i++)
    { if (i>1 && (i-1)%24==0) fprintf(wfile,"\n        ");
      fprintf(wfile,"%3d",p[i]);
    }
    else if (n==3) for (i=1;i<=npt;i++)
    { if (i>1 && (i-1)%18==0) fprintf(wfile,"\n        ");
      fprintf(wfile,"%4d",p[i]);
    }
    else if (n==4) for (i=1;i<=npt;i++)
    { if (i>1 && (i-1)%14==0) fprintf(wfile,"\n        ");
      fprintf(wfile,"%5d",p[i]);
    }
    else if (n==5) for (i=1;i<=npt;i++)
    { if (i>1 && (i-1)%12==0) fprintf(wfile,"\n        ");
      fprintf(wfile,"%6d",p[i]);
    }
    else return(-1);
    fprintf(wfile,"\n");
  }
  return(0);
}

invert(npt,ptr1,ptr2) int npt,*ptr1,*ptr2;
/* permutation ptr1 is inverted and put in ptr2.  */
{ int i;
  for (i=1;i<=npt;i++) ptr2[ptr1[i]]=i;
  return;
}
 
numdigits(n) int n;
{ if (n<10) return(1);
  if (n<100) return(2);
  if (n<1000) return(3);
  if (n<10000) return(4);
  if (n<100000) return(5);
  return(-1);
}

seeknln(rfile) FILE *rfile;
/* The next new line in rfile is found. */
{while (getc(rfile)!='\n'); }
