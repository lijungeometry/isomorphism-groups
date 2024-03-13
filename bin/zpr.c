/* ========================== C MeatAxe =============================
   zpr.c - Print a matrix or permutaion.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zpr.c,v 2.9 1994/03/13 12:37:35 mringe Exp $
 *
 * $Log: zpr.c,v $
 * Revision 2.9  1994/03/13  12:37:35  mringe
 * Maschinenunabhaengiges Dateiformat (Perm)
 *
 * Revision 2.8  1994/02/15  14:09:19  mringe
 * zsetlen() erweitert fuer Permutationen.
 *
 * Revision 2.7  1994/02/15  13:28:01  mringe
 * Benutze zseek() statt fseek().
 *
 * Revision 2.6  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.5  1993/12/08  12:00:41  mringe
 * malloc() prototype.
 *
 * Revision 2.4  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.3  1993/10/26  10:47:35  mringe
 * Compiler Warnings.
 *
 * Revision 2.2  1993/10/21  21:57:35  mringe
 * Permutationen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.17  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.16  1993/10/06  06:09:00  mringe
 * Option -s.
 *
 * Revision 1.15  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.14  1993/08/10  13:46:18  mringe
 * Ueberfluessige klammern beim Perm.-Output weglassen.
 *
 * Revision 1.13  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.12  1993/07/19  17:12:36  mringe
 * Option -G: Benutze MeatAxe-Record.
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.10  1993/06/28  18:13:00  mringe
 * GAP Perm.utationen: Benutze Filenamen fuer die Elemente.
 *
 * Revision 1.9  1993/03/04  16:24:27  mringe
 * GAP-Output verbessert (maximal 80 Zeichen/Zeile).
 *
 * Revision 1.8  1993/03/04  15:53:38  mringe
 * GAP-Output von Permutationen beschleunigt.
 * Benutze utils-Library fuer help().
 *
 * Revision 1.7  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.6  1992/07/22  14:06:38  mringe
 * Use zgen
 *
 * Revision 1.5  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.4  1992/05/26  18:48:37  mringe
 * help() eingebaut.
 *
 * Revision 1.3  1992/05/26  18:07:05  mringe
 * Ausgabe von Permutationen im GAP-Format: Zeilenl"ange korrigiert.

 * Revision 1.2  1992/05/26  17:59:05  mringe
 * Bei Permutationen im GAP-Format keine Fixpunkte ausgeben,
 * da GAP diese am Anfang einer Permutation nicht einlesen
 * kann.
 *
 * Revision 1.1  1992/05/26  17:41:49  mringe
 * Initial revision
 *
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "meataxe.h"

/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */


static void err _PL((int c));
static void prmatrix _PL((void));
static void prgap _PL((void));
static void prperm _PL((void));
static void prmtx _PL((void));
static void setrange _PL((char *c));
static void printsummary _PL((long fl, long nor, long noc));
int main _PL((int argc, char *argv[]));



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long fl, nor, noc;
static FILE *dest = NULL;	/* Output file */
static char *inpname = "stdin";
static FILE *inpfile = NULL;
static long first = 1;		/* Output range (option -r) */
static long last = 10000000;

static char *helptext[] = {
"SYNTAX",
"    zpr [-Gs] [-r <Range>] [<Binfile> [<Textfile>]]",
"",
"OPTIONS",
"    -G   GAP output",
"    -r   Print a range of rows (or permutations). <Range> may be a",
"         single number, or of the form <First>-<Last>.",
"    -s   Print summary only\n",
"",
"FILES",
"    <Binfile>   i  A matrix or permutation in binary format",
"                   Default name: G1",
"    <Textfile>  i  The output in text format",
"                   Default name: T1",
"",
NULL};

static proginfo_t pinfo =
   { "zpr", "Print Permutations Or Matrices", "$Revision: 2.9 $",
      helptext };



/* ------------------------------------------------------------------
   err() - Print error message and exit.
   ------------------------------------------------------------------ */

static void err(c)
int c;

{   fprintf(stderr,"ZPR ERROR - ");
    switch (c)
    {
	case 'f':
	    fprintf(stderr,"FILE I/O ERROR\n");
	    break;
	case 't':
	    fprintf(stderr,"CANNOT PRINT THIS FILE TYPE\n");
	    break;
    }
    exit(EXIT_ERR);
}



/* ------------------------------------------------------------------
   prmatrix() - Print a matrix in standard format.
   ------------------------------------------------------------------ */

static void prmatrix()

{	PTR m1;
	FEL f1;
	long loop1, j1;
	int md, mx, iv;

	if (last > nor) last = nor;
	if (first > last) first = last+1;
	zsetlen(fl,noc);
	m1 = zalloc((long)1);
	md = 4;
	mx = 5;
	if (fl < 100) {md = 3; mx = 25;}
	if (fl < 10) {md = 1; mx = 80;}
	fprintf(dest,"%2d%6ld%6ld%6ld\n",md,fl,last-first+1,noc);
	for (loop1 = 1; loop1 <= last; ++loop1)
	{	if (zreadvec(inpfile,m1,1) != 1) errexit(-1,inpname);
		if (loop1 < first) continue;
		iv = 1;
		for (j1 = 1; j1 <= noc; ++j1)
		{	f1 = zextract(m1,j1);
			switch (md)
			{	case 1:
					fprintf(dest,"%1ld",zftoi(f1));
					break;
				case 3:
					fprintf(dest,"%3ld",zftoi(f1));
					break;
				case 4:
					fprintf(dest,"%8ld",zftoi(f1));
					break;
			}
			if (iv++ >= mx)
			{	fprintf(dest,"\n");
				iv = 1;
			}
		}
		if (iv > 1)
			fprintf(dest,"\n");
	}
}


/* ------------------------------------------------------------------
   prgapmat() - Print a matrix in GAP format.
   ------------------------------------------------------------------ */

static void prgapmat()

{   PTR m1;
    FEL f1;
    FEL gen;
    long loop1, j1;
    int cnt, isprimefield;


    if (last > nor) last = nor;
    if (first > last) first = last+1;

    zsetlen(fl,noc);
    gen = zgen;		/* Generator */
    isprimefield = 1;
    for (j1 = 2; isprimefield && j1 < fl; ++j1)
		if (fl % j1 == 0) isprimefield = 0;
    m1 = zalloc((long)1);
    fprintf(dest,"MeatAxe.Matrix := [\n");
    zseek(inpfile,first);
    for (loop1 = first; loop1 <= last; ++loop1)
    {	
	if (zreadvec(inpfile,m1,1) != 1) errexit(-1,inpname);
	cnt = 0;
	fprintf(dest,"[");
	for (j1 = 1; j1 <= noc; ++j1)
	{   if (cnt > 75)
	    {	fprintf(dest,"\n ");
		cnt = 0;
	    }
	    f1 = zextract(m1,j1);
	    if (isprimefield)
	    {   FEL f2=F_ZERO;
	   	long k=0;
	    	while (f2 != f1)
	    	{   f2 = zadd(f2,gen);
		    ++k;
		}
		fprintf(dest,"%ld",k);
		cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
	    }
	    else
	    {   if (f1 == F_ZERO)
	    	{   fprintf(dest,"0*Z(%ld)",fl);
		    cnt += 5;
		    cnt += fl>9999?5:fl>999?4:fl>99?3:fl>9?2:1;
		}
		else
		{   FEL f2 = gen;
		    long k = 1;
		    while (f2 != f1)
		    {   f2 = zmul(f2,gen);
			++k;
		    }
		    fprintf(dest,"Z(%ld)^%ld",fl,k);
		    cnt += 4;
		    cnt += fl>9999?5:fl>999?4:fl>99?3:fl>9?2:1;
		    cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
		}
	    }
	    if (j1 < noc)
	    {	fprintf(dest,",");
		++cnt;
	    }
	}
	fprintf(dest,"]");
	if (loop1 < last)
	    fprintf(dest,",");
	fprintf(dest,"\n");
    }
    fprintf(dest,"]");
    if (isprimefield)
	fprintf(dest,"*Z(%ld)",fl);
    fprintf(dest,";\n");
}


/* ------------------------------------------------------------------
   prgapperm() - Print a permutation in GAP format.
   ------------------------------------------------------------------ */

#define SIZE(i) ((i)<9?2 : (i)<99?3 : (i)<999?4 : (i)<9999?5 : \
	(i)<99999?6 : (i)>999999?7 : (i)<9999999?8 : 9)


static void prgapperm()

{
    int count;
    long i, k, *p, pos;
    long cycle;
    long *perm;
    char name[250], *c1, *c2;

    /* Make name
       --------- */
    c1 = inpname;
    c2 = name;
    while (isalnum(*c1)) *c2++ = *c1++;
    *c2 = 0;

    if (last > noc) last = noc;
    if (first > last) first = last+1;
    if (first > noc) first = 1;
    perm = (long *) malloc(sizeof(long) * nor);
    if (perm == NULL) errexit(ERR_NOMEM,"zpr: prgapperm()");

    zsetlen((long)-1,nor);
    zseek(inpfile,first);
    fprintf(dest,"MeatAxe.Perms := [\n");
    for (pos = first; pos <= last; ++pos)
    {
	if (zreadlong(inpfile,perm,nor) != nor)
	    errexit(ERR_FILEREAD,inpname);
	fprintf(dest,"    ");
	count = 4 + SIZE(pos);
	p = (long *) perm;
	cycle = 0;
	while (1)
	{
	    /* Find next cycle
	       --------------- */
	    while (cycle < nor && p[cycle] == 0) ++cycle;
	    if (cycle >= nor) break;	/* Done */
	    i = cycle;

	     /* Check if it is a fixed point (we don't
	        print fixed points, GAP doesn't like them)
	        ----------------------------------------- */
	    if (p[i] == i+1)
	    {
		p[i] = 0;	/* Mark it as done */
		continue;
	    }

	    /* Print cycle
	       ----------- */
	    if ((count += SIZE(i)) > 77)
	    {
		fprintf(dest,"\n    (%ld",i+1);
	        count = 5 + SIZE(i);
	    }
	    else
		fprintf(dest,"(%ld",i+1);

	    while (1)
	    {
		k = i;
		i = p[i] - 1;
		p[k] = 0;
		if (p[i] == 0) break;

	        if ((count += SIZE(i)) > 77)
	        {
		    fprintf(dest,",\n     %ld",i+1);
	            count = 4 + SIZE(i);
	        }
	        else
		    fprintf(dest,",%ld",i+1);
	    }
	    fprintf(dest,")");
	    ++count;
	}
	if (pos < last) fprintf(dest,",");
	fprintf(dest,"\n");
    }
    fprintf(dest,"];\n");
}


/* ------------------------------------------------------------------
   prgap() - Print a matrix or permutation in GAP format.
   ------------------------------------------------------------------ */

static void prgap()

{	if (fl == -1)
		prgapperm();
	else if (fl >= 2)
		prgapmat();
	else
		err('t');
}


/* ------------------------------------------------------------------
   prperm() - Print a permutation in standard format.
   ------------------------------------------------------------------ */

static void prperm()

{
    long mono, f1, i, count;
    long *perm;

	if (last > noc) last = noc;
	if (first > last) first = last+1;
	mono = -fl;
    perm = (long *) malloc(sizeof(long) * nor);
    if (perm == NULL) errexit(ERR_NOMEM,"zpr: prperm()");

    fprintf(dest,"%2d%6ld%6ld%6ld\n",(mono == 1) ? 12 : 13,
	mono,nor,last-first+1);
    zsetlen((long)-1,nor);
    zseek(inpfile,first);
    for (count = first; count <= last; ++count)
    {
	if (zreadlong(inpfile,perm,nor) != nor)
	    errexit(ERR_FILEREAD,inpname);
	for (i = 0; i < nor; ++i)
	{
	    f1 = perm[i];
	    if (mono == 1)
		fprintf(dest,"%6ld\n",f1);
	    else
		fprintf(dest,"%6ld%4ld\n",(f1-1)/mono+1,(f1-1)%mono);
	}
    }
}



/* ------------------------------------------------------------------
   prmtx() - Print a matrix or permutation in standard format.
   ------------------------------------------------------------------ */

static void prmtx()

{	if (fl > 1)
		prmatrix();
	else if (fl <= -1)
		prperm();
	else
		err('t');
}


/* ------------------------------------------------------------------
   setrange() - Set output range
   ------------------------------------------------------------------ */

static void setrange(c)
char *c;

{	
    if ((first = getint()) == GETINT_ERR) errexit(ERR_OPTION,"-r");
    if (*opt_text_ptr == 0)
    {	last = first;
    }
    else if (*opt_text_ptr == '-')
    {   ++opt_text_ptr;
        if ((last = getint()) == GETINT_ERR) 
    	    errexit(ERR_OPTION,"-r");
    }
    if (*opt_text_ptr != 0 || first < 1 || last < first)
       	errexit(ERR_OPTION,"-r");
}

/* ------------------------------------------------------------------
   printsummary()
   ------------------------------------------------------------------ */

static void printsummary(fl,nor,noc)
long fl, nor, noc;

{
    printf("%s: ",inpname);
    if (fl == -1)
	printf("%ld Permutation%s of degree %ld\n",noc,
		noc == 1 ? "" : "s", nor);
    else if (fl >= 2)
	printf("%ld x %ld matrix over GF(%ld)\n",nor,noc,fl);
    else
	printf("Unknown file format.\n");

}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    enum {MTX,GAP,SUMMARY} mode = MTX;
    int i;

    dest = stdout;       /* Output file */
    inpfile = stdin;


    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gsr")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': mode = GAP; msg_level = -100; break;
	    case 's': mode = SUMMARY; break;
	    case 'r': setrange(opt_text); break;
	}
    }
    switch (argc - opt_ind)
    {   case 0:
	    if ((dest = os_fopen("T1",FM_CREATE|FM_TEXT)) == NULL)
	    {	perror("T1");
		err('f');
	    }
	    break;
    	case 2:
	    if (strcmp(argv[opt_ind+1],"-"))
		if ((dest = os_fopen(argv[opt_ind+1],FM_CREATE|FM_TEXT))
							== NULL)
		{	perror(argv[opt_ind+1]);
			err('f');
		}
	    /* no break here! */
        case 1:
	    inpname = argv[opt_ind];
	    if ((inpfile = zreadhdr(inpname,&fl,&nor,&noc)) == NULL)
		errexit(-1,inpname);
	    break;
    	default:
	    errexit(ERR_BADUSAGE,"zpr");
    }
    switch (mode)
    {	case MTX:
	    prmtx();
	    break;
	case GAP:
	    prgap();
	    break;
	case SUMMARY:
	    printsummary(fl,nor,noc);
	    break;
    }
    fclose(inpfile);
    fclose(dest);
    return (EXIT_OK);
}


