/* ========================== C MeatAxe =============================
   zcv.c - Convert a matrix or permutation from ASCII (readable)
   format to binary format.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zcv.c,v 2.7 1994/03/13 12:37:35 mringe Exp $
 *
 * $Log: zcv.c,v $
 * Revision 2.7  1994/03/13  12:37:35  mringe
 * Maschinenunabhaengiges Dateiformat (Perm)
 *
 * Revision 2.6  1994/03/08  14:05:01  mringe
 * Suche Input in MTXLIB.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.18  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.17  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.16  1993/08/31  08:02:03  mringe
 * Benutze '=' statt '-' fuer Std.Input.
 *
 * Revision 1.15  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.14  1993/07/28  07:54:47  mringe
 * Freies Format auch bei Mode 5
 *
 * Revision 1.13  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.12  1993/05/18  16:27:36  mringe
 * Bug beim Einlesen des headers korrigiert.
 *
 * Revision 1.11  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.10  1992/07/28  08:35:17  mringe
 * Behendlung negativer Zahlen in Mode 5 korrigiert
 *
 * Revision 1.9  1992/07/23  09:24:33  mringe
 * *** empty log message ***
 *
 * Revision 1.8  1992/07/23  08:39:35  mringe
 * Erlaube Kommentare im Input (`# ...')
 *
 * Revision 1.7  1992/07/22  18:22:01  mringe
 * Some minor changes.
 *
 * Revision 1.6  1992/07/22  14:05:53  mringe
 * Use zchar.
 *
 * Revision 1.5  1992/07/22  10:30:13  mringe
 * Freies Format fuer alle Matrizen, ausser Mode 1
 *
 * Revision 1.4  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.3  1992/07/08  09:38:49  mringe
 * Neu: Mode 6 = Matrix in free format.
 *
 * Revision 1.2  1992/07/07  12:55:06  mringe
 * Fehler im Format bei Mode 4 beseitigt.
 *
 * Revision 1.1  1992/05/26  18:27:13  mringe
 * Initial revision
 *
 */

#define MAXLINE 4000

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void err _PL((int c)) ;
static int readline _PL((void)) ;
static long readlong _PL((void)) ;
static void convmatrix _PL((void));
static void conv3456 _PL((void));
void convperm _PL((void));
void conv1213 _PL((void));
int main  _PL((int argc, char *argv[]));


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long fl;
static int mod;
static long nor,noc;
static FILE *src;			/* Input file */
static FILE *out;			/* Output file */
static char lbuf[MAXLINE] = {0};	/* Input line buffer */
static char *lptr = lbuf;		/* Read pointer */
static char *inpname;
static char *outname = "P2";

static char *helptext[] = {
"SYNTAX",
"    zcv [<Inp> [<Out>]]",
"",
"FILES",
"    <Inp>    i  Text file. <Inp> = `=' means standard input.",
"                Default is `T1'.",
"    <Out>    o  Binary file. Dafault name is `P2'.",
NULL};


static proginfo_t pinfo =
   { "zcv", "Convert to Binary Format", "$Revision: 2.7 $", helptext };





/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(c)
int c;

{
    fprintf(stderr,"ZCV ERROR - ");
    switch (c)
    {
	case 'm':
	    fprintf(stderr,"UNKNOWN MODE\n");
	    break;
	case 'e':
	    fprintf(stderr,"UNEXPECTED EOF ON INPUT FILE\n");
	    break;
	case 'F':
	    errexit(ERR_FILEFMT,inpname);
	case 'r':
	    errexit(ERR_FILEREAD,inpname);
    }
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   readline() - Read next input line, skip empty lines and strip
   	comments. Returns 0 on success, 1 on end-of-file.
   ------------------------------------------------------------------ */

static int readline()

{	char *c;
	int mt = 1;

	while (mt)
	{	lbuf[0] = 0;
		if (feof(src)) return 1;
		fgets(lbuf,sizeof(lbuf),src);
		if (ferror(src)) err('r');
		for (c = lbuf; *c != 0 && *c != '#'; ++c)
			if (!isspace(*c)) mt = 0;
		*c = 0;
	}
	lptr = lbuf;
	return 0;
}


/* ------------------------------------------------------------------
   readlong() - Read integer number
   ------------------------------------------------------------------ */

static long readlong()

{	long l;
	int minus = 0;

	while (isspace(*lptr))
	{	while (isspace(*lptr)) ++lptr;
		if (*lptr == 0)
		{	if (readline())
				err('e');
		}
	}
	if (*lptr == '-') { minus = 1; ++lptr; }
	if (!isdigit(*lptr))
		err('F');
	for (l = 0; isdigit(*lptr); ++lptr)
	{	
		switch (*lptr)
		{   case '0': l = 10 * l + 0; break;
		    case '1': l = 10 * l + 1; break;
		    case '2': l = 10 * l + 2; break;
		    case '3': l = 10 * l + 3; break;
		    case '4': l = 10 * l + 4; break;
		    case '5': l = 10 * l + 5; break;
		    case '6': l = 10 * l + 6; break;
		    case '7': l = 10 * l + 7; break;
		    case '8': l = 10 * l + 8; break;
		    case '9': l = 10 * l + 9; break;
		}
	}
	return minus ? -l : l;
}



/* ------------------------------------------------------------------
   convmatrix() - Convert matrix in fixed format (mode 1)
   ------------------------------------------------------------------ */

static void convmatrix()

{	
	long i, j;
        PTR m1;
	long val = 0;
	int inp;

	if (fl > 9) err('e');
	zsetlen(fl,noc);
	m1 = zalloc((long)1);
	if ((out = zwritehdr(outname,fl,nor,noc)) == NULL)
	    errexit(-1,outname);
	for (i = 1; i <= nor; ++i)
	{	zmulrow(m1,F_ZERO);
		inp = 81;
		for (j = 1; j <= noc; ++j)
		{	if (inp >= 80)	/* read next line */
			{	memset(lbuf,0,sizeof(lbuf));
				if (readline()) err('e');
				inp = 0;
			}
			switch (lbuf[inp++])
			{	case ' ':
				case '0': val = 0; break;
				case '1': val = 1; break;
				case '2': val = 2; break;
				case '3': val = 3; break;
				case '4': val = 4; break;
				case '5': val = 5; break;
				case '6': val = 6; break;
				case '7': val = 7; break;
				case '8': val = 8; break;
				default: err('F');
			}
			if (val > fl) err('F');
			zinsert(m1,j,zitof(val));
		}
		zwritevec(out,m1,1);
	}
	fclose(out);
}


/* ------------------------------------------------------------------
   conv23456() - Convert matrix in free format (mode 3,4,5,6)
   ------------------------------------------------------------------ */

static void conv3456()

{	long i, j;
        PTR m1;
	long val;

	zsetlen(fl,noc);
	m1 = zalloc((long)1);
	if ((out = zwritehdr(outname,fl,nor,noc)) == NULL)
	    errexit(-1,outname);
	for (i = 1; i <= nor; ++i)
	{	zmulrow(m1,F_ZERO);
		for (j = 1; j <= noc; ++j)
		{	val = readlong();
			if (mod == 5)
			{	val %= zchar;
				if (val < 0) val += zchar;
			}
			zinsert(m1,j,zitof(val));
		}
		zwritevec(out,m1,1);
	}
	fclose(out);
}


/* ------------------------------------------------------------------
   convperm() - Convert permutation (mode 2)
   ------------------------------------------------------------------ */

void convperm()		/* mode 2 */

{
    long i, val;
    PTR m1;

    zsetlen(fl,noc);
    m1 = zalloc((long)1);
    if ((out = zwritehdr(outname,fl,nor,noc)) == NULL)
	errexit(-1,outname);
    for (i = 1; i <= nor; ++i)
    {	val = readlong();
	zmulrow(m1,F_ZERO);
	zinsert(m1,val,F_ONE);
	zwritevec(out,m1,1);
    }
    fclose(out);
}


/* ------------------------------------------------------------------
   conv1213() - Convert permutation (mode 12, 13)
   ------------------------------------------------------------------ */

void conv1213()		/* modes 12, 13 */

{
    long nper, i;
    long *buf;
    long kk = 0, f1, f2;


    buf = (long *) malloc(sizeof(long) * nor);
    if (buf == NULL) errexit(ERR_NOMEM,"zcv");
    if ((out = zwritehdr(outname,-fl,nor,noc)) == NULL)
	errexit(-1,outname);

    for (nper = noc; nper != 0; --nper)
    {	
	for (i = 0; i < nor; ++i)
	{
	    switch (mod)
	    {	case 12:
			kk = readlong();
			break;
		case 13:
			f1 = readlong();
			f2 = readlong();
			kk = (f1-1)*fl + f2 + 1;
			break;
	    }
	    buf[i] = kk;
	}
	if (zwritelong(out,buf,nor) != nor)
	{
	    perror(outname);
	    errexit(ERR_FILEWRITE,outname);
	}
    }
    fclose(out);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main (argc, argv)
int argc;
char *argv[];

{	
    char sfl[3], smode[7], snor[7], snoc[7];

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);

    switch (argc-opt_ind)
    {	case 2:
	    outname = argv[opt_ind + 1];
	case 1:
	    if (!strcmp(argv[opt_ind],"="))
	    {
		inpname = "[stdin]";
		src = stdin;
	    }
	    else
	    {
		inpname = argv[opt_ind];
		src = os_fopen(argv[opt_ind],FM_READ|FM_TEXT|FM_LIB);
		if (src == NULL)
		{
		    perror(inpname);
		    errexit(ERR_FILEOPEN,inpname);
		}
	    }
	    break;
	case 0:
	    inpname = "T1";
	    src = os_fopen("T1",FM_READ|FM_TEXT);
	    if (src == NULL)
	    {
		perror("T1");
		errexit(ERR_FILEOPEN,"T1");
	    }
	    break;
	default:
	    errexit(ERR_BADUSAGE,"zcv");
    }

    /* Read the header
       --------------- */
    if (readline()) err('e');
    strncpy(smode,lbuf,2); smode[2] = 0;
    strncpy(sfl,lbuf+2,6); sfl[6] = 0;
    strncpy(snor,lbuf+8,6); snor[6] = 0;
    strncpy(snoc,lbuf+14,6); snoc[6] = 0;
    if (sscanf(lbuf,"%d%ld%ld%ld",&mod,&fl,&nor,&noc) != 4)
    {
        sscanf(smode,"%d",&mod);
        sscanf(sfl,"%ld",&fl);
        sscanf(snor,"%ld",&nor);
        sscanf(snoc,"%ld",&noc);
    }

    if (mod != 1) readline();
    switch (mod)
    {
	case 1:
	    convmatrix();
	    break;
	case 2:
	    convperm();
	    break;
	case 3:
	case 4:
	case 5:
	case 6:
	    conv3456();
	    break;
	case 12:
	case 13:
	    conv1213();
	    break;
	default:
		err('m');
    }
    return (EXIT_OK);
}

