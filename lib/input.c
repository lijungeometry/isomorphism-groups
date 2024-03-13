/* lib.d rfile input.c */
#include <stdio.h>
#include <ctype.h>
#include "defs.h"
#include "input.h"
				
void
format_echocheck(number,rfile,wfile)
	char * number;
	FILE * rfile;
	FILE * wfile;
{
	char * string;
	string = vzalloc2(char,7);
	if (read_next_string(string,6,rfile)==FALSE){
		fprintf(stderr,"\t# Syntax error - no format version.\n");
		exit(2);
	}
	else {
		int i=0;
		int j;
		char * copy;
 
		fprintf(wfile,"Format %s\n",number);
		while (string[i]!='\0' && string[i]!=' ')
			i++;
/* so we see a null or blank character in position i */
		copy = vzalloc2(char,i+1);
		for (j=0;j<i;j++)
			copy[j] = string[j];
		copy[i] = '\0';
		if (strcmp(copy,number)!=0)
			fprintf(stderr,
		"\t# Warning,Expected Format %s, received Format %s\n", number,copy);
		Free_dp((dp)copy); copy=0;
	}
	Free_dp((dp)string); string=0;
}
 
void
format_check(number,rfile)
	char * number;
	FILE * rfile;
{
	char * string;
	string = vzalloc2(char,7);
	if (read_next_string(string,6,rfile)==FALSE){
		fprintf(stderr,"\t# Syntax error - no format version.\n");
		exit(2);
	}
	else {
		int i=0;
		int j;
		char * copy;
 
		while (string[i]!='\0' && string[i]!=' ')
			i++;
/* so we see a null or blank character in position i */
		copy = vzalloc2(char,i+1);
		for (j=0;j<i;j++)
			copy[j] = string[j];
		copy[i] = '\0';
		if (strcmp(copy,number)!=0)
			fprintf(stderr,
		"\t# Warning. Expected Format %s, received Format %s\n", number,copy);
		Free_dp((dp)copy); copy=0;
	}
	Free_dp((dp)string); string=0;
}
		
void
bad_data()
{
	fprintf(stderr,"Bad data.\n");
	exit(2);
}
 
 
boolean
read_next_float(fp,rfile)
	float *fp;
	FILE * rfile;
{
	boolean ans=TRUE;
	int sign=1;
	int c;
	boolean fractional=FALSE;
	float multiplier=0.1;
	*fp=0;
	while ((c=read_char(rfile))!=EOF && !(isdigit(c))){
		if (c=='-')
			sign = -1;
		else if (c=='}'){
			ungetc('}',rfile);
	/* a '}' is used to terminate the input */
			ans = FALSE;
			break;
		}
		else
			sign = 1;
	}
	if (c==EOF)
		bad_data();
	if (ans==TRUE){
		do {
			if (c == '.' && fractional==FALSE)
				fractional = TRUE;
			else {
				if (fractional==FALSE)
					*fp = c-'0' + 10*(*fp);
				else {
					*fp += multiplier*(c-'0');
					multiplier *= 0.1;
				}
			}
			c=getc(rfile);
		} while ((isdigit(c))||(c=='.' && fractional==FALSE));
		ungetc(c,rfile);
	}
	*fp = sign * (*fp);
	return ans;
}
	
boolean read_next_int(kp,rfile)
	int * kp;
	FILE * rfile;
{
	boolean ans=TRUE;
	int sign=1;
	int c;
	*kp=0;
	while ((c=read_char(rfile))!=EOF && !(isdigit(c))){
		if (c=='-')
			sign = -1;
		else if (c=='}' || c==';') {
	/* a '}' or ';' is used to terminate the input */
			ungetc(c,rfile);
			ans = FALSE;
			break;
		}
		else
			sign = 1;
	}
	if (c==EOF)
		bad_data();
	if (ans==TRUE) {
		do {
			*kp = 10*(*kp) + c - '0';
		} while (isdigit(c=getc(rfile)));
		ungetc(c,rfile);
	}
	*kp = sign*(*kp);
	return ans;
}
	
boolean read_next_letter(lp,rfile)
	int * lp;
	FILE * rfile;
{
	boolean ans=TRUE;
	while ((*lp=read_char(rfile))!=EOF && !(isalpha(*lp)))
		if (*lp=='}') {
	/* a '}' is used to terminate the input */
			ungetc(*lp,rfile);
			ans = FALSE;
			break;
		}
	if (*lp==EOF)
		bad_data();
	return ans;
}
 
/* The next complete string is read and the first n or less  non-null characters are
stored in cp. Any space between the last non-null character of the string
and
 the terminating null character of cp is filled
with blank spaces.  Spaces, tabs and
returns, and also matched pairs {} (and all the stuff between them) are skipped
over until other symbols are met. The string is terminated by a tab, space,
return , (,), { or } (which remains unread ).
A string beginning with } causes the function to stop
reading and return false.  The } is then returned to the buffer as if
unread.
*/
boolean
read_next_string(cp,n,rfile)
	char * cp;
	int n;
	FILE * rfile;
{
	int i=0;
	int c;
	boolean ans=TRUE;
	while ((c=read_char(rfile))!=EOF){
		if (i==0){
			if ( c=='}' ) {
				ungetc(c,rfile);
				return  FALSE;
			}
			else if (c=='{'){ /* skip over {}'s */
				int count=1;
				while (count>0){
					if  ((c=read_char(rfile))=='{')
						count++;
					else if (c=='}')
						count--;
					if (c==EOF) break;
				}
			}
			else if (c==' ')
				continue;
			else{
				cp[0]=c;
				i=1;
			}
		}
		else if (i>0){
	if (c==' '||c=='{'||c=='}'||c=='('||c==')') {
					if (c!=' ') ungetc(c,rfile);
					break;
			}
			else {
				if (i<n)
					cp[i++]=c;
			}
		}
	}
	if (c==EOF){
		ans = FALSE;
		i=0; /* if the function's returning false, we'd like to empty the
string */
	}
	while (i<n)
		cp[i++]=' ';
	cp[n]='\0';
	return ans;
}
 
boolean
find_keyword(label,rfile)
	char * label;
	FILE * rfile;
{
	char * string;
	boolean found=TRUE;
	int k=0;
	while (label[k]!='\0')
		k++;
	string = vzalloc2(char,k+1);
	do {
		found=TRUE;
		if (read_next_string(string,k,rfile)==FALSE){
			found = FALSE;
			break;
		}
	} while ((strcmp(string,label)!=0));
	Free_dp((dp)string); string=0;
	return found;
}
 
 
/* The next function should be used in place of getc. It reads and returns
the next character, except that comments (preceded by # and going to the
next new line) and tabs and new lines are replaced by a single space.
*/
int
read_char(rfile)
	FILE * rfile;
{ int n;
  n=getc(rfile);
  /* The input may contain a '\' character which is printed before '\n'
    -- this is particularly true if the input file is generated by GAP */
  if (n=='\\') {
     n = getc(rfile);
     if (n == '\n') n = getc (rfile); 
  }
  if (n=='#') while ((n=getc(rfile))!='\n' && n!=EOF);
  if (n=='\n' || n=='\t') return(' ');
  return(n);
}
 
void
find_char(c,rfile)
	char c;
	FILE * rfile;
{ int n,m;
  char d;
	d=c;
        m= (int)d;
	while ((n=read_char(rfile))!=m){
		if (n==EOF){
			fprintf(stderr,"Unexpected end of file.\n");
			bad_data();
		}
	}
	return;
}
 
boolean
read_next_char(cp,rfile)
  int * cp;
  FILE * rfile;
{
  while ((*cp=read_char(rfile))==' ');
  if (*cp==EOF) bad_data();
  else if (*cp=='}'){ ungetc(*cp,rfile); return FALSE;}
  else return TRUE;
}

