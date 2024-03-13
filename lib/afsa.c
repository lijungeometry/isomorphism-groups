/* lib.d file afsa.c */
#include <stdio.h>
#include "defs.h"
#include "list.h"
#include "word.h"
#include "afsa.h"
extern int gen_array_size;

void
afsa_init(fsap)
afsa * fsap;
{
	fsap->states = 0;
	fsap->symbols = 0;
	fsap->base_symbols = 0;
	fsap->bfs = FALSE;
	fsap->min = FALSE;
	fsap->variables = 1;
	fsap->eos = TRUE;
	fsap->array = 0;
}

void
afsa_clear(fsap)
	afsa* fsap;
{
	int i;
	int ** array=fsap->array;
	for (i=0;i<=(fsap->symbols);i++){
		Free_dp((dp)array[i]);
		array[i]=0;
	}
	Free_dp((dp)array);
	array=0;
}

afsa *
afsa_read(base_alphabetp,file)
	word ** base_alphabetp;
	FILE * file;
{
	afsa * fsap;
	word * base_alphabet=0;
	int * perm=0;
	int symbols=0;
	int base_symbols=0;	
	int ** array=0;
	int i,j;
	int c;
	char * label;
	int letter;
	boolean compressed=FALSE;
	fsap=vzalloc1(afsa);
	afsa_init(fsap);
	label=vzalloc2(char,9);
	find_char('{',file);
	read_next_string(label,8,file);
	if (strcmp(label,"states  ")!=0)
		bad_data();
	else
		read_next_int(&(fsap->states),file);
	read_next_string(label,8,file);
	if (strcmp(label,"symbols ")!=0)
		bad_data();
	else
		read_next_int(&symbols,file);
	if (symbols+1>gen_array_size)
		gen_array_size=symbols+1;
	while (read_next_string(label,8,file)&&strcmp(label,"%       ")!=0){
		if (strcmp(label,"bfs     ")==0)
			fsap->bfs = TRUE;
		else if (strcmp(label,"min     ")==0)
			fsap->min = TRUE;
		else if (strcmp(label,"variable")==0)
			read_next_int(&(fsap->variables),file);
		else if (strcmp(label,"no-eos  ")==0 || strcmp(label,"no_eos  ")==0)
			fsap->eos=FALSE;
		else if (strcmp(label,"base-alp")==0 || strcmp(label,"base_alp")==0){
			find_char('{',file);
			if (base_alphabet==0){
				base_alphabet=vzalloc2(word,gen_array_size);
				for (i=0;i<gen_array_size;i++)
					word_init(base_alphabet+i);
				while (read_next_int(&i,file)){
					read_next_word(base_alphabet+i,file);
					base_symbols++;
				}
			}
			find_char('}',file);
		}	
		else if (strcmp(label,"alphabet")==0){
			if (base_alphabet==0){
				base_alphabet=vzalloc2(word,gen_array_size);
				for (i=0;i<gen_array_size;i++)
					word_init(base_alphabet+i);
				find_char('{',file);
				while (read_next_int(&i,file))
					if (fsap->variables==1){
						read_next_word(base_alphabet+i,file);
						base_symbols++;
					}
				find_char('}',file);
			}
		}
	}
	perm=vzalloc2(int,base_symbols+1);
	if (*base_alphabetp == 0){
		*base_alphabetp = base_alphabet;
		for (i=1;i<=base_symbols;i++)
			perm[i]=i;
	}
	else {
/* now we should compare the base alphabet we've just read with the one we
were expecting. We assume that the symbols used are the same (plus or minus
the $ symbol), but the order may be different. In that case perm tells us how
to permute them. */
		for (i=1;i<=base_symbols;i++){
			for (j=1;j<=base_symbols;j++){
				if ((*base_alphabetp)+j==0 || 
						word_sgn(base_alphabet+i,(*base_alphabetp)+j)==0){
					perm[i]=j;
					break;
				}
			}
		}
		for (i=0;i<gen_array_size;i++)
			word_clear(base_alphabet+i);
		Free_dp((dp)base_alphabet); base_alphabet=0;
	}
	if (fsap->variables>1){
		int * permcopy=perm;
		perm=vzalloc2(int,symbols+2); /* we put in some extra space to avoid
problems when there's no eos, and so symbols is not base_symbols squared */
		for (i=1;i<=base_symbols;i++)
			for (j=1;j<=base_symbols;j++)
				perm[(i-1)*base_symbols + j] = (permcopy[i]-1)*base_symbols
														+ permcopy[j];
		Free_dp((dp)permcopy); permcopy=0;
	}
	array=vzalloc2(int *,symbols+1);
	for (i=0;i<=symbols;i++)
		array[i]=vzalloc2(int ,(fsap->states)+1);
	read_next_string(label,8,file);
	if (strcmp(label,"ctable  ")==0)
		compressed=TRUE;
	j=1;
	while (j<=(fsap->states)){
		(void)read_next_letter(&letter,file);
		if (letter=='N'){
			if ((i=read_char(file))!='g'){
				if (i!=' ') ungetc(i,file);
				array[0][j]=0;
			}
			else {
				read_next_int(&i,file);
				array[0][j]=2+i;
			}
		}
		else if (letter=='A')
			array[0][j]=1;
		else
			fprintf(stderr,"\t# Invalid category\n");
		if (compressed){
			while (read_next_int(&i,file))
				(void)read_next_int((array[perm[i]])+j,file);
			find_char(';',file);
		}
		else 	
			for (i=1;i<=base_symbols;++i) 
				(void)read_next_int((array[perm[i]])+j,file);
		j++;
	} 
	find_char('}',file);
	Free_dp((dp)perm); perm=0;
	Free_dp((dp)label); label=0;
	fsap->symbols=symbols;
	fsap->base_symbols=base_symbols;
	fsap->array=array;
	return fsap;
}
	
void
afsa_print(wfile,fsap)
	FILE * wfile;
	afsa * fsap;
{
	int i,j,k;
	int ** array=fsap->array;
	boolean compressed=FALSE,temp;
	fprintf(wfile,"fsa {\n");
	fprintf(wfile,"\tstates %d\n",fsap->states);
	fprintf(wfile,"\tsymbols %d\n",fsap->symbols);
	if ((fsap->min) == TRUE)
		fprintf(wfile,"\tmin\n");
	if ((fsap->bfs) == TRUE)
		fprintf(wfile,"\tbfs\n");
	fprintf(wfile,"\tvariables %d\n",fsap->variables);
    temp= fsap->eos;
	if (temp == FALSE)
		fprintf(wfile,"\tno_eos\n");
	if ((fsap->variables)==1)
		fprintf(wfile,"\talphabet {");
	else
		fprintf(wfile,"\tbase_alphabet {");
	for (i=1;i<=fsap->base_symbols;i++){
		fprintf(wfile,"%d = ",i);
		gen_print(wfile,i);
		fprintf(wfile," ");
		if (i%10==0)
			fprintf(wfile,"\n");
	}
	fprintf(wfile," }\n");
	
	if ((fsap->variables)==2){
		fprintf(wfile,"\talphabet {");
		k=1;
		for (i=1;i<=fsap->base_symbols;i++){
			for (j=1;j<=fsap->base_symbols;j++){
				if (fsap->eos==FALSE && i==fsap->base_symbols &&
												j==fsap->base_symbols)
					break;
				fprintf(wfile,"%d = ",k);
				fprintf(wfile,"\(");
				gen_print(wfile,i);
				fprintf(wfile,",");
				gen_print(wfile,j);
				fprintf(wfile,") ");
				if (k%5==0)
					fprintf(wfile,"\n");
				k++;
			}
		}
		fprintf(wfile," }\n");
	}
	fprintf(wfile,"\tstart { 1 }\n");
	if (fsap->variables==1)
		fprintf(wfile,"\n%%\natable\n");
	else{
		compressed=TRUE;
		fprintf(wfile,"\n%%\nctable\n");
	}
	for (i=1;i<=fsap->states;i++){
		int entries=0;
		if (array[0][i]==0)
			fprintf(wfile," %3d  N\t",i);
		else if (array[0][i]==1)
			fprintf(wfile," %3d  A\t",i);
		else 
			fprintf(wfile," %3d  Ng%d\t",i,array[0][i]-2);
		for (j=1;j<=fsap->symbols;j++){
			k=array[j][i];
			if (compressed){
				if (k!=0){
					if (entries!=0){
						fprintf(wfile,",");
						if (entries%8==0)
							fprintf(wfile,"\n         \t");
					}
					fprintf(wfile," %d > %d",j,k);
					entries++;
				}
			}
			else{
				if (j!=1 && j%15==1)
					fprintf(wfile,"\n\t      ");
				fprintf(wfile," %3d",k);
			}
		}
		fprintf(wfile,";\n");
	}
	fprintf(wfile,"}\n");
}

void
twoafsa_init(fsap)
twoafsa * fsap;
{
	fsap->states = 0;
	fsap->symbols=0;
	fsap->base_symbols = 0;
	fsap->bfs = FALSE;
	fsap->min = FALSE;
	fsap->eos = TRUE;
	fsap->array = 0;
}

void
twoafsa_clear(fsap)
	twoafsa* fsap;
{
	int i,j;
	int *** array=fsap->array;
	for (i=0;i<=(fsap->base_symbols);i++){
		for (j=0;j<=(fsap->base_symbols);j++){
			Free_dp((dp)array[i][j]);
			array[i][j]=0;
		}
		Free_dp((dp)array[i]);
		array[i]=0;
	}
	Free_dp((dp)array);
	array=0;
}

twoafsa *
twoafsa_read(base_alphabetp,file)
	word ** base_alphabetp;
	FILE * file;
{
	twoafsa * fsap;
	word * base_alphabet=0;
	int * perm=0;
	int symbols=0;
	int base_symbols=0;	
	int *** array=0;
	int i,j,k;
	int c;
	char * label;
	int letter;
	boolean compressed=FALSE;
	fsap=vzalloc1(twoafsa);
	twoafsa_init(fsap);
	label=vzalloc2(char,9);
	find_char('{',file);
	read_next_string(label,8,file);
	if (strcmp(label,"states  ")!=0)
		bad_data();
	else
		read_next_int(&(fsap->states),file);
	read_next_string(label,8,file);
	if (strcmp(label,"symbols ")!=0)
		bad_data();
	else
		read_next_int(&symbols,file);
	if (symbols+1>gen_array_size)
		gen_array_size=symbols+1;
	base_alphabet=vzalloc2(word,gen_array_size);
	for (i=0;i<gen_array_size;i++)
		word_init(base_alphabet+i);
	while (read_next_string(label,8,file)&&strcmp(label,"%       ")!=0){
		if (strcmp(label,"bfs     ")==0)
			fsap->bfs = TRUE;
		if (strcmp(label,"min     ")==0)
			fsap->min = TRUE;
		else if (strcmp(label,"no-eos  ")==0 ||strcmp(label,"no_eos  ")==0)
			fsap->eos=FALSE;
		else if (strcmp(label,"base-alp")==0 ||strcmp(label,"base_alp")==0){
			find_char('{',file);
			while (read_next_int(&i,file)){
				read_next_word(base_alphabet+i,file);
				base_symbols++;
			}
			find_char('}',file);
		}	
	}
	perm=vzalloc2(int,base_symbols+1);
	if (*base_alphabetp == 0){
		*base_alphabetp = base_alphabet;
		for (i=1;i<=base_symbols;i++)
			perm[i]=i;
	}
	else {
/* now we should compare the base alphabet we've just read with the one we
were expecting. We assume that the symbols used are the same (plus or minus
the $ symbol), but the order may be different. In that case perm tells us how
to permute them. */
		for (i=1;i<=base_symbols;i++){
			for (j=1;j<=base_symbols;j++){
				if ((*base_alphabetp)+j==0 || 
						word_sgn(base_alphabet+i,(*base_alphabetp)+j)==0){
					perm[i]=j;
					break;
				}
			}
		}
		for (i=0;i<gen_array_size;i++)
			word_clear(base_alphabet+i);
		Free_dp((dp)base_alphabet); base_alphabet=0;
	}
	array=vzalloc2(int **,symbols+1);
	for (i=0;i<=base_symbols;i++){
		array[i]=vzalloc2(int* ,base_symbols+1);
		for (j=0;j<=base_symbols;j++)
			array[i][j]=vzalloc2(int,(fsap->states)+1);
	}
	read_next_string(label,8,file);
	if (strcmp(label,"ctable  ")==0)
		compressed=TRUE;
	k=1;
	while (k<=(fsap->states)){
		(void)read_next_letter(&letter,file);
		if (letter=='N'){
			if ((i=read_char(file))!='g'){
				if (i!=' ') ungetc(i,file);
				array[0][0][k]=0;
			}
			else {
				read_next_int(&i,file);
				array[0][0][k]=2+i;
			}
		}
		else if (letter=='A')
			array[0][0][k]=1;
		else
			fprintf(stderr,"\t# Invalid category\n");
		if (compressed){
			while (read_next_int(&i,file)){
				(void)read_next_int
	((array[perm[(int)((i-1)/base_symbols) +1]] [perm[(i-1)%base_symbols+1]])+k,
																	file);
			}
			find_char(';',file);
		}
		else {	
			for (i=1;i<=base_symbols;++i) 
				for (j=1;j<=base_symbols;j++){
					if (fsap->eos==FALSE && perm[i]==base_symbols 
							&& perm[j]==base_symbols)
						break;
					(void)read_next_int((array[perm[i]][perm[j]])+k,file);
				}
		}
		k++;
	} 
	find_char('}',file);
	Free_dp((dp)perm); perm=0;
	Free_dp((dp)label); label=0;
	fsap->symbols=symbols;
	fsap->base_symbols=base_symbols;
	fsap->array=array;
	return fsap;
}
	
void
twoafsa_print(wfile,fsap)
	FILE * wfile;
	twoafsa * fsap;
{
	int i,j,k,l;
	int *** array=fsap->array;
	boolean  temp;
	fprintf(wfile,"fsa {\n");
	fprintf(wfile,"\tstates %d\n",fsap->states);
	fprintf(wfile,"\tsymbols %d\n",fsap->symbols);
	if (fsap->min)
		fprintf(wfile,"\tmin\n");
	if (fsap->bfs)
		fprintf(wfile,"\tbfs\n");
	fprintf(wfile,"\tvariables 2\n");
	temp=fsap->eos;
	if (temp == FALSE)	
		fprintf(wfile,"\tno_eos\n");
	fprintf(wfile,"\tbase_alphabet {");
	for (i=1;i<=fsap->base_symbols;i++){
		fprintf(wfile,"%d = ",i);
		gen_print(wfile,i);
		fprintf(wfile," ");
		if (i%10==0)
			fprintf(wfile,"\n");
	}
	fprintf(wfile," }\n");
	
	fprintf(wfile,"\talphabet {");
	k=1;
	for (i=1;i<=fsap->base_symbols;i++){
		for (j=1;j<=fsap->base_symbols;j++){
			if (fsap->eos==FALSE && i==fsap->base_symbols
				 	&& j==fsap->base_symbols)
				break;
			fprintf(wfile,"%d = ",k);
			fprintf(wfile,"\(");
			gen_print(wfile,i);
			fprintf(wfile,",");
			gen_print(wfile,j);
			fprintf(wfile,") ");
			if (k%5==0)
				fprintf(wfile,"\n");
			k++;
		}
	}
	fprintf(wfile," }\n");
	fprintf(wfile,"\tstart { 1 }\n");
	fprintf(wfile,"\n%%\nctable\n");
	for (i=1;i<=fsap->states;i++){
		int entries=0;
		if (array[0][0][i]==0)
			fprintf(wfile," %3d  N\t",i);
		else if (array[0][0][i]==1)
			fprintf(wfile," %3d  A\t",i);
		else 
			fprintf(wfile," %3d  Ng%d\t",i,array[0][0][i]-2);
		for (j=1;j<=fsap->base_symbols;j++){
			for (k=1;k<=fsap->base_symbols;k++){
				if (fsap->eos==FALSE && j==fsap->base_symbols
					 	&& k==fsap->base_symbols)
					break;
				l=array[j][k][i];
				if (l!=0){
					if (entries!=0){
						fprintf(wfile,",");
						if (entries%8==0)
							fprintf(wfile,"\n         \t");
					}
					fprintf(wfile,"%d > %d",(fsap->base_symbols)*(j-1)+k,l);
					entries++;
				}
			}
		}
		fprintf(wfile,";\n");
	}
	fprintf(wfile,"}\n");
}

afsa * afsa_eosdelete(fsap)
	afsa * fsap;
{
	afsa * newfsap;
	int ** array=fsap->array;
	int ** new_array;
	int i,k,l,m;
	int eosstate=0;
	for (k=1;k<=fsap->states;k++)
		if ((eosstate=array[fsap->symbols][k])>0)
			break;
	newfsap = vzalloc1(afsa);
	afsa_init(newfsap);
	newfsap->bfs = fsap->bfs;
	newfsap->min = fsap->min;
	newfsap->variables = fsap->variables;
	newfsap->eos = FALSE;
	if (eosstate!=0)
		newfsap->states = fsap->states - 1;
	else
		newfsap->states = fsap->states;
	if (fsap->variables == 1)
		newfsap->base_symbols = fsap->base_symbols -1;
	else
		newfsap->base_symbols = fsap->base_symbols;
	newfsap->symbols = fsap->symbols - 1;
	new_array = vzalloc2(int *,newfsap->symbols + 1);
	for (i=0;i<=newfsap->symbols;i++)
		new_array[i]=vzalloc2(int,newfsap->states + 1);
	if (eosstate!=0){
		l=1;
		for (k=1;k<=fsap->states;k++){
			if (k!=eosstate){
				if (array[fsap->symbols][k]>0)
					new_array[0][l]=ACCEPTSTATE;
				for (i=1;i<=newfsap->symbols;i++){
					if ((m=array[i][k])<eosstate)
						new_array[i][l]=m;
					else
						new_array[i][l]=m-1;
				}
			l++;
			}
		}
	}
	else {
/* its unlikely that eosstate is zero, but possible */
		for (k=1;k<=fsap->states;k++)
			for (i=1;i<=newfsap->symbols;i++)
				new_array[i][k] = array[i][k];
	}
	newfsap->array=new_array;
	return newfsap;
}

afsa * afsa_eosadd(fsap)
	afsa * fsap;
{
	afsa * newfsap;
	int ** array=fsap->array;
	int ** new_array;
	int i,j,k,l,top;
	int eosstate=0;
	assert(fsap->eos==FALSE);
	newfsap = vzalloc1(afsa);
	afsa_init(newfsap);
	newfsap->bfs = fsap->bfs;
	newfsap->min = fsap->min;
	newfsap->variables = fsap->variables;
	newfsap->eos = TRUE;
	newfsap->states = fsap->states + 1;
	if (fsap->variables == 1)
		newfsap->base_symbols = fsap->base_symbols +1;
	else
		newfsap->base_symbols = fsap->base_symbols;
	newfsap->symbols = fsap->symbols + 1;
	new_array = vzalloc2(int *,newfsap->symbols + 1);
	for (i=0;i<=newfsap->symbols;i++)
		new_array[i]=vzalloc2(int,newfsap->states + 1);
	l=1;
	top=1;
	for (k=1;k<=fsap->states;k++){
		for (i=1;i<=fsap->symbols;i++){
			if ((j=array[i][k])>top){
				if (eosstate==0)
					top=j;
				else  
	/* eosstate=top+1, all the states above top need to be renumbered */
					j++;
			}
			new_array[i][l]=j;
		}
		if (array[0][k]==ACCEPTSTATE){
			new_array[0][l]=NONACCEPTSTATE;
			if (eosstate==0){
				eosstate=top+1;
				new_array[0][eosstate]=ACCEPTSTATE;
				new_array[newfsap->symbols][eosstate]=eosstate;
			}
			new_array[newfsap->symbols][l]=eosstate;
		}
		l++;
		if (l==eosstate)
			l++;
	}
	if (eosstate==0) /* unlikely but possible if the fsa has an empty
language */
		(newfsap->states)--;
	newfsap->array=new_array;
	return newfsap;
}

int
afsa_num_negstates(fsap)
	afsa * fsap;
{
	int k=0;
	int i,j;
	int ** array=fsap->array;
	for (i=1;i<=fsap->symbols;i++)
		for (j=1;j<=fsap->states;j++)
			if (array[i][j]<k)
				k=array[i][j];
	return (-k);
}

