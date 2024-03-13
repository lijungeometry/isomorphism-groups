/* lib.d file storage.c */
/*This file is for all functions related to storage. There are also a
number of storage allocation routines which are \#defined in defs.webh.
|store_ptrs| is the number of pointers allocated by the system free
store allocation functions.
*/
#include <stdio.h>
#include "defs.h"
#include "list.h"
#include "word.h"
#include "input.h"

char * malloc_value = 0;
int store_ptrs=0;

void
Free_dp(p)
	dp p;
{
	assert(p);
	store_ptrs--;
	free((char *)p);
}

int
max(x,y)
	int x,y;
{
	return (x > y? x : y);
}

int
min(x,y) 
	int x,y;
{
	return (x < y? x : y);
}

