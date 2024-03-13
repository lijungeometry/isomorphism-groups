/* lib file list.c */
#include <stdio.h> 
#include <ctype.h> 
#include "defs.h"
#include "list.h"

#define SIGNATURE 0
#define COPY 1
#define CREATE 2
#define KILL 3
#define PRINT 4
static boolean list_ordered_insert();
static boolean list_fifo_insert();
static boolean list_lifo_insert();

void 
list_init(lp,fntab,type) 
	list *lp; 
	elt_fntab *fntab; 
        char type;
{ 
	lp->first.next = NULL; 
	lp->first.dp = NULL; 
	lp->lastp = &(lp->first); /* the first link is also the last */
	lp->fntab = fntab; 
        lp->type = type;
} 

void 
list_reset(lp) 
  list * lp;
{ 
	link * p; 
	link * q; 
	auto void (*killfnp) PARMS ((dp)); 

	p = lp->first.next;/* |p| is the pointer from the first link to the
						next (i.e. the second) one. */ 
	while(p) { /* while there is still a second link */ 
	
		q = p; /* |q| is the second link */
		killfnp = lp->fntab->kill; 
		p= p->next; /* move |p| to point to the third link */
		(*killfnp)(q->dp);/* free the space occupied by the data  */  
		q->dp=0;
		Free_dp((dp)q); q = 0; /*free the space occupied by the link*/
  /* (p could become a null pointer) */ 
	}
	lp->first.next = NULL ;/* first is the only link now */ 
	lp->lastp = &(lp->first); /*the first link is also the last link*/
} 


void
list_clear(lp)
list * lp;
{
  list_reset(lp);
}
/* Initialize a |list_traverser| to make it traverse a particular list, 
whose address is given in the second argument.  A |list_traverser| attached 
to a list consists of a pointer (|upto|) to a link in
that list together with the function table of the list . At initialization
the pointer is set to point to the first link in the list .
The list traverser doesn't need to be assigned space, since it is always set
equal to something which already exists.
*/
void list_traverser_init(ltp,lp) 
	list_traverser* ltp; 
	list* lp; 
{ 
	ltp->upto = &(lp->first); 
	ltp->fntab = lp->fntab; 
} 

void
list_traverser_clear(ltp)
  list_traverser * ltp;
{
}

/* |list_next()| is used when traversing a list . It
returns |TRUE| if there is a next element on the list (i.e. if the
|list_traverser| points to a link which is not the last link),
and |FALSE| otherwise. If there is a next item, the |list_traverser| is
advanced one item and the next item is put into the
address given by the datapointer |out| (using the copy function in the
relevant function table).  The copy function clears |*out| 
and assigns the more space to it as necessary, so there is no need to
reset it if it has been used before. |*out| should however be
initialized. 

The procedure does not change the list being traversed.
*/
boolean 
list_next(ltp,out) 
	list_traverser * ltp; 
	dp  out; 
{ 

	if (ltp->upto->next ==  NULL) /*the traverser is set to point at 
the last link */ 
		return FALSE; 
	ltp->upto = ltp->upto->next; /*otherwise move the traverser to point
to the next link */ 
	(*(ltp->fntab->copy))(ltp->upto->dp,out); 
/* copy the data in the next link into the address given by |out| */
	return TRUE; 
} 

/* The Boolean valued function |list_empty()| returns |TRUE| if the list
 is empty, |FALSE| otherwise without changing the list in any way.
*/
boolean list_empty(lp)
  list * lp;
{ 
	if (lp->first.next == NULL) 
		return TRUE; 
	else 
		return FALSE; 
} 

boolean list_insert(lp,in)
  list * lp;
  dp in;
{
  switch(lp->type){
    case ORDERED:
      return list_ordered_insert(lp,in);
      break;
    case FIFO:
      return list_fifo_insert(lp,in);
      break;
    case LIFO:
      return list_lifo_insert(lp,in);
      break;
  }
}
    
static boolean
list_ordered_insert(lp,in) 
        list * lp;
	dp in; 
{ 
	link *p; 
	auto int (*sgnfnp) PARMS ((dp,dp)); 
	auto void (*copyfnp) PARMS ((dp,dp)); 
	boolean ans = TRUE;
	assert(in); 
	sgnfnp = lp->fntab->signature; 
	copyfnp = lp->fntab->copy; 
	p = &(lp->first); 
	while ((p->next) &&  (*sgnfnp)(in,p->next->dp) == -1) 
				p = p->next; 
	if ((p->next)&& (*sgnfnp)(in,p->next->dp) == 0) { 
			/* something similar enough to |*in| has been found */
		(*copyfnp)(in,p->next->dp);
		ans = FALSE;
	} 
	else { 
		link * linkp; 
		auto data_ptr (*createfnp) PARMS ((VOID)); 
	
		createfnp = lp->fntab->create; 
		copyfnp = lp->fntab->copy; 
		linkp = valloc1(link); 
		linkp->next = p->next; /* what was the "next link" becomes the
				link after the new link. If we're at the end of
				the list this is a null pointer.*/ 
		linkp->dp = (*createfnp)(); 
		(*copyfnp)(in,linkp->dp);  
		if (p->next == 0) 
			lp->lastp = linkp; /* if we were at the end of the
	list, the new link is now the last link */ 
	p->next = linkp;/*the new link is now the "next link" after |*p|*/
	} 
	return ans;
} 


static boolean
list_fifo_insert(lp,in) 
	list * lp; 
	dp in; 
{ 
	auto void (*copyfnp) PARMS ((dp,dp)); 
	link *p; 
	link * linkp; 
	auto data_ptr (*createfnp) PARMS ((VOID)); 
	assert(in); 
	p = lp->lastp; 

	createfnp = lp->fntab->create; 
	copyfnp = lp->fntab->copy; 
	linkp = valloc1(link); 
	linkp->next = p->next; /* what was the "next link" becomes the
			link after the new link. If we're at the end of
			the list this is a null pointer.*/ 
	linkp->dp = (*createfnp)(); 
	(*copyfnp)(in,linkp->dp);  
	if (p->next == 0) 
		lp->lastp = linkp; /* if we were at the end of the
list, the new link is now the last link */ 
	p->next = linkp;/*the new link is now the "next link" after |*p|*/
	return TRUE;
} 

static boolean
list_lifo_insert(lp,in) 
	list * lp; 
        dp in;
{ 
	auto void (*copyfnp) PARMS ((dp,dp)); 
	link *p; 
	link * linkp; 
	auto data_ptr (*createfnp) PARMS ((VOID)); 
	assert(in); 
	p = &(lp->first); 

	createfnp = lp->fntab->create; 
	copyfnp = lp->fntab->copy; 
	linkp = valloc1(link); 
	linkp->next = p->next; /* what was the "next link" becomes the
			link after the new link. If we're at the end of
			the list this is a null pointer.*/ 
	linkp->dp = (*createfnp)(); 
	(*copyfnp)(in,linkp->dp);  
	if (p->next == 0) 
		lp->lastp = linkp; /* if we were at the end of the
list, the new link is now the last link */ 
	p->next = linkp;/*the new link is now the "next link" after |*p|*/
	return TRUE;
} 

/*
Only used for ordered lists */
boolean list_delete(lp,in) 
	list * lp;
        dp in; 
{ 
	link *p; 
	auto int (*sgnfnp) PARMS ((dp,dp)); 

	assert(in); 
        assert(lp->type==ORDERED);
	sgnfnp = lp->fntab->signature;
	p = &(lp->first); 
	while ((p->next) &&  (*sgnfnp)(in,p->next->dp) == -1) 
				p = p->next; 
	if ((p->next)&&(*sgnfnp)(in,p->next->dp) == 0) { 
		/* something similar enough to |*in| has been found */
		link * temp; 
		auto void (*killfnp) PARMS ((dp)); 
		killfnp = lp->fntab->kill; 
		temp = p->next; 
		(*killfnp)(temp->dp);  
		temp->dp=0;
		if (temp->next == 0) /*we're deleting the last link*/ 
			lp->lastp = p; 
		p->next = temp->next;/* what was the link after next becomes
			the next link */ 
		Free_dp((dp)temp); temp = 0; 
		return TRUE; 
	} 
	else return FALSE; 
} 


boolean 
list_get_first(lp,out) 
	list * lp;
        dp  out; 
{ 
	if (lp->first.next != NULL){ /* there is a second link */
		auto void (*copyfnp) PARMS ((dp,dp)); 
		copyfnp = lp->fntab->copy; 
		(*copyfnp)(lp->first.next->dp,out); 
		return TRUE; 
	} 
	else return FALSE; 
} 


boolean 
list_delget_first(lp,out) 
	list * lp; 
	dp  out; 
{ 
	link * linkp; 
	if ((linkp = lp->first.next) != NULL){ /* there is a second link */
		auto void (*copyfnp) PARMS ((dp,dp)); 
		auto void (*killfnp) PARMS ((dp)); 
		copyfnp = lp->fntab->copy; 
		(*copyfnp)(linkp->dp,out); 
		killfnp = lp->fntab->kill; 
		(*killfnp)(linkp->dp); 
		linkp->dp=0;
		lp->first.next = linkp->next; 
		Free_dp((dp)linkp); linkp = 0; 
		if (lp->first.next == NULL) /*empty dp*/ 
			lp->lastp = &(lp->first); 
		return TRUE; 
	} 
	else return FALSE; 
 
} 




/* Print for ordinary use.
The data in each link after the first  (where there is no data stored)
is printed out using the 
print function from the relevant function table. If a link contains no
data the message "NULL data" is printed.
List items are now separated by commas. Set the flag `comma' TRUE if you
want a comma after the last item.
*/
void 
list_print(wfile,lp,comma) 
	FILE * wfile;
	list * lp; 
        boolean comma;
{ 
	link * p; 
	auto void (*printfnp) PARMS ((FILE *,dp)); 
 
	assert(lp); 
	printfnp = lp->fntab->print; 
	p = (lp->first).next; 
	while(p) { 
		if (p->dp) {
			(*printfnp)(wfile,p->dp); 
                        if (p->next || comma) fprintf(wfile,",\n"); 
			else fprintf(wfile,"\n"); 
		}
		else
			fprintf(wfile," NULL data "); 
		p = p->next; 
	} 
} 


/*
The original contents of |*newLLp| are destroyed.
The two lists are assumed to be of the same type, i.e. they contain the
same kind of data and use the same function table.
*/
void 
list_cpy(oldlp, newlp) 
	list * oldlp, *newlp; 
{ 
	link * oldp, *newp; 

	list_reset(newlp); 
	oldp = &(oldlp->first); 
	newp = &(newlp->first); 
	while (oldp->next) { 
		link *linkp; 
		auto data_ptr (*createfnp) PARMS ((VOID)); 
		auto void (*copyfnp) PARMS ((dp,dp)); 
		linkp = valloc1(link); 

		createfnp = newlp->fntab->create; 
		copyfnp = oldlp->fntab->copy; 
		linkp->dp = (*createfnp)(); 
		if (oldp->next->dp)
			(*copyfnp)(oldp->next->dp,linkp->dp);  
		linkp->next = 0; 
		newp->next = linkp; 	
		newlp->lastp = linkp; 
		oldp = oldp->next; 
		newp = newp->next; /* move forward in both lists */
	} 
} 


void 
list_mv(oldp, newp) 
	list * oldp, * newp; 
{ 
	assert(oldp != newp);
	assert (oldp->fntab == newp->fntab);
        newp = oldp;
	
} 
 
