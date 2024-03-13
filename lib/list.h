/* lib file list.h */
#define ORDERED 0
#define FIFO 1
#define LIFO 2
/* The structure of a |list|. We start off with some casts which will
be useful in element function tables.
*/
typedef	int		(*DPDP2I)			PARMS((dp,dp));
typedef	void	(*DPDP2V)			PARMS((dp,dp));
typedef	dp		(*V2DP)				PARMS((VOID));
typedef	void	(*DP2V)				PARMS((dp));
typedef	int		(*DP2I)				PARMS((dp));
typedef	void	(*DPI2V)			PARMS((dp, int));
typedef dp		(*DP2DP)			PARMS((dp));

typedef struct {
	int		(*signature)		PARMS((dp,dp));
	void	(*copy)				PARMS((dp,dp));
	dp		(*create)			PARMS((VOID));
	void	(*kill)				PARMS((dp));
	void	(*print)			PARMS((FILE *,dp));
}
elt_fntab; /*function table for manipulating the elements of a list*/
typedef dp		(*E2DP)				PARMS((elt_fntab*));

typedef struct link { 
	struct link * next;
	data_ptr dp;
} link;

/* The structure of a |list|. */
typedef struct list {
        char type; /* ORDERED, FIFO or LIFO */
	elt_fntab *fntab;/*address of table of functions which manipulate 
		structures stored in list*/
		/*Functions in the table need to be cast to return
			the correct kind of value. */
	link first;
	link * lastp;
} list;

/* The structure used to traverse a linked list.  */
typedef struct {
	elt_fntab * fntab; /*the function table of the list itself*/
	link * upto; /*space is not specially  allocated for this pointer. It
		is designed to move along the list, pointing at what is there
		already*/
} list_traverser;

extern void 
list_init PARMS((list *, elt_fntab *, char type));

extern void list_reset PARMS(( list *lp)); 
extern void list_clear PARMS(( list *lp)); 

/* Initialize a |list_traverser| to make it traverse a particular
list, whose address is given in the second argument.
*/
extern void list_traverser_init PARMS((list_traverser * ltp, list * lp)); 
extern void list_traverser_clear PARMS((list_traverser * ltp)); 

extern boolean list_next PARMS ((list_traverser*,dp));

/* The Boolean valued function |list_empty()| returns |TRUE| if the list is
empty, |FALSE| otherwise without changing the list in any way.
*/
extern boolean list_empty PARMS(( list * lp));

/* Is there an item on the list which is equal (by the lights of the
appropriate signature function) to the item pointed to by |in|? The
procedure returns |TRUE| or |FALSE|. If the answer is |TRUE|, then
|out| is made to point to a copy of the item on the list. |*out| does
not need to be a 100\% faithful copy of |*in|. It may only be a 95\%
faithful copy. Whether it is guaranteed to be
a faithful copy, or is not guaranteed to be faithful, depends on the
signature function and the copy function associated to this particular
list.
*/

/* |list_insert()| inserts an item into a list.
|lp| is the address of a list and |in| is the address of
the item. It returns |TRUE| if an item is inserted.
*/
extern boolean list_insert PARMS ((list*,dp));

/* Delete an item from the list if the signature function says it is
``equal'' to the item pointed to by |in|.
Returns |TRUE| if an item has been deleted, and |FALSE| otherwise.
|lp| is the address of a list and |in| is the address of
the item.
*/
extern boolean list_delete PARMS ((list*,dp));

/* Get a copy of the first item on the list, and put it into |*out|. The list
does not change. The procedure returns |FALSE| if the list is empty,
and |TRUE| otherwise.
|lp| is the address of a list and |in| is the address of
the item.
*/
extern boolean list_get_first PARMS(( list *lp,dp in));	

/* Get a copy of the first item of the list and put it into |*out|.
Remove the first item from the list.  Returns |TRUE| or |FALSE|, 
in the same way as |list_get_first()|.  |lp| is the address of a list 
and |in| is the address of the item.  */
extern boolean list_delget_first PARMS ((list*,dp out));

/* Print for ordinary use.
*/
extern void list_print PARMS ((FILE * wfile, list * lp, boolean comma));

/*Copy a list. The list |*oldp| is untouched and |*newp| is a copy of
|*oldp|. The original contents of |*newp| are destroyed.
Warning Do not try to copy a list to itself.
*/
extern void list_cpy PARMS(( list* oldp, list* newp));
/* The list *oldp is renamed *newp, and the original list *newp is destroyed
*/

