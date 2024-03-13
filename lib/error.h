/* lib.d file error.h */
/* If the compilation is done with the DEBUG flag turned on, then any
assert statement which is false leads to a message with filename, line
number and then a core dump.
*/
/*extern abort();*/
extern void error();
#ifdef DEBUG
#define assert(EX) {if ((EX)); else { \
	fprintf(stderr,"\nThe following assertion failed\n %s \nfilename: %s ",\
		"EX", __FILE__);\
	fprintf(stderr," line number %d\naborted\n", __LINE__) ;\
	fflush(stderr);\
	fflush(stdout);\
	abort();}}

#else
#define assert()
#endif
