# lib.d  makefile 
# 22. 7. 91.
SHELL=/bin/sh
CFLAGS= -O2 $(DEBUG)
LIB = isomlib.a
LIBOBJ = afsa.o input.o list.o permiofns.o permfns.o reduce.o storage.o word.o
$(LIB): $(LIBOBJ)
	ar cr $(LIB) $(LIBOBJ)
	if  which ranlib > /dev/null 2> /dev/null; then ranlib $(LIB); fi
.c.o:
	${CC} -c $(CFLAGS) $*.c 
clean:	
	\rm *.o
