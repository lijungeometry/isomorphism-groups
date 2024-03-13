SHELL = /bin/sh
CFLAGS= -O2 -I/usr/local/include 
.c.o:
	${CC} -c $(CFLAGS) lil.c 
${BINDIR}/lil: lil.o  
	${CC} $(CFLAGS)  -o ${BINDIR}/lil lil.o 
