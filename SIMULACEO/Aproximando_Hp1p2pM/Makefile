NOMBREPROGRAMA=ahjoint
CFLAGS=-c -Wall

all: $(NOMBREPROGRAMA) 

$(NOMBREPROGRAMA): extras.o main.o 
	gcc -o $(NOMBREPROGRAMA) extras.o main.o -lpdsdatafunc -lpdsba -lpdsra -lpdsrv -lpdsm -lm -lpdsit

main.o: main.c extras.h
	gcc -o main.o main.c $(CFLAGS)

extras.o: extras.c extras.h
	gcc -o extras.o extras.c $(CFLAGS) 
clean:
	rm -rf *.o $(NOMBREPROGRAMA) 

