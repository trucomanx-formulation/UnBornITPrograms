NOMBREPROGRAMA=ceo-binary-simulator
CFLAGS=-c -Wall

all: $(NOMBREPROGRAMA) mutua
	rm -rf *.o

################################################################################

mutua: extrasbvector.o ceobinary.o extrascoin.o extras.o mutua.o 
	gcc -o mutua extrasbvector.o ceobinary.o extrascoin.o extras.o mutua.o -lpdsrv -lpdsra -lpdsba -lpdsdatafunc -lm
mutua.o: mutua.c extras/extras.h
	gcc -o mutua.o mutua.c $(CFLAGS)

################################################################################
$(NOMBREPROGRAMA): extrasbvector.o ceobinary.o extrascoin.o extras.o main.o 
	gcc -o $(NOMBREPROGRAMA) extrasbvector.o ceobinary.o extrascoin.o extras.o main.o -lpdsrv -lpdsra -lpdsba -lpdsdatafunc -lm

main.o: main.c extras/extras.h
	gcc -o main.o main.c $(CFLAGS)

ceobinary.o: extras/ceobinary.c extras/ceobinary.h
	gcc -o ceobinary.o extras/ceobinary.c $(CFLAGS)

extrascoin.o: extras/extrascoin.c extras/extrascoin.h
	gcc -o extrascoin.o extras/extrascoin.c $(CFLAGS)

extrasbvector.o: extras/extrasbvector.c extras/extrasbvector.h
	gcc -o extrasbvector.o extras/extrasbvector.c $(CFLAGS)

extras.o: extras/extras.c extras/extras.h
	gcc -o extras.o extras/extras.c $(CFLAGS)

clean:
	rm -rf *.o $(NOMBREPROGRAMA) mutua 

