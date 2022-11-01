CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11 -g

Esercizio_10.exe : main.o random.o function.o
	$(CC) function.o random.o main.o -o Esercizio_10.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
function.o : function.cpp function.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio_10.exe seed.out
erase:
	sh clean.sh
