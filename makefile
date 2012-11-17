CC = g++
CFLAGS = -o -m64 -lm -Wall -pthread

install: seqtools2

clean:
	rm *.o seqtools2

seqtools2: main.cpp main.h general_utils.o
	$(CC) $(CFLAGS) -o seqtools2 main.cpp general_utils.o -lz

general_utils.o: general_utils.cpp general_utils.h
	 $(CC) $(CFLAGS) -c general_utils.cpp




