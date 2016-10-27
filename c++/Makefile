#/******************************************************************************
# * Version: 1.0.1
# * Last modified on: 27	October, 2016 
# * Developers: Michael G. Epitropakis, Xiaodong Li.
# *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
# *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
# * ***************************************************************************/
CFLAGS = -g -Wall -Weffc++ -I./ 
CC = g++
LDFLAGS = 
FILES = cec2013.o cfunction.o rand2.o

all: $(FILES) main.cpp
	$(CC) $(CFLAGS) main.cpp $(LDFLAGS) $(FILES)

plots: $(FILES) plots.cpp
	$(CC) $(CFLAGS) plots.cpp $(FILES) -o plots.out $(LDFLAGS)

cec2013.o: cec2013.h cec2013.cpp cfunction.h cfunction.cpp
	$(CC) $(CFLAGS) -c cec2013.cpp

cfunction.o: cfunction.h cfunction.cpp
	$(CC) $(CFLAGS) -c cfunction.cpp

rand2.o: rand2.c rand2.h
	$(CC) $(CFLAGS) -c rand2.c

clean:
	-rm -f *.o *~ core a.out 
