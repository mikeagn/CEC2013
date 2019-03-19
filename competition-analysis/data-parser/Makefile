# ******************************************************************************
# * Version: 1.0
# * Last modified on: 21 January, 2013 
# * Developers: Michael G. Epitropakis, Xiaodong Li.
# *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
# *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
# * ****************************************************************************
CFLAGS = -g -Wall -std=c++11 -I./ 
CC = g++
LDFLAGS = 
FILES = cec2013.o cfunction.o rand2.o params.o
RANDFILE=randfile

all: competition parse-static parse-dynamic


competition: $(FILES) main.cpp
ifeq ($(wildcard $(RANDFILE)),"")
	-ln -s /dev/urandom randfile
endif
	$(CC) $(CFLAGS) main.cpp $(LDFLAGS) $(FILES)

tex: $(FILES) create_tex.cpp
	$(CC) $(CFLAGS) create_tex.cpp $(LDFLAGS) $(FILES) -o ctex

parse-static: $(FILES) parse_submission_file.cpp
	$(CC) $(CFLAGS) parse_submission_file.cpp $(LDFLAGS) $(FILES) -o parse-static

parse-dynamic: $(FILES) parse_submission_file_dynamic.cpp
	$(CC) $(CFLAGS) parse_submission_file_dynamic.cpp $(LDFLAGS) $(FILES) -o parse-dynamic

plots: $(FILES) plots.cpp
	$(CC) $(CFLAGS) plots.cpp $(FILES) -o plots $(LDFLAGS)

cec2013.o: cec2013.h cec2013.cpp cfunction.h cfunction.cpp
	$(CC) $(CFLAGS) -c cec2013.cpp

cfunction.o: cfunction.h cfunction.cpp
	$(CC) $(CFLAGS) -c cfunction.cpp

rand2.o: rand2.c rand2.h
	$(CC) $(CFLAGS) -c rand2.c

params.o: params.h params.cpp
	$(CC) $(CFLAGS) -c params.cpp

clean:
	-rm -rf *.o *~ core a.out parse-dynamic parse-static ctex *.dSYM
