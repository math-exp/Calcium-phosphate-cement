# Time-stamp: <2022-10-31 17:41:08 shige>
# Makefile

CC = gcc
#CFLAGS = -O4
#CFLAGS = -g -static

BIN = paste
SRC = $(BIN).c head_$(BIN).c ../tools/matutil.c
OBJ = $(BIN).o head_$(BIN).o matutil.o

LIBS = -lm

all: $(BIN)

debug: CFLAGS= -ggdb

debug: $(BIN)

$(BIN) : $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) -o $(BIN) $(LIBS)

$(OBJ) : $(SRC)
	$(CC) $(CFLAGS) -c $(SRC) 

clean :
#	cp *.c *.h ../data.prm ../do.sh Makefile ../backup
	cp *.c *.h ../do.sh Makefile ../backup
	rm -f $(OBJ) $(BIN)
	rm -f *~ *# ../*~
