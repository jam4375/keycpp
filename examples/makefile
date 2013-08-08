IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -std=c++11 -O3

ODIR=obj
LDIR =../lib

LIBS=-lm -lkeycpp -llapack

_DEPS = keycpp.h Matrix.h kiss_fft.h gnuplot_i.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS)) ../lib/libkeycpp.a

_OBJ = example.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS) ; $(CC) -c -o $@ $< $(CFLAGS)

example: $(OBJ) ; rm -f /tmp/gnuplot* ; $(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS)

.PHONY: clean

clean: ; rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ ./example
