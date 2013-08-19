IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -std=c++11 -O3

ODIR=obj
LDIR =

_DEPS = Figure.h Spline.h keycpp.h Matrix.h kiss_fft.h gnuplot_i.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gnuplot_i.o kiss_fft.o keycpp.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS) ; $(CC) -c -o $@ $< $(CFLAGS)

libkeycpp: $(OBJ) ; ar crf ../lib/libkeycpp.a $^

.PHONY: clean

clean: ; rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ ../lib/libkeycpp.a core $(INCDIR)/*~
