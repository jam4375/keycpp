IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -std=c++11 -O3

ODIR=obj
LDIR =

LIBS=-lm -llapack

_DEPS = Figure.h Spline.h keycpp.h Matrix.h kiss_fft.h gnuplot_i.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gnuplot_i.o kiss_fft.o keycpp.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS) ; $(CC) -c -o $@ $< $(CFLAGS)

libkeycpp: $(OBJ) ; ar -x /usr/lib/lapack/liblapack.a; ar -x /usr/lib/libblas.a; ar -x /usr/lib/gcc/x86_64-linux-gnu/4.7/libgfortran.a; ar crf ../lib/libkeycpp.a $^ *.o; rm *.o

.PHONY: clean

clean: ; rm -f *.o $(ODIR)/*.o *~ core $(INCDIR)/*~ ../lib/libkeycpp.a core $(INCDIR)/*~
