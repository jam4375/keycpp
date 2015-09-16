IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -std=c++11 -O3 -g -Wall -Wextra -pedantic

ODIR=obj
LDIR =../lib

ifeq ($(lib),mkl)
        LIBS=-lm -L/$(MKLROOT)/lib/intel64 -lmkl_rt
else
        LIBS=-lm -lblas -llapack -larpack
endif

_DEPS = Figure.h Spline.h keycpp.h vector_k.h Matrix.h znaupd.h kiss_fft.h gnuplot_i.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = example.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS) ; $(CC) -c -o $@ $< $(CFLAGS)

example: $(OBJ) ; $(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS)

.PHONY: clean

clean: ; rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ ./example
