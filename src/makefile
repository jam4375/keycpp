IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -std=c++11 -O3 -Wall -Weffc++ -pedantic  \
    -pedantic-errors -Wextra  -Wall -Wcast-align \
    -Wcast-qual  -Wchar-subscripts  -Wcomment -Wconversion \
    -Wdisabled-optimization \
    -Werror -Wfloat-equal  -Wformat  -Wformat=2 \
    -Wformat-nonliteral -Wformat-security  \
    -Wformat-y2k \
    -Wimport  -Winit-self \
    -Winvalid-pch   \
    -Wunsafe-loop-optimizations  -Wlong-long -Wmissing-braces \
    -Wmissing-field-initializers -Wmissing-format-attribute   \
    -Wmissing-include-dirs -Wmissing-noreturn \
    -Wpacked -Wparentheses  -Wpointer-arith \
    -Wredundant-decls -Wreturn-type \
    -Wsequence-point  -Wshadow -Wsign-compare  -Wstack-protector \
    -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default \
    -Wswitch-enum -Wtrigraphs  -Wuninitialized \
    -Wunknown-pragmas  -Wunreachable-code -Wunused \
    -Wunused-function  -Wunused-label  -Wunused-parameter \
    -Wunused-value  -Wunused-variable  -Wvariadic-macros \
    -Wvolatile-register-var  -Wwrite-strings

ODIR=obj
LDIR =

_DEPS = Figure.h Spline.h keycpp.h vector_k.h Matrix.h kiss_fft.h gnuplot_i.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gnuplot_i.o kiss_fft.o keycpp.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS) ; $(CC) -c -o $@ $< $(CFLAGS)

libkeycpp: $(OBJ) ; ar crf ../lib/libkeycpp.a $^

.PHONY: clean

clean: ; rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ ../lib/libkeycpp.a core $(INCDIR)/*~
