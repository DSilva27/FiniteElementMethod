IDIR = ../include
CC = g++
CFLAGS = -I$(IDIR)

ODIR = ../tools

_DEPS = DblInt.h finite_element.h LinAlg.h LInt.h triangle.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = DblInt.cpp LinAlg.cpp LInt.cpp functions_triang.cpp
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: main.cpp functions.cpp $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
