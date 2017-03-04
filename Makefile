hd = $(HOME)/lib
LIB = -lm -L${hd} -lcutil -fopenmp


CC = gcc
CFLAGS = -O2 -fopenmp

OBJS1721 = wp_LS_weight.o  qromo.o midpnt.o polint.o meshlink2.o nbrsfind2.o time.o sort2.o
wp_LS_weight_omp:	$(OBJS1721)
	$(CC) -o $@ $(OBJS1721) $(LIB)
	cp -f $@ $(HOME)/exec/$@

clean:
	rm -f *.o
