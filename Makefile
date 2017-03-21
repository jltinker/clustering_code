hd = $(HOME)/lib
LIB = -lm -L${hd} -lcutil -fopenmp


CC = gcc
CFLAGS = -O2 -fopenmp

OBJS1721 = wp_LS_weight.o  qromo.o midpnt.o polint.o meshlink2.o nbrsfind2.o time.o sort2.o
wp_LS_weight_omp:	$(OBJS1721)
	$(CC) -o $@ $(OBJS1721) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS1722 = ximulti_LS_weight_omp.o  qromo.o midpnt.o polint.o meshlink2.o nbrsfind2.o time.o sort2.o
ximulti_LS_weight_omp:	$(OBJS1722)
	$(CC) -o $@ $(OBJS1722) $(LIB)
	cp -f $@ $(HOME)/exec/$@

OBJS1723 = wp_LS_weight_xcorr_omp.o  qromo.o midpnt.o polint.o meshlink2.o nbrsfind2.o time.o sort2.o
wp_LS_weight_xcorr_omp:	$(OBJS1723)
	$(CC) -o $@ $(OBJS1723) $(LIB)
	cp -f $@ $(HOME)/exec/$@

clean:
	rm -f *.o
