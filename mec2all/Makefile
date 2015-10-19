FC = gfortran

OBJ = mec2all

SUBS =  an2dsr.o an2mom.o antpin.o dsrin.o truth.o iyesno.o cstring.o \
        fmreps.o mt_in.o printx.o pttpin.o  trpl2v.o v2trpl.o getstring.o \
		 timdat.o

all: ${OBJ}

%.o: %.f90
		$(FC) -c $^

mec2all:%: %.o  $(SUBS)
		$(FC) $^ -o $@

clean:
		rm -f $(OBJ) *.o
