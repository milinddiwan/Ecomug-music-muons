FC      = gfortran
FFLAGS  = -O2 -Wall -Wunused

TARGET  = muon_range

OBJS    = randn.o      \
          corset.o     \
          corgen.o     \
          music-sr.o   \
          muon_range.o

ECOMUG_OBJS = ecomug.o test_ecomug.o

all: $(TARGET) test_ecomug

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

test_ecomug: $(ECOMUG_OBJS)
	$(FC) $(FFLAGS) -o $@ $(ECOMUG_OBJS)

randn.o:      randn.f90
corset.o:     corset.f90
corgen.o:     corgen.f90
music-sr.o:   music-sr.f90
muon_range.o: muon_range.f90 music-sr.o
ecomug.o:     ecomug.f90
test_ecomug.o: test_ecomug.f90 ecomug.o

# corgen depends on corset (uses the Cholesky factor it produces)
corgen.o: corset.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: all clean
clean:
	rm -f $(OBJS) $(ECOMUG_OBJS) $(TARGET) test_ecomug ecomug_mod.mod
