FC      = gfortran
FFLAGS  = -O2 -Wall -Wunused

TARGET  = muon_range

OBJS    = randn.o      \
          corset.o     \
          corgen.o     \
          music-sr.o   \
          muon_range.o

ECOMUG_OBJS = randn.o ecomug.o test_ecomug.o

PROP_OBJS   = randn.o corset.o corgen.o music-sr.o ecomug.o muon_prop_200m.o

all: $(TARGET) test_ecomug muon_prop_200m

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

test_ecomug: $(ECOMUG_OBJS)
	$(FC) $(FFLAGS) -o $@ $(ECOMUG_OBJS)

muon_prop_200m: $(PROP_OBJS)
	$(FC) $(FFLAGS) -o $@ $(PROP_OBJS)

randn.o:      randn.f90
corset.o:     corset.f90
corgen.o:     corgen.f90
music-sr.o:   music-sr.f90
muon_range.o: muon_range.f90 music-sr.o
ecomug.o:          ecomug.f90
test_ecomug.o:     test_ecomug.f90 ecomug.o
muon_prop_200m.o:  muon_prop_200m.f90 ecomug.o music-sr.o

# corgen depends on corset (uses the Cholesky factor it produces)
corgen.o: corset.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: all clean
clean:
	rm -f $(OBJS) $(ECOMUG_OBJS) $(PROP_OBJS) $(TARGET) \
	      test_ecomug muon_prop_200m ecomug_mod.mod
