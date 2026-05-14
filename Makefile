FC      = gfortran
FFLAGS  = -O2 -Wall -Wunused

TARGET  = muon_range

OBJS    = randn.o      \
          corset.o     \
          corgen.o     \
          music-sr.o   \
          muon_range.o

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

randn.o:      randn.f90
corset.o:     corset.f90
corgen.o:     corgen.f90
music-sr.o:   music-sr.f90
muon_range.o: muon_range.f90 music-sr.o

# corgen depends on corset (uses the Cholesky factor it produces)
corgen.o: corset.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
