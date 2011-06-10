CXX=g++
MPI=no

CFLAGS= -Wall -g -I$(LIBDIR) $(OTHERCFLAGS) # -O3
#CFLAGS=-Wall -g -O3
LDFLAGS = -L$(LIBDIR) -l$(LIBNAME) -ldes -lm -lstdc++
#OTHER_CFLAGS = -I$(LIBDIR)

SOURCES = ParticleCommunity.cpp ParticleIntegrator.cpp FitnessLandscape.cpp \
	Genotype.cpp LValueWriter.cpp LParameters.cpp LSimulation.cpp
OBJS = $(SOURCES:.cpp=.o)

EXE = hiv
LIBNAME = adap-dyn
LIBDIR = ../../adap-dyn
LIBSTEM = lib$(LIBNAME).a
LIB = $(LIBDIR)/$(LIBSTEM)

default: $(EXE)

$(EXE): $(OBJS) $(LIB)
	$(CXX) $(CFLAGS) -o $(EXE) $(OBJS) $(LDFLAGS)

$(LIB): dummy
	cd $(LIBDIR) && $(MAKE) 'CFLAGS=$(CFLAGS)' $(LIBSTEM)

test-des: test-des.o $(OBJS) $(LIB)
	$(CXX) $(CFLAGS) -o test-des test-des.o $(OBJS) $(LDFLAGS)

test-rand: test-rand.o $(OBJS) $(LIB)
	$(CXX) $(CFLAGS) -o test-rand test-rand.o $(OBJS) $(LDFLAGS)

dot:
	@sh -c 'for d in out/0/*.dot; do echo dot -Tps -o out/$${d#out/0/}.ps $$d; dot -Tps -o out/$${d#out/0}.ps $$d; done'

#dot: out/landscape.ps out/quasispecies.ps
#out/landscape.ps: out/0/landscape.dot
#	dot -Tps -o out/landscape.ps out/0/landscape.dot
#out/0/landscape.dot: 

out/quasispecies.ps: out/0/quasispecies.dot
	dot -Tps -o out/quasispecies.ps out/0/quasispecies.dot

run: $(EXE) clear
	$(EXE)

clear:
	$(RM) -r out/*

clean:
	$(RM) -r $(ALLOBJS) $(EXE) $(FLOTSAM)

tags TAGS:
	etags *.{h,cpp}

fresh: clean libclean default

libclean: dummy
	+(cd $(LIBDIR); $(MAKE) clean)

dummy:

#.cc.o:
#	$(CXX) $(CFLAGS) -c -o $@ $<

#.cpp.o:
%.o : %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

#ALLOBJS = $(OBJS) $(OTHER_OBJS) $(RAND_OBJS)
ALLOBJS = *.o

FLOTSAM = core test-des *.idb *.pdb *.exe *.ilk *.obj *~ *\# #out/*

# header file dependencies (surely incomplete)
Genotype.o::Genotype.h LParameters.h
FitnessLandscape.o::FitnessLandscape.h Genotype.h LParameters.h
ParticleIntegrator.o::ParticleIntegrator.h FitnessLandscape.h LParameters.h \
	Genotype.h
ParticleCommunity.o::ParticleCommunity.h Genotype.h
LValueWriter.o::LValueWriter.h LParameters.h Genotype.h
LParameters.o::LParameters.h Genotype.h
main.o::LParameters.h Genotype.h
test-des.o::LParameters.h Genotype.h
