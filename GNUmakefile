CXX=g++
MPI=no

CFLAGS= -Wall -g -I$(LIBDIR) -I$(LIBDIR)/Displays $(OTHERCFLAGS) # -O3
#CFLAGS=-Wall -g -O3
LDFLAGS = -L$(LIBDIR) -l$(LIBNAME) -L$(VXLDIR)/lib -lvnl_algo -lvnl -lvcl -lnetlib -lssl -lm -lstdc++
#OTHER_CFLAGS = -I$(LIBDIR)
VXLDIR ?= ../vxl

SOURCES = Collective.cpp FitnessLandscape.cpp Individual.cpp \
	BitString.cpp LParameters.cpp ConsensusSimulation.cpp LOutputController.cpp
OBJS = $(SOURCES:.cpp=.o)

EXE = consensus
LIBNAME = adap-dyn
LIBDIR ?= ../adap-dyn
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

#out/quasispecies.ps: out/0/quasispecies.dot
#	dot -Tps -o out/quasispecies.ps out/0/quasispecies.dot

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

%.o : %.cpp
	$(CXX) $(CFLAGS) -MD -c -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.o : %.c
	$(CC) $(CFLAGS) -MD -c $< -o $@
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

#ALLOBJS = $(OBJS) $(OTHER_OBJS) $(RAND_OBJS)
ALLOBJS = *.o

FLOTSAM = core test-des *.idb *.pdb *.exe *.ilk *.obj *~ *\# #out/*

PXXS = $(SOURCES:.cpp=.P)
PFILES = $(PXXS:.c=.P)
-include $(PFILES)

