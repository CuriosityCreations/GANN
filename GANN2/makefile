CTYPE = g++
CFLAGS = -O3 -static

TOPH = DCGA.h classErr.h
GAH = $(TOPH) classSelectableBase.h \
	classEvolvable.h classChromosome.h classGene.h

SMALLOBJS = DCGA.o
GAOBJS = classGene.o classChromosome.o classEvolvable.o classSelectableBase.o
PROGOBJS = partitions.o OGA.o MLP.o oCscore.o

ALLOBJS = $(PROGOBJS) $(SEQOBJS) $(GAOBJS) $(SMALLOBJS)

GAGAexec: $(ALLOBJS)
	$(CTYPE) $(CFLAGS) $(ALLOBJS)

partitions.o: partitions.cpp OGA.h MLP.h oCscore.h $(GAH)
	$(CTYPE) -c $(CFLAGS) partitions.cpp

OGA.o: OGA.cpp OGA.h MLP.h oCscore.h $(GAH)
	$(CTYPE) -c $(CFLAGS) OGA.cpp

MLP.o: MLP.cpp MLP.h oCscore.h $(GAH)
	$(CTYPE) -c $(CFLAGS) MLP.cpp

oCscore.o: oCscore.cpp oCscore.h $(GAH)
	$(CTYPE) -c $(CFLAGS) oCscore.cpp

classSelectableBase.o: classSelectableBase.cpp classSelectableBase.h \
	classEvolvable.h classChromosome.h classGene.h $(TOPH)
	$(CTYPE) -c $(CFLAGS) classSelectableBase.cpp

classEvolvable.o: classEvolvable.cpp classEvolvable.h \
	classChromosome.h classGene.h $(TOPH)
	$(CTYPE) -c $(CFLAGS) classEvolvable.cpp

classChromosome.o: classChromosome.cpp classChromosome.h \
	classGene.h $(TOPH)
	$(CTYPE) -c $(CFLAGS) classChromosome.cpp

classGene.o: classGene.cpp classGene.h $(TOPH)
	$(CTYPE) -c $(CFLAGS) classGene.cpp

DCGA.o: DCGA.cpp DCGA.h
	$(CTYPE) -c $(CFLAGS) DCGA.cpp
