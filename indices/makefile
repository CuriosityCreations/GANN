CTYPE = g++
CFLAGS = -O3 -static

INDexec: indices.o classGene.o classChromosome.o classBaseRep.o
	$(CTYPE) $(CFLAGS) indices.o classGene.o classChromosome.o classBaseRep.o

indices.o: indices.cpp indices.h classChromosome.h classGene.h classBaseRep.h
	$(CTYPE) -c $(CFLAGS) indices.cpp

classChromosome.o: classChromosome.cpp classChromosome.h classGene.h classBaseRep.h
	$(CTYPE) -c $(CFLAGS) classChromosome.cpp

classGene.o: classGene.cpp classGene.h classBaseRep.h
	$(CTYPE) -c $(CFLAGS) classGene.cpp

classBaseRep.o: classBaseRep.cpp classBaseRep.h
	$(CTYPE) -c $(CFLAGS) classBaseRep.cpp
