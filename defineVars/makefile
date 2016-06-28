CTYPE = g++
CFLAGS = -O3 -static
DefineVars: DefineVars.o
	$(CTYPE) $(CFLAGS) DefineVars.o

DefineVars.o: DefineVars.cpp DefineVars.h
	$(CTYPE) -c $(CFLAGS) DefineVars.cpp
