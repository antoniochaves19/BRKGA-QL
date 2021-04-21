all: runTest

CCC = g++
CCFLAGS = -O3 -openmp -Wall -std=c++11
TARGETDIR=.

OBJS2 = \
	$(TARGETDIR)/BRKGA_QL.o \
        $(TARGETDIR)/Read.o \
	$(TARGETDIR)/Decoder.o \
        $(TARGETDIR)/LocalSearch.o

$(TARGETDIR)/runTest: $(OBJS2)
	$(CCC) $(CCFLAGS) -o $(TARGETDIR)/runTest $(OBJS2)
	
$(TARGETDIR)/Read.o: Read.h Read.cpp
	$(CCC) $(CCFLAGS) -c Read.cpp -o $(TARGETDIR)/Read.o

$(TARGETDIR)/Decoder.o: Decoder.h Decoder.cpp
	$(CCC) $(CCFLAGS) -c Decoder.cpp -o $(TARGETDIR)/Decoder.o

$(TARGETDIR)/LocalSearch.o: LocalSearch.h LocalSearch.cpp
	$(CCC) $(CCFLAGS) -c LocalSearch.cpp -o $(TARGETDIR)/LocalSearch.o

$(TARGETDIR)/BRKGA_QL.o: BRKGA_QL.h BRKGA_QL.cpp
	$(CCC) $(CCFLAGS) -c BRKGA_QL.cpp -o $(TARGETDIR)/BRKGA_QL.o

test: runTest
	./runTest

clean:
	$(RM) \
	$(TARGETDIR)/BRKGA_QL.o \
	$(TARGETDIR)/Read.o \
	$(TARGETDIR)/Decoder.o \
   	$(TARGETDIR)/LocalSearch.o 
     
