# C++ compiler
CXX=clang++
# GSL >=2.0 required
GSLCFLAGS=`gsl-config --cflags`
GSLLIBS=`gsl-config --libs`
# C++ compiler flags -I/usr/local/opt/llvm/include -fopenmp
CXXFLAGS =-Xpreprocessor -fopenmp -Wall -std=c++11 -O3 -lstdc++ -I/usr/local/opt/llvm/include -lomp
# ini file parser files
INISRC =\
		 inih/ini.c \
		 inih/INIReader.cpp
INIINC=-I./inih
# source files of KoMPoST
KOMPOSTSRC =\
		 src/Main.cpp \
		 src/EventInput.cpp \
		 src/BackgroundEvolution.cpp \
		 src/GreensFunctions.cpp \
		 src/KineticEvolution.cpp \
		 src/ScalingVariable.cpp

all: 
	$(CXX) -o KoMPoST.exe $(KOMPOSTSRC) $(INISRC) $(INIINC) $(CXXFLAGS) $(GSLCFLAGS) $(GSLLIBS) 

.PHONY: clean
clean:
	rm KoMPoST.exe
