# Makefile

FLAGS=-c -std=c++17 -lm

all: moments mass stiff elem lib

lib:
	ar -r BernsteinFEM.a *.o

clean:
	rm -f *.o

moments: moments1d moments2d

moments1d:
	g++ $(FLAGS) Moments1D.cpp

moments2d: moments2dtri moments2dquad

moments2dtri:
	g++ $(FLAGS) Moments2DTri.cpp

moments2dquad:
	g++ $(FLAGS) Moments2DQuad.cpp

mass: mass1d mass2d

mass1d:
	g++ $(FLAGS) MassM1D.cpp

mass2d:	mass2dtri mass2dquad

mass2dtri:
	g++ $(FLAGS) MassM2DTri.cpp

mass2dquad:
	g++ $(FLAGS) MassM2DQuad.cpp

stiff: stiff1d stiff2d

stiff1d:
	g++ $(FLAGS) StiffM1D.cpp

stiff2d: stiff2dtri stiff2dquad

stiff2dtri:
	g++ $(FLAGS) StiffM2DTri.cpp

stiff2dquad:
	g++ $(FLAGS) StiffM2DQuad.cpp

elem: elem1d elem2d

elem1d:
	g++ $(FLAGS) Elem1D.cpp

elem2d:	elem2dtri elem2dquad

elem2dtri:
	g++ $(FLAGS) Elem2DTri.cpp

elem2dquad:
	g++$(FLAGS) Elem2DQuad.cpp