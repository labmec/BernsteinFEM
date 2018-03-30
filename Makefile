# Makefile

OBJ= -c
FLAGS= -g -std=c++17 -lm

all: moments mass stiff elem lib

lib:
	ar -r BernsteinFEM.a *.o

clean:
	rm -f *.o

moments: moments1d moments2d

moments1d:
	g++ $(OBJ) $(FLAGS) Moments1D.cpp

moments2d: moments2dtri moments2dquad

moments2dtri:
	g++ $(OBJ) $(FLAGS) Moments2DTri.cpp

moments2dquad:
	g++ $(OBJ) $(FLAGS) Moments2DQuad.cpp

mass: mass1d mass2d

mass1d:
	g++ $(OBJ) $(FLAGS) MassM1D.cpp

mass2d:	mass2dtri mass2dquad

mass2dtri:
	g++ $(OBJ) $(FLAGS) MassM2DTri.cpp

mass2dquad:
	g++ $(OBJ) $(FLAGS) MassM2DQuad.cpp

stiff: stiff1d stiff2d

stiff1d:
	g++ $(OBJ) $(FLAGS) StiffM1D.cpp

stiff2d: stiff2dtri stiff2dquad

stiff2dtri:
	g++ $(OBJ) $(FLAGS) StiffM2DTri.cpp

stiff2dquad:
	g++ $(OBJ) $(FLAGS) StiffM2DQuad.cpp

elem: elem1d elem2d

elem1d:
	g++ $(OBJ) $(FLAGS) Elem1D.cpp

elem2d:	elem2dtri elem2dquad

elem2dtri:
	g++ $(OBJ) $(FLAGS) Elem2DTri.cpp

elem2dquad:
	g++$(OBJ) $(FLAGS) Elem2DQuad.cpp