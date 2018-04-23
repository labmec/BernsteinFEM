# Makefile

# compiler options
OBJ= -c
FLAGS= -g -std=c++17 -lm -Wall -O3

# include directories
INC_ARMA=-I../Armadillo/include
INC_MOM=$(INC_ARMA) -IMoments/
INC_MASS=$(INC_ARMA) -IMass/
INC_STIFF=$(INC_ARMA) -IStiffness/
INC_ELEM=$(INC_ARMA) -IElements/
INC_QUADRA=$(INC_ARMA) -IQuadra/
INC_ALL=$(INC_MOM) $(INC_MASS) $(INC_STIFF) $(INC_QUADRA) $(INC_ELEM)

# library name
LIB=BernsteinFEM.a

all: quadra moments mass stiff elem
	$(MAKE) lib
	$(MAKE) clean

lib:
	ar -r $(LIB) *.o

clean:
	rm -f *.o

quadra:
	g++ $(OBJ) $(FLAGS) Quadra/JacobiGaussNodes.cpp 

moments: moments1d moments2d

moments1d:
	g++ $(OBJ) $(FLAGS) $(INC_QUADRA) Moments/Moments1D.cpp

moments2d: moments2dtri moments2dquad

moments2dtri:
	g++ $(OBJ) $(FLAGS) $(INC_QUADRA) Moments/Moments2DTri.cpp

moments2dquad:
	g++ $(OBJ) $(FLAGS) $(INC_QUADRA) Moments/Moments2DQuad.cpp

mass: mass1d mass2d

mass1d:
	g++ $(OBJ) $(FLAGS) $(INC_MOM) Mass/MassM1D.cpp

mass2d:	mass2dtri mass2dquad

mass2dtri:
	g++ $(OBJ) $(FLAGS) $(INC_MOM) Mass/MassM2DTri.cpp

mass2dquad:
	g++ $(OBJ) $(FLAGS) $(INC_MOM) Mass/MassM2DQuad.cpp

stiff: stiff1d stiff2d

stiff1d:
	g++ $(OBJ) $(FLAGS) $(INC_MOM) Stiffness/StiffM1D.cpp

stiff2d: stiff2dtri stiff2dquad

stiff2dtri:
	g++ $(OBJ) $(FLAGS) $(INC_MOM) Stiffness/StiffM2DTri.cpp

stiff2dquad:
	g++ $(OBJ) $(FLAGS) $(INC_MOM) Stiffness/StiffM2DQuad.cpp

elem: elem1d elem2d

elem1d:
	g++ $(OBJ) $(FLAGS) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem1D.cpp

elem2d:	elem2dtri elem2dquad

elem2dtri:
	g++ $(OBJ) $(FLAGS) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem2DTri.cpp

elem2dquad:
	g++ $(OBJ) $(FLAGS) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem2DQuad.cpp