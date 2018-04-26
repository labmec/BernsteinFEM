# Makefile

# compiler options
OBJ= -c
FLAGS= -std=c++17 -lm -Wall

# choose either you want debug or optmization on compilation
OPT_DBG = -g #-O3

# include directories
INC_MOM=-IMoments/
INC_MASS=-IMass/
INC_STIFF=-IStiffness/
INC_ELEM=-IElements/
INC_QUADRA=-IQuadra/
INC_ALL=$(INC_MOM) $(INC_MASS) $(INC_STIFF) $(INC_QUADRA) $(INC_ELEM)

# library name
LIB=BernsteinFEM.a

all: quadra moments mass stiff elem
	$(MAKE) lib

lib:
	ar -r $(LIB) *.o

clean:
	rm -f *.o

quadra:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) Quadra/JacobiGaussNodes.cpp 

moments: moments1d moments2d

moments1d:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_QUADRA) Moments/Moments1D.cpp

moments2d: moments2dtri moments2dquad

moments2dtri:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_QUADRA) Moments/Moments2DTri.cpp

moments2dquad:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_QUADRA) Moments/Moments2DQuad.cpp

mass: mass1d mass2d

mass1d:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MOM) Mass/MassM1D.cpp

mass2d:	mass2dtri mass2dquad

mass2dtri:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MOM) Mass/MassM2DTri.cpp

mass2dquad:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MOM) Mass/MassM2DQuad.cpp

stiff: stiff1d stiff2d

stiff1d:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MOM) Stiffness/StiffM1D.cpp

stiff2d: stiff2dtri stiff2dquad

stiff2dtri:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MOM) Stiffness/StiffM2DTri.cpp

stiff2dquad:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MOM) Stiffness/StiffM2DQuad.cpp

elem: elem1d elem2d

elem1d:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem1D.cpp

elem2d:	elem2dtri elem2dquad

elem2dtri:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem2DTri.cpp

elem2dquad:
	g++ $(OBJ) $(FLAGS)  $(OPT_DBG) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem2DQuad.cpp