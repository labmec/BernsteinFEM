# Makefile

# compiler options
OBJ= -c
GFLAGS= -ggdb -std=c++17 -lm -Wall
OFLAGS= -std=c++17 -lm -Wall -O3

# object files
ELEM_OBJS= Elem1D.o Elem2DQuad.o Elem2DTri.o
MASS_OBJS= MassM1D.o MassM2DTri.o MassM2DQuad.o
QUADRA_OBJS= JacobiGaussNodes.o
MOM_OBJS= Moments1D.o Moments2DQuad.o Moments2DTri.o
STIFF_OBJS= StiffM1D.o StiffM2DTri.o StiffM2DQuad.o
DER_OBJS= QuadD.o TriD.o
OBJECTS= $(ELEM_OBJS) $(MASS_OBJS) $(QUADRA_OBJS) $(MOM_OBJS) $(STIFF_OBJS) $(DER_OBJS)

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

.PHONY : all
all: $(LIB)

BernsteinFEM.a: $(OBJECTS)
	ar -r $(LIB) $(OBJECTS)

clean:
	rm -f $(OBJECTS)

JacobiGaussNodes.o: Quadra/JacobiGaussNodes.cpp Quadra/JacobiGaussNodes.h
	g++ $(OBJ) $(GFLAGS) Quadra/JacobiGaussNodes.cpp

Moments1D.o: Moments/Moments1D.cpp Moments/Moments.h
	g++ $(OBJ) $(GFLAGS) $(INC_QUADRA) Moments/Moments1D.cpp

Moments2DTri.o: Moments/Moments2DTri.cpp Moments/Moments.h
	g++ $(OBJ) $(GFLAGS) $(INC_QUADRA) Moments/Moments2DTri.cpp

Moments2DQuad.o: Moments/Moments2DQuad.cpp Moments/Moments.h
	g++ $(OBJ) $(GFLAGS) $(INC_QUADRA) Moments/Moments2DQuad.cpp

MassM1D.o: Mass/MassM1D.cpp Mass/MassM.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) Mass/MassM1D.cpp

MassM2DTri.o: Mass/MassM2DTri.cpp Mass/MassM.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) Mass/MassM2DTri.cpp

MassM2DQuad.o: Mass/MassM2DQuad.cpp Mass/MassM.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) Mass/MassM2DQuad.cpp

StiffM1D.o: Stiffness/StiffM1D.cpp Stiffness/StiffM.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) Stiffness/StiffM1D.cpp

StiffM2DTri.o: Stiffness/StiffM2DTri.cpp Stiffness/StiffM.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) Stiffness/StiffM2DTri.cpp

StiffM2DQuad.o: Stiffness/StiffM2DQuad.cpp Stiffness/StiffM.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) Stiffness/StiffM2DQuad.cpp

QuadD.o: Derivatives/QuadD.cpp Derivatives/Derivatives.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) $(INC_QUADRA) Derivatives/QuadD.cpp

TriD.o: Derivatives/TriD.cpp Derivatives/Derivatives.h
	g++ $(OBJ) $(GFLAGS) $(INC_MOM) $(INC_QUADRA) Derivatives/TriD.cpp

Elem1D.o: Elements/Elem1D.cpp Elements/Elem.h
	g++ $(OBJ) $(GFLAGS) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem1D.cpp

Elem2DTri.o: Elements/Elem2DTri.cpp Elements/Elem.h
	g++ $(OBJ) $(GFLAGS) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem2DTri.cpp

Elem2DQuad.o: Elements/Elem2DQuad.cpp Elements/Elem.h
	g++ $(OBJ) $(GFLAGS) $(INC_MASS) $(INC_STIFF) $(INC_MOM) $(INC_QUADRA) Elements/Elem2DQuad.cpp