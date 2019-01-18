# Makefile

# compiler options
OBJ= -c
GFLAGS= -ggdb -std=c++11 -lm -Wall
OFLAGS= -std=c++11 -lm -Wall -O3
CXX= g++
CXX_FLAGS= $(GFLAGS) $(INC_ALL)

# object files
ELEM_OBJS= LinearEl.o QuadrilateralEl.o TriangularEl.o CubeEl.o TetrahedronEl.o
MASS_OBJS= MassM1D.o MassM2DTri.o MassM2DQuad.o BMass.o
QUADRA_OBJS= JacobiGaussNodes.o
MOM_OBJS= BMoment.o Moments1D.o Moments2DQuad.o Moments2DTri.o Moments3DCube.o Moments3DTetra.o
STIFF_OBJS= StiffM1D.o StiffM2DTri.o StiffM2DQuad.o BStiff.o
DER_OBJS= QuadD.o TriD.o
OBJECTS= $(MASS_OBJS) $(QUADRA_OBJS) $(MOM_OBJS) $(STIFF_OBJS) $(DER_OBJS) $(ELEM_OBJS)

# include directories
INC_ARMA=-I../Armadillo/include
INC_MOM=-IMoments/
INC_MASS=-IMass/
INC_STIFF=-IStiffness/
INC_ELEM=-IElements/
INC_QUADRA=-IQuadra/
INC_DERIV=-IDerivatives/
INC_ALL=$(INC_MOM) $(INC_MASS) $(INC_STIFF) $(INC_QUADRA) $(INC_ELEM) $(INC_DERIV) $(INC_ARMA)

# library name
LIB=BernsteinFEM.a


.PHONY : all
# default makes the static library (shared library is not available yet)
all: $(LIB)

# Make lib from all the object files
BernsteinFEM.a: $(OBJECTS)
	ar -r $(LIB) $(OBJECTS)

clean:
	rm -f $(OBJECTS)

# generic rule for the object files
%.o : */%.cpp
	$(CXX) $(CXX_FLAGS) -c $<