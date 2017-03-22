CC		= g++
#CC		= /usr/bin/clang++
PROG	= TumorDyn
#PROG	= HCL_Solve_test
PROG_O	= TumorDyn.o
EXC		= ./bin

VPATH	= ./main ./src

IDIR	= ./include
LIBS	= -llapack -lblas
CFLAGS	= -I$(IDIR) -O3 -Wall 
#CFLAGS	= -I$(IDIR) -O3 -Wall -g
#CFLAGS = -I$(IDIR) -O3 -Wall -Wno-overloaded-virtual
#CFLAGS  = -pg -I$(IDIR)

ODIR	= ./obj
PDIR	= ./main

export GLIBCXX_FORCE_NEW

_LDEPS = Data_Type.h LeastSquares_Solve.h LowDim_Alg.h HighDim_Alg.h\
			   Mat_Entry.h Alg_Solve.h Comput_Choices.h Alloc_Utils.h\
				 Quad_Coef.h 
_DEPS = TDInput.h TDSolve.h
DEPS	= $(patsubst %,$(IDIR)/%,$(_DEPS))


_OBJ   = LeastSquares_Solve.o LowDim_Alg.o HighDim_Alg.o\
				 Mat_Entry.o Alg_Solve.o Comput_Choices.o Alloc_Utils.o\
				 Quad_Coef.o TDInput.o TDSolve.o
OBJ		= $(patsubst %,$(ODIR)/%,$(_OBJ)) $(ODIR)/$(PROG_O)

$(PROG):$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o:%.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PDIR)/%.o:%.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: cleanall
cleanall: 
	rm $(PROG) $(PDIR)/$(PROG_O) $(ODIR)/*.o

.PHONY: clean
clean:
	rm $(PROG) $(PDIR)/$(PROG_O)

.PHONY: cleandata
cleandata:
	rm output*
