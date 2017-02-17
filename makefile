###########################################

FLIBS = ../../wmsc-lucamazzone/lectures/day5/code_day5/pardiso_example/lib/libpardiso500-GNU481-X86-64.so -llapack -lblas -fopenmp -lpthread -lm

GFORTRAN = gfortran -O3 -g 


##############################
all :  mainfile.exec
#############################


library.o : library.f95
	$(GFORTRAN)$ -c $< -o $@
	
solution_lib.o : solution_lib.f95
	$(GFORTRAN)$ -c $< -o $@
		
param.o : params.f95
	$(GFORTRAN)$ -c $< -o $@
			
main.o : MAIN.f95 param.o library.o solution_lib.o
	$(GFORTRAN)$  -c $< -o $@ 

				
#############################
# executable
mainfile.exec: main.o library.o solution_lib.o param.o
	$(GFORTRAN)$  $^ -o $@  $(FLIBS)
#############################








#FFLAGS = -g 
#
#Place of Sparse BLAS objects
#SB_LIB = -L../SOFTWARE -lSparseBLAS_GNU
#Place of Sparse BLAS modules
#SB_INCL = -I../SOFTWARE
#Place of numeric libraries
#SYS_LIB = 
# PARDISO
#FLIBS = -L../PARDISO ../PARDISO/libpardiso500-MACOS-X86-64.dylib -llapack -lblas -fopenmp -lpthread -lm


#FC = gfortran -O3 -g
#LD = $(FC)
#CF = -fopenmp -lpthread -lm
###############################################################################
#EXEC_F = mainfile
#
#OBJS =\
#library.o solution_lib.o params.o main.o\
###############################################################################

###############################################################################
#EXEC_F: $(OBJS) 
#	$(LD) -o $(EXEC_F) $(CF)  $(OBJS) $(SB_LIB) $(FLIBS)

# $(LDFLAGS)
###############################################################################
#main.o: ../SOFTWARE/libSparseBLAS_GNU.a\
#	library.o solution_lib.o params.o\
###############################################################################


#library.o : library.f95
#	$(FC) $(CF)  $(SB_INCL) -c $*.f95

#solution_lib.o : solution_lib.f95
#	$(FC) $(CF)  $(SB_INCL) -c $*.f95

#params.o : params.f95
#	$(FC) $(CF)  $(SB_INCL) -c $*.f95

#main.o : MAIN.f95 params.o library.o solution_lib.o
#	$(FC) $(CF)  $(SB_INCL) -c $*.f95

#library.o : library.f95
#		$(GFORTRAN)$ -c $< -o $@

#solution_lib.o : solution_lib.f95
#		$(GFORTRAN)$ -c $< -o $@

#param.o : params.f95
#		$(GFORTRAN)$ -c $< -o $@

#main.o : MAIN.f95 param.o library.o solution_lib.o
#		$(GFORTRAN)$  -c $< -o $@ 


#############################
# executable
#mainfile.exec: main.o library.o solution_lib.o param.o
#		$(GFORTRAN)$  $^ -o $@
#############################

clean :
		rm -f *.exec *.o *.mod *.out *.f95.*
