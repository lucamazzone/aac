###########################################

FLIBS = ../../wmsc-lucamazzone/lectures/day5/code_day5/pardiso_example/lib/libpardiso500-GNU481-X86-64.so -llapack -lblas -fopenmp -lpthread -lm

GFORTRAN = gfortran -O3 -g 


##############################
all :  mainfile.exec aggregator.exec
#############################


library.o : library.f95
	$(GFORTRAN)$ -c $<  -o $@
	
solution_lib.o : solution_lib.f95
	$(GFORTRAN)$ -c $< -o $@
		
param.o : params.f95
	$(GFORTRAN)$ -c $< -o $@
			
main.o : MAIN.f95 param.o library.o solution_lib.o
	$(GFORTRAN)$  -c $< -fopenmp -o $@ 
	
aggregator.o : Aggregator.f95
	$(GFORTRAN)$  -c $< -o $@ 

				
#############################
# executables
mainfile.exec: main.o library.o solution_lib.o param.o
	$(GFORTRAN)$  $^ -o $@  $(FLIBS)
aggregator.exec: aggregator.o
	$(GFORTRAN)$  $^ -o $@  $(FLIBS)
#############################



clean :
		rm -f *.exec *.o *.mod *.out *.f95.* *.txt
