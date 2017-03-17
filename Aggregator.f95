program Aggregator

use parameters
use solution_lib
use library

implicit none

double precision :: intvecmat(2,Zsize)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(10002, file="intvecmat.txt")! read in values
read(10002,*) intvecmat
close(10002)
















end program Aggregator