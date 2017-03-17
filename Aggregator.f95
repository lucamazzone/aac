program Aggregator

!use parameters
!use solution_lib
!use library

implicit none

double precision :: intvecmat(2,10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(10002, file="intvectors.txt")! read in values
read(10002,*) intvecmat
close(10002)


print*, intvecmat(1,:)
















end program Aggregator