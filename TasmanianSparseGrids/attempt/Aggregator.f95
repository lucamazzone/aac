program Aggregator

use params
use solution_lib
use library

implicit none

double precision :: intvecmat(snum,Zsize), distribution1(vecinterp,Zsize*vecinterp), distribution2(vecinterp,Zsize*vecinterp)
double precision :: rhomat(momnum*Zsize,snum)
!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


open(10003, file="distribution1.txt") 
read(10003,*) distribution1
close(10003)

open(10004, file="distribution2.txt")
read(10004,*) distribution2
close(10004)

open(10002, file="intvectors.txt")! read in values
read(10002,*) intvecmat
close(10002)

open(10001, file="rhomatrix.txt")
read(10001,*) rhomat
close(10001)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! print*, intvecmat(1,:)
















end program Aggregator