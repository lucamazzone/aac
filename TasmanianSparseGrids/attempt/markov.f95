program markov
integer, parameter :: start = 1
integer, parameter ::  samples = 36
integer, parameter ::  t = 30
double precision, allocatable ::  chains(:,:)
integer, parameter :: I = 8585739
double precision :: pr_mat_z(5,5),z0(5)
integer, parameter :: snum = 2
double precision, parameter :: rhosigma =0.75 , phi = 0.07, nstdevz =1.0, musigma = 0.1,rhoz = 0.7
double precision :: logS(2), SS(2),Sprob(2,2)
integer, parameter :: Zsize =10 
double precision :: Zprob(Zsize,Zsize), logz(Zsize),zeta(Zsize,snum)
integer :: aggregate



call tauchen(snum,rhosigma,phi,nstdevz,Sprob,logS)
print*, Sprob
SS = exp(log(musigma) + logS(:))
do aggregate = 1,2
   call tauchen(Zsize,rhoz,SS(aggregate),nstdevz,Zprob,logz)
   zeta(:,aggregate) = exp(logz)
end do


print*, 'zeta', zeta(:,1)
print*, 'zeta', zeta(:,2)


allocate(chains(t-1,samples))

Sprob(1,1) =   0.87158036240791925     
Sprob(2,1) =   0.12841963759208072      
Sprob(2,2) =   0.12841963759208075      
Sprob(1,2) =   0.87158036240791925

!call srand(I)

call markovchain(chains,Sprob,start,samples,t)

open(unit=10006,file='chainsmat.txt',ACTION="write",STATUS="new")
do kkk=1,t-1
write(10006,'(*(F14.7))')(real(chains(kkk,jjj) ),jjj=1,samples)
end do
close(10006)

contains

	subroutine markovchain(chains,Sprob,start,samples,t)

        integer, intent(in) :: start,samples,t
        double precision, intent(out) :: chains(t-1,samples)
        double precision, intent(in)  :: Sprob(2,2)
        double precision :: X(t-1,1), triangular(2,2),cum(2,2),ppi(2+1,1),ppiz(2,1)
        double precision :: s(2,1),chain(t-1,1)
        integer :: state(t-1,2),V(2,1),iii,jjjj,kkk
	integer,parameter :: seed = 6477876

       	call srand(seed)
	print*, rand(),rand()        
	do iii = 1,samples
        
        do jjj = 1,t-1
        X(jjj,1) = rand()
        end do
        
        s = 0
        s(start,1) = 1
        triangular = 1.0
        triangular(2,1) = 0
        !ppi(1,1) = 0.0
        
        V(1,1) = 1
        V(2,1) = 2
        
        cum = matmul(Sprob,triangular)
!        print*, 'cum', cum


        state(1,:) = s(:,1)
        do kkk = 1,t-1
        ppiz = matmul(Sprob,s)   ! cum
!	print*, ppiz, 'ppiz'
        ppi(:,1) = [0.d0,ppiz(1,1),1.d0]
	!print*, 'ppi', ppi
        
        	!print*, 'x val'  , X(kkk,1)  !, 'big pp', ppi(jjj+1,1), 'small pp', ppi(jjj,1)
                !if ((X(kkk,1) .LT. ppi(jjj+1,1)) .AND. (X(kkk,1) .GT. ppi(jjj,1) )) then
                !print*, X(kkk,1), ppi(2,1),ppi(3,1)
                if (X(kkk,1) .GT. ppi(2,1) .AND. X(kkk,1) .LT. ppi(3,1) ) then
                !print*, 'x val', X(kkk,1), 'big  pp', ppi(jjj+1,1), 'small pp', ppi(jjj,1)
                !state(kkk,:) = s(:,1)
                s(2,1) = 1
                s(1,1) = 0
                else
                s(1,1) = 1
                s(2,1) = 0
                end if

            state(kkk,:) = s(:,1)
		chains(kkk,iii) = dot_product(state(kkk,:),V(:,1))
		!print*, chains(kkk,iii)
        end do
	!print*, 'chain', chains(:,iii)
	!chain(1,1) = start
        !chain = matmul(state(1,:),V(:,1))
	!chain(1,1)    =	DOT_PRODUCT(state(1,:), V(:,1))
	!print*, chain(1,1)
	!print*, chain
        !chains(:,iii) = chain(:,1)
        end do
        

	!print*, chains
	end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function normcdf(x,mu,sigma)
implicit none

!input/output declarations
!x: the input value for the normal cdf
!mu: the mean of the normal distribution
!sigma: the standard deviation of the normal distribution
double precision :: x,mu,sigma

!other declarations
double precision :: z

!standardized value ~ N(0,1)
z = (x-mu)/sigma

normcdf =  0.5 * erfcc( ( -1.0 * z ) / sqrt(2.0) )

end function normcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function erfcc(x)
implicit none

!input/output declarations
double precision :: x

!other declarations
double precision :: t,z

z = abs(x)
t = 1.0 / (1.0 + 0.5*z)

erfcc = t * exp(-z * z-1.26551223+t*(1.00002368+t*(0.37409196+&
    t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+&
    t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))

if (x.lt.0.0) erfcc = 2.0 - erfcc
end function erfcc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linspace(z,x,y,n)
    implicit none

    !n = the dimension of the output vector
    !z = the n x 1 output vector, with equally spaced points between x and y
    !x = the minimum of the linear grid
    !y = the maximum of the linear grid


    !input/output declarations
    integer :: n
    double precision :: z(n),x,y

    !local declarations
    integer :: i
    double precision :: d

    d = (y-x)/dble(n-1)
    z(1) = x

    do i = 2,n-1
        z(i) = z(i-1) + d
    end do

    z(n) = y

    return

end subroutine linspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)
implicit none

!input/output declarations
integer, intent(in) :: znum
double precision, intent(in) :: rhoz,sigmaz,nstdevz
double precision, intent(out) :: pr_mat_z(znum,znum),z0(znum)

!other declarations
integer :: zct,zprimect
double precision :: zmin,zmax,gridinc,stdev

!determine end points of the grid (log space)
stdev =  ((sigmaz**2.0)/(1-rhoz**2.0))**0.5
zmin = - nstdevz * stdev
zmax = nstdevz * stdev

!insert points into the grid (log space)
call linspace(z0,zmin,zmax,znum)
gridinc = z0(2)-z0(1)

!loop over z states
do zct=1,znum

    !insert transition matrix middle rows
    do zprimect=2,(znum-1)
        pr_mat_z(zct,zprimect) = &
            normcdf(z0(zprimect)+gridinc/2.0,rhoz*z0(zct),sigmaz) - &
            normcdf(z0(zprimect)-gridinc/2.0,rhoz*z0(zct),sigmaz)
    end do !zct

    !first interval and last interval take the rest of the weight
    pr_mat_z(zct,1) = normcdf(z0(1)+gridinc/2.0,rhoz*z0(zct),sigmaz)
    pr_mat_z(zct,znum) = 1.0 - normcdf(z0(znum)-gridinc/2.0,rhoz*z0(zct),sigmaz)

end do !zct

!round the transition matrix
do zct=1,znum
    pr_mat_z(zct,:) = pr_mat_z(zct,:)/sum(pr_mat_z(zct,:))
end do !zct

!convert grid back to z-space
!z0 = exp(z0)

end subroutine tauchen



end program 
