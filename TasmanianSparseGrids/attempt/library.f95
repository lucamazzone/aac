module library
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine creategrid(lmax,bmax,stepl,stepb,vecsize,lgrid,bgrid)
implicit none

!input/output declarations
integer, intent(in) :: vecsize
double precision, intent(in) :: bmax,lmax,stepl,stepb
double precision, intent(out) :: lgrid(vecsize,vecsize), bgrid(vecsize,vecsize)

!other declarations
double precision :: lvector(vecsize),bvector(vecsize)
integer :: ii


    call linspace(lvector,stepl,lmax,vecsize)
    !vecsize = the dimension of the output vector
    !lvector = the nx1 output vector, w/ equally spaced points b/w  x and y
    !stepl  = the minimum of the linear grid
    !lmax  = the maximum of the linear grid

    do ii=1,vecsize
       lgrid(:,ii)  = lvector
    end do 

    call linspace(bvector,stepb,bmax,vecsize)
    !vecsize  = the dimension of the output vector
    !bvector = the nx1 output vector, w/ equally spaced points b/w  x and y
    !stepb  = the minimum of the linear grid
    !bmax = the maximum of the linear grid

    do ii=1,vecsize
       bgrid(:,ii)  = bvector
    end do

    bgrid = transpose(bgrid)
    
  end subroutine creategrid



subroutine qsimpweightsnodes(a,b,n,weights,nodes)
implicit none

!a = lower limit of integration
!b = upper limit of integration
!n = number of intervals, must be even.  There will be n+1 weights/nodes
!weights = n+1 x 1 vector of integration weights
!nodes = n+1 x 1 vector integration nodes

!input/output declarations
integer, intent(in)  :: n
double precision, intent(in)  :: a,b
double precision, intent(out)  :: weights(n+1),nodes(n+1)

!other declarations
integer :: ct

weights(:) = 1.0

do ct=1,n+1
    nodes(ct) = a + ((b-a)/dble(n)) * (dble(ct)-1.0)
end do !ct

do ct=2,n
    if (mod(ct,2)==0) then 
        weights(ct) = 4.0
    else 
        weights(ct) = 2.0
    end if
end do !ct

weights = weights * ( (b-a) / ( 3.0 * dble(n) ) )

end subroutine qsimpweightsnodes

end module library
