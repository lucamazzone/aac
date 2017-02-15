program MAIN 
  use params
  use library
  use solution_lib
!  USE BLAS_SPARSE

  implicit none
  
!!!!!!!!!!!!!
integer, parameter :: snum = 2
double precision  :: logS(snum), Sprob(snum,snum), SS(snum)
double precision  :: logz(Zsize), Zprob(Zsize,Zsize), zeta(Zsize,snum)
double precision  :: s(Zsize,1),s_alt(Zsize,1)
integer :: zct,aggregate,curr_state
double precision :: lgrid(vecsize,vecsize),bgrid(vecsize,vecsize)
double precision :: lgrid_int(vecinterp,vecinterp),bgrid_int(vecinterp,vecinterp)
double precision :: weights(nsimp+1), nodes(nsimp+1)
double precision :: weights_b(nsimp+1), nodes_b(nsimp+1)
double precision :: qfun(vecsize,vecsize)
double precision :: Ybig,Cons,Cons_1,Nbig_1,Ybig_1,zeta_tomorrow(Zsize,1)

!!!!!!!!!!!!!

integer :: pt(64)
integer :: mtype, solver,error
integer :: iparm(64)
real(kind=8) :: dparm(64)


!!!!!!!!!!!!!

!! construct stochastic process

call tauchen(snum,rhosigma,phi,nstdevz,Sprob,logS)
SS = exp(log(musigma) + logS(:))

do aggregate = 1,2
   call tauchen(Zsize,rhoz,SS(aggregate),nstdevz,Zprob,logz)
   zeta(:,aggregate) = exp(logz)
end do

!! obtain stationary distribution

do zct = 1,Zsize
s(zct,1) = 1.0/Zsize
s_alt(zct,1) = 0.0
end do

do while ( maxval(abs(s-s_alt)) .GT.  0.00001)
   s_alt = s
   s = matmul(transpose(Zprob(:,:)),s_alt(:,:))
end do

!! build grids

call creategrid(lmax,bmax,stepl,stepb,vecsize,lgrid,bgrid)
call creategrid(lmax,bmax,stepl,stepb,vecinterp,lgrid_int,bgrid_int)

!! generate Simpson nodes
   
call qsimpweightsnodes(stepl,lmax,nsimp,weights,nodes)
call qsimpweightsnodes(stepb,bmax,nsimp,weights_b,nodes_b)

call  q_fun(qfun,Ybig,Cons,Cons_1,Nbig_1,Ybig_1,lgrid,&
 bgrid,zeta_tomorrow,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi)


   mtype = 11
    solver = 0
    

      call pardisoinit(pt, mtype, solver, iparm, dparm, error)
  !    IF (error .NE. 0) THEN
  !      IF (error.EQ.-10 ) WRITE(*,*) 'No license file found'
  !      IF (error.EQ.-11 ) WRITE(*,*) 'License is expired'
  !      IF (error.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
  !      STOP
  !    ELSE
  !      WRITE(*,*) '[PARDISO]: License check was successful ... '
  !    END IF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!! SUBROUTINES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains

      subroutine valuefun(vecsize,Zsize,politics,Zprob,outt)

  integer,intent(in) :: vecsize,Zsize
  double precision, intent(in) :: politics(vecsize,vecsize,Zsize),Zprob(Zsize,Zsize) 
  integer, intent(out) :: outt
  double precision:: polco(vecsize*vecsize)
  integer :: G, ISTAT,N,i,A,j,Zsq
  integer, DIMENSION(:), ALLOCATABLE :: INDX,JNDX
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: VAL

  N=vecsize**2
  Zsq=Zsize**2
   ALLOCATE (VAL(N*Zsq),INDX(N*Zsq),JNDX(N*Zsq))

   
   CALL duscr_begin( N*Zsize, N*Zsize, G, ISTAT)

   do j=1,Zsize
    polco = reshape(politics(:,:,j),(/1/))         
    do i = 1,N*Zsize      
       INDX(i+(j-1)*N) = i
       JNDX(i+(j-1)*N) = polco(i)+(j-1)*N
       VAL(i+(j-1)*N) = Zprob((i-1)/j +1,j)
   end do
   end do

   CALL duscr_insert_entries(G, VAL, INDX, JNDX, ISTAT)

    ! ut = [];

 !     for curr_state = 1:5;

!    q = qq(:,:,curr_state);
!    x(:,:,curr_state) = zeta1(curr_state)*Ybig^(1/gamma)*lgrid.^(alpha-alpha/gamma)
!    - Nbig^chi*Ybig^eta*lgrid - bgrid;
!    ut0 = kappa*(x(:,:,curr_state) +
!    q(politics(:,:,curr_state)).*bgrid(politics(:,:,curr_state)));
!    ut = [ut, ut0];
!            end

!    ciao =   reshape(ut,vecsize,vecsize,5);
!    tv = (speye(5*vecsize^2) - Q*beta*trans)\ciao(:);

    !  CALL USPG(A,blas_num_nonzeros,outt)
 
    
 
 end subroutine valuefun

    end program MAIN

 
