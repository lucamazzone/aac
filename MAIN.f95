  program MAIN 
  use params
  use library
  use solution_lib
  use omp_lib


  implicit none
  
!!!!!!!!!!!!!
integer, parameter :: snum = 2
integer :: zct,aggregate,curr_state
integer :: iii,jjj,kkk,rc

double precision  :: logS(snum), Sprob(snum,snum), SS(snum)
double precision  :: logz(Zsize), Zprob(Zsize,Zsize), zeta(Zsize,snum)
double precision  :: s(Zsize,1),s_alt(Zsize,1)

double precision :: lgrid(vecsize,vecsize),bgrid(vecsize,vecsize)
double precision :: lgrid_int(vecinterp,vecinterp),bgrid_int(vecinterp,vecinterp)
double precision :: weights(nsimp+1), nodes(nsimp+1)
double precision :: weights_b(nsimp+1), nodes_b(nsimp+1)
double precision :: qfun(vecsize,vecsize),Q, v(vecsize,vecsize,Zsize),qq(vecsize,vecsize,Zsize)
double precision :: Ybig,Cons,Cons_1,Nbig_1,Ybig_1,zeta_tomorrow(Zsize,1),Nbig,wage

!!

integer :: iter,maxiter,nn=Zsize*vecsize**2,politics(Zsize*vecsize**2)
double precision :: value0(vecsize**2,Zsize),expv0(vecsize**2,Zsize),vvalue(Zsize*vecsize**2),zeta1(Zsize)
double precision :: vvalue0(Zsize*vecsize**2),l_grid(vecsize**2),b_grid(vecsize**2),q_q(vecsize**2)
double precision :: obj(vecsize**2),epsilon





!integer :: polco(Zsize*vecsize**2),ia(1+Zsize*vecsize**2), nonzeros
!integer,dimension(:), allocatable :: ja
!real*8,dimension(:), allocatable :: a

!integer :: pt(64)
!integer :: iparm(64)
!real*8 :: dparm(64)
!integer :: mtype,solver,error,nn,lpol,bpol,phase,msglvl,mnum=1,idum,nrhs,maxfct=1
!real*8 :: x(vecsize,vecsize,Zsize),ut0(vecsize,vecsize,Zsize),vv(Zsize*vecsize**2), ut(Zsize*vecsize**2),zeta1(Zsize)
!real*8 :: ddum



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


!do kkk=1,Zsize
!   do jjj=1,vecsize
!      do iii=1,vecsize
!      politics0(iii,jjj,kkk) = (iii + (jjj-1)*vecsize)   
!      end do
!   end do
!end do


!!!!!!! experiments for valfun

 Cons= 0.6
 Nbig = 0.7
 Ybig = 0.72
 Nbig_1 = 0.66
 Ybig_1 = Ybig
 Cons_1 = 0.5
 zeta1 = zeta(:,1)
 wage = (Cons**eta)*(Nbig**chi)
 
 
 
 
do curr_state = 1,Zsize
call  q_fun(qfun,Ybig,Cons,Cons_1,Nbig_1,Ybig_1,lgrid,&
&  bgrid,zeta1,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi,Q)
qq(:,:,curr_state) = qfun(:,:)
end do


value0 = 1.0
maxiter = 10000
l_grid = reshape(lgrid,(/vecsize**2/))
b_grid = reshape(bgrid,(/vecsize**2/))
epsilon = 20.0


do iter=1,maxiter


do curr_state = 1,Zsize
expv0(:,curr_state) = matmul(Zprob(curr_state,:),transpose(value0(:,:)))
end do
    if (epsilon<0.001) then
	do iii=1,nn
	    kkk = (iii+vecsize**2-1)/(vecsize**2)
	    q_q = reshape(qq(:,:,kkk),(/vecsize**2/))
	    jjj = mod(iii-1,vecsize**2)+1
obj = kappa*( zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-wage*l_grid(jjj)-b_grid(jjj) + b_grid*q_q)+ &
	    & beta*(1-kappa)*Q*expv0(:,kkk)
	    politics(iii) =  maxloc(obj,1)
	end do
	exit
    else
	do iii=1,nn
	    kkk = (iii+vecsize**2-1)/(vecsize**2)
	    q_q = reshape(qq(:,:,kkk),(/vecsize**2/))
	    jjj = mod(iii-1,vecsize**2)+1
obj = kappa*( zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-wage*l_grid(jjj)-b_grid(jjj) + b_grid*q_q)+ &
	    & beta*(1-kappa)*Q*expv0(:,kkk)
	    vvalue(iii)  = maxval(obj)
	end do
    end if
	vvalue0 = reshape(value0,(/nn/))
	epsilon = norm2(vvalue0-vvalue)
	value0 = reshape(vvalue,(/vecsize**2,Zsize/))
end do


do curr_state = 1,3
!$OMP PARALLEL PRIVATE(kkk,q_q,jjj,obj) SHARED(qq,zeta1,Ybig,l_grid,b_grid)

!$OMP DO 
do iii=1,nn
    kkk = (iii+vecsize**2-1)/(vecsize**2)
    q_q = reshape(qq(:,:,kkk),(/vecsize**2/))
    jjj = mod(iii-1,vecsize**2)+1
!obj = kappa*(zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-wage*l_grid(jjj)-b_grid(jjj)+b_grid*q_q) +&
!& 	 beta*(1-kappa)*Q*expv0(:,kkk)
    obj = objectif(jjj,kkk,kappa,gamma,alpha,beta,zeta1,Ybig,wage,Q,q_q,l_grid,b_grid,expv0,vecsize,Zsize)
    vvalue(iii) = maxval(obj)
end do
!$OMP END DO

!$OMP END PARALLEL

print*, curr_state
if (curr_state>1)then
print*, 'exiting!'
exit
end if
end do 

print*, 'parallel V', vvalue(1:10)
print*, 'sequential V', vvalue0(1:10)

!do kkk = 1,Zsize
!    do jjj = 1,vecsize
!	do iii = 1,vecsize
!	    x(iii,jjj,kkk) = zeta1(kkk)*(Ybig**(1/gamma))*lgrid(iii,jjj)**(alpha-alpha/gamma)- &
!	    & wage*lgrid(iii,jjj) - bgrid(iii,jjj)
!	    lpol = mod(politics(iii,jjj,kkk)-1,vecsize)+1
!	    bpol = (politics(iii,jjj,kkk)+vecsize-1)/vecsize
!	    ut0(iii,jjj,kkk) = kappa*x(iii,jjj,kkk)+qq(lpol,bpol,kkk)*bgrid(lpol,bpol)
!	end do
!    end do
!end do

!ut = reshape(ut0,(/Zsize*vecsize**2/))






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!! SUBROUTINES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains
    
    
     function objectif(jjj,kkk,kappa,gamma,alpha,beta,zeta1,Ybig,wage,Q,q_q,l_grid,b_grid,expv0,vecsize,Zsize)
     integer, intent(in) :: kkk,jjj,vecsize,Zsize
     double precision, intent(in) :: kappa, gamma, alpha, beta, Q, Ybig, wage,expv0(vecsize**2,Zsize)
     double precision, intent(in) :: l_grid(vecsize**2),b_grid(vecsize**2),q_q(vecsize**2),zeta1(Zsize)
     double precision :: objectif(vecsize**2)
     
     objectif = kappa*(zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-wage*l_grid(jjj)-&
     & b_grid(jjj) + b_grid*q_q) + beta*(1-kappa)*Q*expv0(:,kkk)
     
     end function
     
     

      subroutine valuefun(vecsize,Zsize,politics,Zprob,Nbig,Ybig,Cons,Q,qq,zeta,wage,beta,kappa,&
        &   gamma,alpha,bgrid,lgrid,v)
integer, intent(in) :: vecsize, Zsize, politics(vecsize,vecsize,Zsize)
double precision, intent(in) :: Nbig, Ybig, Cons,Q,zeta(Zsize),wage,Zprob(Zsize,Zsize)
double precision, intent(in) :: beta,kappa,gamma,alpha,lgrid(vecsize,vecsize),bgrid(vecsize,vecsize)
double precision, intent(out) ::  v(vecsize,vecsize,Zsize)


!! other declarations
double precision :: T(Zsize*vecsize**2,Zsize*vecsize**2), EYE(Zsize*vecsize**2,Zsize*vecsize**2)
double precision :: BUM(Zsize*vecsize**2,Zsize*vecsize**2),MUB(Zsize*vecsize**2,Zsize*vecsize**2)
double precision :: x(vecsize,vecsize,Zsize), ut0(vecsize,vecsize,Zsize), ut(Zsize*vecsize**2)
double precision :: zeta1(Zsize), tv(Zsize*vecsize**2),res(Zsize*vecsize**2),qq(vecsize,vecsize,Zsize)
integer::  polco(vecsize**2), iii,jjj,kkk,lpol,bpol,iter

! pardiso declarations
integer :: pt(64)
integer :: iparm(64)
integer :: dparm(64)
integer :: solver, mtype, error


!!!!!!!!!!!!!
solver = 0
mtype = 11

T = 0.0
EYE = 0.0

do kkk=1,Zsize
   polco = reshape(politics(:,:,kkk),(/vecsize**2/))
   
  do jjj=1,Zsize
      do iii=1,vecsize**2
         T( iii + (kkk-1)*vecsize**2 , polco(iii)+ (jjj-1)*vecsize**2 ) = &
        &      Zprob(kkk,jjj)
      end do
  end do
end do


 do iii=1,vecsize**2
    EYE(iii,iii)= 1.0
 end do

 BUM = EYE - (1-kappa)*Q*beta*T

 
 do curr_state = 1,Zsize
    do jjj = 1,vecsize       
       do iii = 1,vecsize
    
 x(iii,jjj,curr_state) = zeta(curr_state)*(Ybig**(1/gamma))*lgrid(iii,jjj)**(alpha-alpha/gamma)-&
    &  wage*lgrid(iii,jjj) - bgrid(iii,jjj)
 lpol =  mod(politics(iii,jjj,curr_state)-1,vecsize)+1
 bpol =  (politics(iii,jjj,curr_state)+vecsize-1)/vecsize
 ut0(iii,jjj,curr_state) = kappa*x(iii,jjj,curr_state)+qq(lpol,bpol,curr_state)*bgrid(lpol,bpol)

       end do
    end do
end do


ut = reshape(ut0,(/Zsize*vecsize**2/))

!! initialize pardiso
call pardisoinit(pt,mtype,solver,iparm,dparm,error)

iparm(3) = 1

!! check 
! call pardiso_chkmatrix(mtype,n,a,ia,ja,error)

!tv = -1.0
!call seidel(3,Zsize*vecsize**2,BUM,ut,1.0_8,tv,res,iter,rc)
!v = reshape(tv,(/vecsize,vecsize,Zsize/))


!call inverse(BUM,MUB,Zsize*vecsize**2)



 end subroutine valuefun

    end program MAIN

 
