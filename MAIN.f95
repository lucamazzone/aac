  program MAIN 
  use params
  use library
  use solution_lib


  implicit none
  
!!!!!!!!!!!!!
integer, parameter :: snum = 2
integer :: zct,aggregate,curr_state
integer :: iii,jjj,kkk,rc
integer :: politics(vecsize,vecsize,Zsize),politics0(vecsize,vecsize,Zsize)

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

integer :: polco(Zsize*vecsize**2),ia(1+Zsize*vecsize**2), nonzeros
integer,dimension(:), allocatable :: ja
real*8,dimension(:), allocatable :: a

integer :: pt(64)
integer :: iparm(64)
real*8 :: dparm(64)
integer :: mtype,solver,error,nn,lpol,bpol,phase,msglvl,mnum=1,idum,nrhs,maxfct=1
real*8 :: x(vecsize,vecsize,Zsize),ut0(vecsize,vecsize,Zsize),vv(Zsize*vecsize**2), ut(Zsize*vecsize**2),zeta1(Zsize)
real*8 :: ddum



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

do curr_state = 1,Zsize
call  q_fun(qfun,Ybig,Cons,Cons_1,Nbig_1,Ybig_1,lgrid,&
&  bgrid,zeta_tomorrow,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi,Q)
qq(:,:,curr_state) = qq(:,:)
end do

do kkk=1,Zsize
   do jjj=1,vecsize
      do iii=1,vecsize
      politics0(iii,jjj,kkk) = (iii + (jjj-1)*vecsize)   
      end do
   end do
end do


!!!!!!! experiments for valfun
Q = 1.0

polco =reshape(politics0(:,:,:),(/Zsize*vecsize**2/))


ia(1) = 1

do jjj = 1,Zsize
    do iii = 1,vecsize**2
	if (polco(iii+ (jjj-1)*vecsize**2).EQ.iii)  then
	ia( iii+1 + (jjj-1)*vecsize**2  ) = ia(iii + (jjj-1)*vecsize**2 ) + 1 + Zsize - 1
	else 
	ia( iii+1 + (jjj-1)*vecsize**2  ) = ia(iii + (jjj-1)*vecsize**2 ) + 1 + Zsize
	end if
    end do
end do

nonzeros = ia(1+Zsize*vecsize**2)-1

allocate( ja(nonzeros))
allocate( a(nonzeros))



do iii=1,Zsize*vecsize**2
    if      ( mod(iii-1,vecsize**2)+1  < polco(iii) ) then
	ja(ia(iii)) = iii + (jjj-ia(iii))*vecsize**2
	do jjj = ia(iii)+1,ia(iii+1)-1
	    ja(jjj) = polco(iii) + (jjj-ia(iii))*vecsize**2
	end do
    else if ( mod(iii-1,vecsize**2)+1 > polco(iii) ) then
	do jjj = ia(iii),ia(iii+1)-2
	    ja(jjj) = polco(iii) + (jjj-ia(iii))*vecsize**2
	end do
	ja(ia(iii+1)-1) = iii + (jjj-ia(iii))*vecsize**2
    else
	do jjj = ia(iii),ia(iii+1)-1
	    ja(jjj) = polco(iii) + (jjj-ia(iii))*vecsize**2
	end do
    end if
end do



do iii=1,Zsize*vecsize**2
    if      (mod(iii-1,vecsize**2)+1 < polco(iii)) then
	a(ia(iii)) = 1.0
	do jjj = ia(iii)+1,ia(iii+1)-1
	    a(jjj) = -beta*Q*(1-kappa)*Zprob ( (iii+vecsize**2-1)/(vecsize**2),jjj-ia(iii)-1)
	end do
    else if (mod(iii-1,vecsize**2)+1 > polco(iii)) then
	do jjj = ia(iii),ia(iii+1)-2
	    a(jjj) = -beta*Q*(1-kappa)*Zprob( (iii+vecsize**2-1)/(vecsize**2),jjj-ia(iii))
	end do
	a(ia(iii+1)-1) = 1.0
    else
	do jjj= ia(iii),ia(iii+1)-1
	    a(jjj) = -beta*Q*(1-kappa)*Zprob( (iii+vecsize**2-1)/(vecsize**2),jjj-ia(iii))
	end do
	a(ia(iii)+(iii+vecsize**2-1)/(vecsize**2)-1) = &
	& 1-beta*Q*(1-kappa)*Zprob((iii+vecsize**2-1)/(vecsize**2),(iii+vecsize**2-1)/(vecsize**2))
    end if
end do



solver = 0
mtype = 11
nn = Zsize*vecsize**2

call pardisoinit(pt,mtype,solver,iparm,dparm,error)

iparm(3)=1

call pardiso_chkmatrix(mtype,nn,a,ia,ja,error)



 qq = beta
 Cons= 0.6
 Nbig = 0.7
 Ybig = 0.72
 zeta1 = zeta(:,1)
 !tv = 0.0
 wage = (Cons**eta)*(Nbig**chi)

do kkk = 1,Zsize
    do jjj = 1,vecsize
	do iii = 1,vecsize
	    x(iii,jjj,kkk) = zeta1(kkk)*(Ybig**(1/gamma))*lgrid(iii,jjj)**(alpha-alpha/gamma)- &
	    & wage*lgrid(iii,jjj) - bgrid(iii,jjj)
	    lpol = mod(politics(iii,jjj,kkk)-1,vecsize)+1
	    bpol = (politics(iii,jjj,kkk)+vecsize-1)/vecsize
	    ut0(iii,jjj,kkk) = kappa*x(iii,jjj,kkk)+qq(lpol,bpol,kkk)*bgrid(lpol,bpol)
	end do
    end do
end do

ut = reshape(ut0,(/Zsize*vecsize**2/))


!call pardiso_chkvec(nn,1,ut,error)

 
!call  pardiso_printstats(mtype,nn,a,ia,ja,1,ut,error)


phase = 11  ! only reordering and symbolic factorization
msglvl = 1  ! with statistical information

!call pardiso(pt,maxfct,mnum,mtype,phase,nn,a,ia,ja,idum,nrhs,iparm,msglvl,ddum,ddum,error,dparm)

phase = 22  ! only factorization

!call pardiso(pt,maxfct,mnum,mtype,phase,nn,a,ia,ja,idum,nrhs,iparm,msglvl,ddum,ddum,error,dparm)

phase = 33 ! only solve
iparm(8) = 1 ! max number of iterative refinement steps
iparm(2) = 0 ! try

!call pardiso(pt,maxfct,mnum,mtype,phase,nn,a,ia,ja,idum,nrhs,iparm,msglvl,ut,vv,error,dparm)

!phase = 33
!iparm(8) = 1
!iparm(12) = 1


!call pardiso(pt,maxfct,mnum,mtype,phase,nn,a,ia,ja,idum,nrhs,iparm,msglvl,ut,vv,error,dparm)

phase = -1
!call pardiso(pt,maxfct,mnum,mtype,phase,nn,ddum,idum,idum,idum,nrhs,iparm,msglvl,ddum,ddum,error,dparm)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!! SUBROUTINES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains

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

 
