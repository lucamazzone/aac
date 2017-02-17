  program MAIN 
  use params
  use library
  use solution_lib
!  use invert_numeric
  use invert_direct


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
 bgrid,zeta_tomorrow,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi,Q)
!! remember this has to be called Zsize times, one for each curr_state

do kkk=1,Zsize
   do jjj=1,vecsize
      do iii=1,vecsize
      politics0(iii,jjj,kkk) = iii + (jjj-1)*vecsize   
      end do
   end do
end do

!!!!!!! experiment for valfun
 Q = 1.0
 qq = beta
 Cons= 0.6
 Nbig = 0.7
 Ybig = 0.72
 !zeta1 = zeta(:,1)
 !tv = 0.0
 
 wage = (Cons**eta)*(Nbig**chi)
 


 ! external DGETRI

! call valuefun(vecsize,Zsize,politics0,Zprob,Nbig,Ybig,Cons,Q,qq,zeta,wage,beta,kappa,&
!           gamma,alpha,bgrid,lgrid,v)

 
!print*, v

   


    

  !    call pardisoinit(pt, mtype, solver, iparm, dparm, error)
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

      subroutine valuefun(vecsize,Zsize,politics,Zprob,Nbig,Ybig,Cons,Q,qq,zeta,wage,beta,kappa,&
           gamma,alpha,bgrid,lgrid,v)
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


!!!!!!!!!!!!!
T = 0.0
EYE = 0.0

do kkk=1,Zsize
   polco = reshape(politics(:,:,kkk),(/vecsize**2/))
   
  do jjj=1,Zsize
      do iii=1,vecsize**2
         T( iii + (kkk-1)*vecsize**2 , polco(iii)+ (jjj-1)*vecsize**2 ) = &
              Zprob(kkk,jjj)
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
      wage*lgrid(iii,jjj) - bgrid(iii,jjj)
 lpol =  mod(politics(iii,jjj,curr_state)-1,vecsize)+1
 bpol =  (politics(iii,jjj,curr_state)+vecsize-1)/vecsize
 ut0(iii,jjj,curr_state) = kappa*x(iii,jjj,curr_state)+qq(lpol,bpol,curr_state)*bgrid(lpol,bpol)

       end do
    end do
end do


ut = reshape(ut0,(/Zsize*vecsize**2/))


! OMEGA Relaxation coefficient = 1.0
!tv = -1.0
!call seidel(3,Zsize*vecsize**2,BUM,ut,1.0_8,tv,res,iter,rc)
!v = reshape(tv,(/vecsize,vecsize,Zsize/))


call inverse(BUM,MUB,Zsize*vecsize**2)



 end subroutine valuefun

    end program MAIN

 
