module solution_lib
implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 subroutine  q_fun(qfun,Ybig,Cons,Cons_1,Nbig_1,Ybig_1,lgrid,&
    bgrid,zeta_tomorrow,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi,Q)

! inputs and outputs
integer, intent(in) :: vecsize, curr_state,Zsize
double precision, intent(in) :: Ybig, Cons, Cons_1, Nbig_1, Ybig_1, zeta_tomorrow(Zsize,1), &
     lgrid(vecsize,vecsize), bgrid(vecsize,vecsize), Zprob(Zsize,Zsize), &
     alpha,beta,gamma,eta,chi
double precision, intent(out) :: qfun(vecsize,vecsize),Q
! other declarations
double precision :: bmax, p(vecsize,vecsize,Zsize), profit(vecsize,vecsize,Zsize),&
     repayment(vecsize,vecsize,Zsize),wage
integer :: ii,jj,kk

bmax = bgrid(vecsize,vecsize)
Q =(Cons/Cons_1)**(eta)
wage = Nbig_1**chi * Ybig_1**eta



do ii = 1,Zsize

   p(:,:,ii) = zeta_tomorrow(ii,1)*(Ybig_1/lgrid**alpha)**(1/gamma)
   do jj = 1,vecsize
      do kk = 1,vecsize
         profit(kk,jj,ii) = p(kk,jj,ii)*(lgrid(kk,jj)**alpha) - wage*lgrid(kk,jj)
      end do
   end do
repayment(:,:,ii) = (profit(:,:,ii)-bgrid)/bmax;

end do



do ii = 1,Zsize
   do jj = 1,vecsize
      do kk = 1,vecsize
         
         if (repayment(kk,jj,ii) .LT. 0.0) then
         repayment(kk,jj,ii) = 0.0
         elseif (repayment(kk,jj,ii) .GT. 1.0) then
            repayment(kk,jj,ii) = 1.0
         end if   
      end do
   end do
end do


qfun = repayment(:,:,1)*0.0

do ii = 1,Zsize;  ! 10  = # of discretized idyiosinc shocks
qfun = qfun + beta*Q*Zprob(ii,curr_state)*repayment(:,:,ii);
end do

!qprice(qprice>1) = 1;
!qprice(qprice<0) = 0;


do jj = 1,vecsize
   do kk = 1,vecsize         
         if (qfun(kk,jj) .LT. 0.0) then
         qfun(kk,jj) = 0.0
         elseif (qfun(kk,jj) .GT. 1.0) then
         qfun(kk,jj) = 1.0
         end if   
   end do
end do

    
end subroutine  q_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module solution_lib



