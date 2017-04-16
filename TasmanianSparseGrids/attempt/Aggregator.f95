module Aggregator

! change end program into end module at the end
use params
use solution_lib
use library

!double precision :: points(3)
!double precision :: pred(3)
!integer :: iii,jjj,kkk,rc,aggregate

!points(1) = 0.72
!points(2) = 0.68
!points(3) = 0.15
!aggregate = 1

!call mapping_inverse(points,aggregate,pred)




contains
subroutine mapping_inverse(points,aggregate,pred)

implicit none

!double precision, intent(in) :: points(2),mzero,Tol
!integer, intent(in) :: aggregate
!double precision, intent(out) :: vals(3)
double precision, intent(in) :: points(3)
double precision, intent(out) :: pred(3)
integer, intent(in) :: aggregate
double precision :: vals(3),threshold,mzero

!! other declarations
double precision :: intvecmat(snum,Zsize), distribution1(vecinterp,Zsize*vecinterp), distribution2(vecinterp,Zsize*vecinterp)
double precision :: rhomat(momnum*Zsize,snum), momentsmat(Zsize*momnum,snum),Tol
!!!
integer :: zct,curr_state,loop,agg
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

integer :: iter,maxiter,nn=Zsize*vecsize**2,politics(Zsize*vecsize**2,1),policcc(1)
double precision :: value0(vecsize**2,Zsize),expv0(vecsize**2,Zsize),vvalue(Zsize*vecsize**2),zeta1(Zsize)
double precision :: vvalue0(Zsize*vecsize**2),l_grid(vecsize**2),b_grid(vecsize**2),q_q(vecsize**2)
double precision :: obj(vecsize**2),epsilon,value(vecsize,vecsize,Zsize),cums(Zsize)

!!

double precision :: labpol(Zsize*vecsize**2), debpol(Zsize*vecsize**2),lab_pol(vecsize,vecsize,Zsize),deb_pol(vecsize,vecsize,Zsize)
double precision :: labpol_int(vecinterp,vecinterp,Zsize),debpol_int(vecinterp,vecinterp,Zsize),zi(vecinterp)
double precision :: polprimewgt1(vecinterp,vecinterp,Zsize),polprimewgt2(vecinterp,vecinterp,Zsize)
double precision :: v_int(vecinterp,vecinterp,Zsize)
integer ::  polprimeind1(vecinterp,vecinterp,Zsize),polprimeind2(vecinterp,vecinterp,Zsize)

!!

double precision:: V_e(vecinterp,Zsize),l_y(Zsize),coeffs(2),result,marginals,size_active,implied_consumption,C_high,C_low
integer :: entering(Zsize),v_entry(Zsize),polprimeind3(vecinterp,Zsize)
double precision :: labpol_ent(vecinterp,Zsize),labentry(Zsize),polprimewgt3(vecinterp,Zsize),Nagg,Y_agg,Bagg,err_cons
double precision :: labdist(vecinterp,Zsize), dist(vecinterp,vecinterp,Zsize),nprimesimp(nsimp+1,Zsize),nprime(vecinterp,Zsize)
double precision :: N_1,Y_1,epsiloun,C_pred,intvec(Zsize),Pint,distr(vecinterp,vecinterp*Zsize)
double precision :: momstoremat(Zsize,momnum),rhomatrix(Zsize,momnum),intvector(Zsize)

!!

double precision :: Nref, Nshift, N_prime, Y_prime, zval, wgt(nsimp+1,Zsize), F_k, prodentry(Zsize),nprimeval
integer :: kct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!open(10003, file="distribution1.txt") 
!read(10003,*) distribution1
!close(10003)

!open(10004, file="distribution2.txt")
!read(10004,*) distribution2
!close(10004)

open(10002, file="intvectors.txt")! read in values
read(10002,*) intvecmat
close(10002)

open(10001, file="rhomatrix.txt")
read(10001,*) rhomat
close(10001)

open(10001, file="momentsmat.txt")
read(10001,*) momentsmat
close(10001)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!aggregate = 1


intvector = intvecmat(aggregate,:)
momstoremat = reshape(momentsmat(:,aggregate),(/Zsize,momnum/))
rhomatrix = reshape(rhomat(:,aggregate),(/Zsize,momnum/))
!points(1) = 0.66
!points(2) = 0.64
mzero = points(3)   ! 0.15
Tol = 0.01


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call tauchen(snum,rhosigma,phi,nstdevz,Sprob,logS)
SS = exp(log(musigma) + logS(:))
do agg = 1,2
   call tauchen(Zsize,rhoz,SS(agg),nstdevz,Zprob,logz)
   zeta(:,agg) = exp(logz)
end do



! zeta(:,2) = Zsize*(zeta(:,2))/sum(zeta(:,2))
!! obtain stationary distribution

do zct = 1,Zsize
s(zct,1) = 1.0/Zsize
s_alt(zct,1) = 0.0
end do

do while ( maxval(abs(s-s_alt)) .GT.  0.00001)
   s_alt = s
   s = matmul(transpose(Zprob(:,:)),s_alt(:,:))
end do

!! obtain cumulative distribution

cums(1) = s(1,1)
do iii=2,Zsize
cums(iii) = cums(iii-1) + s(iii,1)
end do

!! build grids

call creategrid(lmax,bmax,stepl,stepb,vecsize,lgrid,bgrid)
call creategrid(lmax,bmax,stepl,stepb,vecinterp,lgrid_int,bgrid_int)

!! generate Simpson nodes
   
call qsimpweightsnodes(stepl,lmax,nsimp,weights,nodes)
call qsimpweightsnodes(stepb,bmax,nsimp,weights_b,nodes_b)

!!!

threshold = Tol/5
epsiloun = 1.0

Nbig = points(1)
Ybig = points(2)
pred(1) = 0.7
pred(2) = 0.7
pred(3) = Ybig
zeta1 = zeta(:,aggregate)

print*, zeta1

loop = 0

do while(epsiloun > threshold) !! Loop over expected future aggregates

    C_low = 0.4
    C_high = 1.0
    C_pred = pred(3)
    N_1 = pred(1)
    Y_1 = pred(2)
    Cons_1 = pred(3)*Y_1/Ybig
    print*, 'forecast for N',  N_1

    loop=loop+1
    print*, 'loop is' , loop

    if (loop .GT. 15) then
      threshold = Tol/2
    else if (loop .GT. 30) then
      threshold = Tol
    else if (loop .GT. 50) then
      threshold = pred(2)/10
    end if
    
    do while( abs(C_high-C_low) .GT. 0.01 )  !! Golden search for market eq
    
      Cons=   0.5*C_low + 0.5*C_high
      wage = (Cons**eta)*(Nbig**chi)
      do curr_state = 1,Zsize
      call  q_fun(qfun,Ybig,Cons,Cons_1,N_1,Y_1,lgrid,&
      &  bgrid,zeta1,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi,Q)
      qq(:,:,curr_state) = qfun(:,:)
      end do
      print*, 'Q', Q

      
      value0 = 1.0
      maxiter = 10000
      l_grid = reshape(lgrid,(/vecsize**2/))
      b_grid = reshape(bgrid,(/vecsize**2/))
      epsilon = 20.0
      
      do iter=1,maxiter     !! Loop for solution of firm problem
        
    		do curr_state = 1,Zsize
	        expv0(:,curr_state) = matmul(Zprob(curr_state,:),transpose(value0(:,:)))
	        end do
        
    		do iii=1,nn
	        kkk = (iii+vecsize**2-1)/(vecsize**2)
	        q_q = reshape(qq(:,:,kkk),(/vecsize**2/))
	        jjj = mod(iii-1,vecsize**2)+1
	        obj = objectif(jjj,kkk,kappa,gamma,alpha,beta,zeta1,Ybig,wage,Q,q_q,l_grid,b_grid,expv0,vecsize,Zsize)
	        vvalue(iii) = maxval(obj)
	        policcc = maxloc(obj)
	        politics(iii,1) = policcc(1)
	        end do
        
        vvalue0 = reshape(value0,(/nn/))
	epsilon = norm2(vvalue0-vvalue)
		    if (epsilon<0.0001)then
		      exit
		    else
		      value0 = reshape(vvalue,(/vecsize**2,Zsize/))
		    end if
      end do !! end of firm problem loop
      
     
      do iii=1,nn
	      labpol(iii) = lgrid(mod(politics(iii,1)-1,vecsize)+1,1)
	      debpol(iii) = bgrid(1,(politics(iii,1)+vecsize-1)/vecsize)
      end do
      
      lab_pol = reshape(labpol,(/vecsize,vecsize,Zsize/))
      deb_pol = reshape(debpol,(/vecsize,vecsize,Zsize/))
      value = reshape(vvalue,(/vecsize,vecsize,Zsize/))
      
      print*, 'lab_pol', sum(lab_pol(15,:,7))/vecsize
      do kkk = 1,Zsize
	    do jjj = 1,vecinterp
	      call pwl_interp_2d(vecsize,vecsize,lgrid(:,1),bgrid(1,:),lab_pol(:,:,kkk),vecinterp,lgrid_int(:,1),bgrid_int(:,jjj),zi)
	      labpol_int(:,jjj,kkk) = zi
	      call pwl_interp_2d(vecsize,vecsize,lgrid(:,1),bgrid(1,:),deb_pol(:,:,kkk),vecinterp,lgrid_int(:,1),bgrid_int(:,jjj),zi)
	      debpol_int(:,jjj,kkk) = zi
	      call pwl_interp_2d(vecsize,vecsize,lgrid(:,1),bgrid(1,:),value(:,:,kkk), vecinterp,lgrid_int(:,1),bgrid_int(:,jjj),zi)
	      v_int(:,jjj,kkk) = zi
	    end do
      end do
      print*, 'lab_int', sum(labpol_int(15,:,7))/vecinterp
     
!      call convertpolicy3(polprimeind1,polprimewgt1,labpol_int,lgrid_int(:,1))
!      call convertpolicy3(polprimeind2,polprimewgt2,debpol_int,bgrid_int(1,:))
      entering = 0
      
      do kkk=1,Zsize
        V_e(:,kkk) = -kappa*csi + &
        & (1-kappa)*beta*Q*matmul(Zprob(kkk,:),transpose(reshape(v_int(:,1,:),(/vecinterp,Zsize/) )))
        v_entry(kkk) = maxloc(V_e(:,kkk),1)
        l_y(kkk) =  V_e(v_entry(kkk),kkk)
          if (l_y(kkk)>0) then
          entering(kkk) = 1
          end if
      end do
      print*, 'value of (attempted) entrants', l_y(9:10)
      marginals = 0.0
      if (sum(entering(2:Zsize)-entering(1:Zsize-1)) >  0)  then   
      print*, 'entering', entering 
    	call coeff(coeffs,Zsize,2,l_y,v_entry)
	    result = -coeffs(1)/coeffs(2)
	    marginals = marginals_entering(Zsize,result,v_entry,cums)
	    if (isnan(marginals)) then
	    marginals = 0.0
	    end if
      end if
      
      do kkk=1,Zsize
	    labentry(kkk) = lgrid_int(v_entry(kkk),1)*entering(kkk)
      end do
     
      print*, 'labentry', labentry(10),'marginals',marginals

      
      labpol_ent = 1.0
      labpol_ent(1,:) = labentry
      size_active = 1.0-mzero+mzero*dot_product(s(:,1),entering) + marginals
      implied_consumption = Ybig - mzero*csi*dot_product(s(:,1),entering) - marginals*csi

      if (Cons-implied_consumption .LT. 0.0) then
	    C_low = Cons
      elseif (Cons-implied_consumption .GT. 0.0) then
	    C_high = Cons
      end if
    print*, 'implied consumption', implied_consumption
    !print*, 'c_high= ', C_high
    !print*, 'c_low= ', C_low 
    
    end do   !! end of mkt clearing loop
    
call convertpolicy3(polprimeind1,polprimewgt1,labpol_int,lgrid_int(:,1))
call convertpolicy3(polprimeind2,polprimewgt2,debpol_int,bgrid_int(1,:))
    
   do kkk=1,Zsize
   	prodentry(kkk) = labentry(kkk)**(alpha*(gamma-1)/gamma)
   end do
   
   call transform_simp(nprimesimp,nprime,Zprob,polprimeind1,polprimewgt2,lgrid_int,nodes) 
   print*, 'nprimesimp', sum(nprimesimp(:,8))/(nsimp+1)
   
   Nref = dot_product(s(:,1),momstoremat(:,1))
   Nshift = N_1/Nref
   
   N_prime = 0.0
   Y_prime = 0.0
   
  
   
   do kct = 1,nsimp+1  ! ,nsimp+1  !! do
   do zct = 1,Zsize !,Zsize  !! do
   zval = zeta1(zct)
   F_k = Fk(nodes(kct),zct,rhomatrix,momstoremat)
   wgt(kct,zct) = (s(zct,1)*weights(kct)*F_k)/intvector(zct)
   end do
   end do
   wgt = wgt/sum(wgt)
    do kct = 1,nsimp+1
    do zct = 1,Zsize
	zval = zeta1(zct)
	N_prime = N_prime + nprimesimp(kct,zct)*wgt(kct,zct)
	Y_prime  = Y_prime +  wgt(kct,zct)*zval*nprimesimp(kct,zct)**(alpha*(gamma-1)/gamma)
    end do
    end do

   
   N_prime = N_prime*(1-mzero) + mzero*dot_product(s(:,1),labentry)
   Y_prime = ((1-mzero)*Y_prime + mzero*dot_product(s(:,1),prodentry))**(gamma/(gamma-1))
   
   
 epsiloun =( abs(N_prime - pred(1)) + abs(Y_prime - pred(2)) + abs(implied_consumption- pred(3)))/3
print*,'epsilon is', epsiloun
print*, 'error on lab', abs(N_prime-pred(1))
print*, 'error on prod', abs(Y_prime-pred(2))
print*, 'error on cons', abs(implied_consumption-pred(3))
print*, 'N_prime is', N_prime
print*, 'Y_prime is', Y_prime
    pred(1) = N_prime*0.15 + pred(1)*0.85
    pred(2) = Y_prime*0.15 + pred(2)*0.85
    pred(3) = implied_consumption*0.25 + pred(3)*0.75
   
end do  !! end of expectations loop


   
    
      









end subroutine mapping_inverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!! SUBROUTINES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
     function objectif(jjj,kkk,kappa,gamma,alpha,beta,zeta1,Ybig,wage,Q,q_q,l_grid,b_grid,expv0,vecsize,Zsize)
     implicit none
     integer, intent(in) :: kkk,jjj,vecsize,Zsize
     double precision, intent(in) :: kappa, gamma, alpha, beta, Q, Ybig, wage,expv0(vecsize**2,Zsize)
     double precision, intent(in) :: l_grid(vecsize**2),b_grid(vecsize**2),q_q(vecsize**2),zeta1(Zsize)
     double precision :: objectif(vecsize**2)
     integer :: iii
     
     objectif = kappa*(zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-wage*l_grid(jjj)-&
&		 b_grid(jjj) + b_grid*q_q) + beta*(1-kappa)*Q*expv0(:,kkk)
     
     
     do iii=1,vecsize**2
	if (    kappa*(zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-&
&	wage*l_grid(jjj) - b_grid(jjj) + b_grid(iii)*q_q(iii))  <   0) then
	 objectif(iii) = -(100.0 + iii/100.0)
         end if
     end do
     
     end function
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	function marginals_entering(Zsize,result,v_entry,cums)
	implicit none
	integer, intent(in) :: Zsize,v_entry(Zsize)
	double precision, intent(in) :: result,cums(Zsize)
	integer :: low, high
	double precision :: marginals_entering
	
	
	low = 1
	high = Zsize
	
	if (result > v_entry(high)) then
	    marginals_entering = 0
	elseif (result < v_entry(low)) then
	    marginals_entering = 0
	else
	    do while(v_entry(low+1)<result)
		low = low+1
	    end do
	    
	    do while(v_entry(high-1)>result) 
		high = high-1
	    end do
	    
	 marginals_entering = ( cums(high)-cums(low))* (result - v_entry(low))/(v_entry(high)-v_entry(low))
	 
	 end if
	 
	 end function marginals_entering


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine transform_simp(nprimesimp,nprime,Zprob,polprimeind1,polprimewgt2,&
& lgrid_int,nodes)
implicit none
integer, intent(in) :: polprimeind1(vecinterp,vecinterp,Zsize)
double precision, intent(in) :: polprimewgt2(vecinterp,vecinterp,Zsize),nodes(nsimp+1)
double precision, intent(in) :: Zprob(Zsize,Zsize),lgrid_int(vecinterp,vecinterp)
double precision, intent(out) :: nprime(vecinterp,Zsize),nprimesimp(nsimp+1,Zsize)
!! other declarations
integer :: ind1,zprime,kkk,jjj,iii
double precision :: wgt2,nprimesimpa(nsimp+1)


nprime = 0.0

do kkk=1,Zsize
    do jjj=1,vecinterp
	do iii=1,vecinterp
	    ind1 = polprimeind1(iii,jjj,kkk)
	    wgt2 = polprimewgt2(iii,jjj,kkk)
		nprime(iii,kkk) = nprime(iii,kkk)+lgrid_int(ind1,1)*wgt2/vecinterp 
		nprime(iii,kkk) = nprime(iii,kkk) + lgrid_int(ind1,1)*(1-wgt2)/vecinterp
	end do
    end do
end do

do kkk=1,Zsize
    call pwl_value_1d(vecinterp,lgrid_int(:,1),nprime(:,kkk),nsimp+1,nodes,nprimesimpa)
    nprimesimp(:,kkk) = nprimesimpa
end do

end subroutine transform_simp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



double precision function Fk(kval,zct,rhomat,momstoremat)
implicit none
!this function evaluates the density-proportional function
!Fk = exp(rho_1 * (k-m^zct_1)+....+rho_momuse * ((k - m^zct_1)^momuse - m^zct_momuse) )
double precision,intent (in) :: kval,rhomat(Zsize,momnum),momstoremat(Zsize,momnum)
integer :: zct,momct
double precision :: rhovec(momnum),momvec(momnum)

rhovec = rhomat(zct,:)
momvec = momstoremat(zct,:)
Fk = rhovec(1)*(kval - momvec(1))

do momct=2,momnum
    Fk = Fk + rhovec(momct) * ( (kval-momvec(1))**dble(momct) - momvec(momct) )
end do !momct

Fk = exp(Fk)

end function Fk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pwl_value_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_VALUE_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  yi(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do
  
  return
end subroutine pwl_value_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
subroutine hunt(xx,n,x,jlo)
    implicit none
    
    !xx = the n x 1 table of values that you're comparing x to
    !n = the dimension of xx
    !x = the value which you're interested in
    !jlo = (on input) the guess for the integer such that x is in between xx(jlo) and xx(jlo+1)
    !jlo = (on output) the value of the integer such that x is in between xx(jlo) and xx(jlo+1)
    
    !input/output declarations
    integer :: jlo,n
    double precision  :: x,xx(n)
    
    !local declarations
    integer :: inc,jhi,jm
    logical :: ascnd
    
    !determine if table is ascending
    ascnd = xx(n).ge.xx(1)
    
    !in case input guess isn't useful, for robustness
    if (jlo.le.0.or.jlo.gt.n) then 
        jlo = 0
        jhi = n+1
        goto 3
    endif
    
    inc=1 !initialize the hunting increment
    
    !hunt up
    if (x.ge.xx(jlo).eqv.ascnd) then
1       jhi = jlo+inc
        if (jhi.gt.n) then
            jhi = n+1
        else if (x.ge.xx(jhi).eqv.ascnd) then
            jlo = jhi
            inc = inc+inc
            goto 1
        end if    
    !hunt down        
    else
        jhi = jlo
2       jlo = jhi - inc
        if (jlo.lt.1) then
            jlo = 0
        else if (x.lt.xx(jlo).eqv.ascnd) then
            jhi = jlo
            inc = inc+inc
            goto 2
        end if
        
        
        
    endif
    
    !now, hunt is done, begin the bisection phase
3   if (jhi-jlo.eq.1) then
        if (x.eq.xx(n)) jlo=n-1
        if (x.eq.xx(1)) jlo=1
        return
    end if
    jm = (jhi + jlo)/2
    if (x.ge.xx(jm).eqv.ascnd) then
        jlo = jm
    else 
        jhi = jm
    end if
    goto 3
    
end subroutine hunt
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
subroutine convertpolicy3(polprimeind,polprimewgt,policy,grid)
implicit none
double precision, intent(in) :: policy(vecinterp,vecinterp,Zsize),grid(vecinterp)
double precision, intent(out) :: polprimewgt(vecinterp,vecinterp,Zsize)
integer, intent(out) :: polprimeind(vecinterp,vecinterp,Zsize)
integer :: zct,kct,ind,kkk,jjj,iii
double precision :: wgt,primeval
  
do kkk  = 1,Zsize
	do jjj = 1,vecinterp
    		do iii = 1,vecinterp
primeval = policy(jjj,iii,kkk)
ind = 1
call hunt(grid,vecinterp,primeval,ind) 
if (ind<1) then
    wgt = 0.0
    ind = 1
else if (ind > 0 .AND. ind < vecinterp ) then
    wgt = (primeval - grid(ind))/(grid(ind+1)-grid(ind))
else if (ind > vecinterp-1) then
    wgt = 1.0
    ind = vecinterp-1
end if
    polprimeind(jjj,iii,kkk) = ind
    polprimewgt(jjj,iii,kkk) = wgt
        end do
    end do 
end do
end subroutine convertpolicy3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module Aggregator
