module Aggregator
implicit none

  double precision, parameter :: alpha = 0.7 !labor productivity
  double precision, parameter :: beta = 0.985 ! hh discount factor => r = 2,67% 
  double precision, parameter :: eta = 1 ! hh IES
  double precision, parameter :: chi = 0.5  ! Fritsch labor elasticity
  double precision, parameter :: gamma = 7.7 ! elasticity of subs b/w goods
  double precision, parameter :: rhoz = 0.7 ! serial corr of idiosync shocks
  double precision, parameter :: rhosigma = 0.85 ! serial corr of unc shocks
  double precision, parameter :: kappa = 0.4 ! Jensen effect
  double precision, parameter :: phi = 0.13 ! std of unc shocks
  double precision, parameter :: musigma = 0.18 ! mean of unc process
  double precision, parameter :: csi = 0.5 ! entry costs

  double precision, parameter :: nstdevz = 1.0  ! stuff for tauchen

  integer, parameter :: Zsize = 10 ! colsize of idiosync shock matrix
  integer, parameter :: vecsize = 40 ! colsize of lgrid and bgrid
  integer, parameter :: vecinterp = 150 ! size of interpolated version
  integer, parameter :: snum = 2 ! high, low uncertainty
  
  integer, parameter :: nsimp = 100 ! number of Simpson quadnodes (ntb even!)
  integer, parameter :: momnum = 5 ! number of reference moments to be used

  ! parameters for the grid of the firm problem

  double precision, parameter :: lmax = 1.3871 ! maxval of labor grid
  double precision, parameter :: bmax  = 0.5 ! maxval of debt grid
  double precision, parameter :: stepl = lmax/vecsize ! distance b/w gridpoints
  double precision, parameter :: stepb = bmax/vecsize ! distance b/w gridpoints

! change end program into end module at the end
!use params
!use solution_lib
!use library

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

subroutine  q_fun(qfun,Ybig,Cons,Cons_1,Nbig_1,Ybig_1,lgrid,&
&    bgrid,zeta_tomorrow,Zprob,curr_state,vecsize,Zsize,alpha,beta,gamma,eta,chi,Q)

! inputs and outputs
integer, intent(in) :: vecsize, curr_state,Zsize
double precision, intent(in) :: Ybig, Cons, Cons_1, Nbig_1, Ybig_1, zeta_tomorrow(Zsize,1), &
&     lgrid(vecsize,vecsize), bgrid(vecsize,vecsize), Zprob(Zsize,Zsize), &
&     alpha,beta,gamma,eta,chi
double precision, intent(out) :: qfun(vecsize,vecsize),Q
! other declarations
double precision :: bmax, p(vecsize,vecsize,Zsize), profit(vecsize,vecsize,Zsize),&
&     repayment(vecsize,vecsize,Zsize),wage
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

function r8vec_bracket5 ( nd, xd, xi )

!*****************************************************************************80
!
!! R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
!
!  Discussion:
!
!    We assume XD is sorted.
!
!    If XI is contained in the interval [XD(1),XD(N)], then the returned 
!    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
!
!    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
!
!    This code implements a version of binary search which is perhaps more
!    understandable than the usual ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data values.
!
!    Input, real ( kind = 8 ) XD(N), the sorted data.
!
!    Input, real ( kind = 8 ) XD, the query value.
!
!    Output, integer ( kind = 4 ) R8VEC_BRACKET5, the bracket information.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) b
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r8vec_bracket5
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi

  if ( xi < xd(1) .or. xd(nd) < xi ) then

    b = -1

  else

    l = 1
    r = nd

    do while ( l + 1 < r )
      m = ( l + r ) / 2
      if ( xi < xd(m) ) then
        r = m
      else
        l = m
      end if
    end do

    b = l

  end if

  r8vec_bracket5 = b

  return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is intended to be the largest
!    representable real value.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.79769313486231571D+308

  return
end function r8_huge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine pwl_interp_2d( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NXD, NYD, the number of X and Y data values.
!
!    Input, real ( kind = 8 ) XD(NXD), YD(NYD), the sorted X and Y data.
!
!    Input, real ( kind = 8 ) ZD(NXD,NYD), the Z data.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the coordinates of the
!    interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nxd
  integer ( kind = 4 ) nyd

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) det
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dxi
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  real ( kind = 8 ) dyi
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!  real ( kind = 8 ) r8_huge
!  integer ( kind = 4 ) r8vec_bracket5
  real ( kind = 8 ) xd(nxd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nyd)
  real ( kind = 8 ) yi(ni)
  real ( kind = 8 ) zd(nxd,nyd)
  real ( kind = 8 ) zi(ni)

  do k = 1, ni

    i = r8vec_bracket5 ( nxd, xd, xi(k) )
    if ( i == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if

    j = r8vec_bracket5 ( nyd, yd, yi(k) )
    if ( j == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if

    if ( yi(k) < yd(j+1) &
    &  + ( yd(j) - yd(j+1) ) * ( xi(i) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

      dxa = xd(i+1) - xd(i)
      dya = yd(j)   - yd(j)

      dxb = xd(i)   - xd(i)
      dyb = yd(j+1) - yd(j)

      dxi = xi(k)   - xd(i)
      dyi = yi(k)   - yd(j)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

    else

      dxa = xd(i)   - xd(i+1)
      dya = yd(j+1) - yd(j+1)

      dxb = xd(i+1) - xd(i+1)
      dyb = yd(j)   - yd(j+1)

      dxi = xi(k)   - xd(i+1)
      dyi = yi(k)   - yd(j+1)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)

    end if

  end do

  return

end subroutine pwl_interp_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coeff(coeffs,n,k,l_y,x_vec)
implicit none 
integer, intent(in) :: n,k,x_vec(n)
double precision,intent(in) :: l_y(n,1)
double precision ::  x_mat(n,k), x_prod(k,k), x_inv(k,k),predict(k,1)
double precision, intent(out) :: coeffs(k,1)
integer :: I
double precision :: linear(n)

linear = (/ (I, I = 1, n) /)
do I = 1,n
 x_mat(I,1) = 1.
end do
x_mat(1:n,k) = x_vec
 
x_prod = matmul(transpose(x_mat),x_mat)
call inverse(x_prod,x_inv,k)

predict =  matmul(transpose(x_mat),l_y)  

coeffs = matmul(x_inv,predict)


end subroutine coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Aggregator
