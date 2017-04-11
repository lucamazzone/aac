program Aggregator

! change end program into end module at the end
use params
use solution_lib
use library

!contains
!subroutine mapping_inverse()

implicit none

!double precision, intent(in) :: points(2),mzero,Tol
!integer, intent(in) :: aggregate
!double precision, intent(out) :: vals(3)
double precision :: points(2),mzero,Tol,vals(3)
integer :: aggregate
double precision :: pred(3),threshold

!! other declarations
double precision :: intvecmat(snum,Zsize), distribution1(vecinterp,Zsize*vecinterp), distribution2(vecinterp,Zsize*vecinterp)
double precision :: rhomat(momnum*Zsize,snum), momentsmat(Zsize*momnum,snum)
!!!
integer :: zct,curr_state,loop
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
aggregate = 1
! intvector = intvecmat(:,aggregate)
! momstoremat = reshape(momentsmat(:,:,aggregate),(/Zsize,momnum/))
! rhomatrix = reshape(rhomat(:,:,aggregate),(/Zsize,momnum/))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
zeta1 = zeta(:,aggregate)
loop = 0

do while( epsiloun > threshold) !! Loop over expected future aggregates

    C_low = 0.4
    C_high = 1.1
    C_pred = pred(3)
    Cons_1 = pred(3)*Y_1/Ybig
    N_1 = pred(1)
    Y_1 = pred(2)

    loop=loop+1
    
    if (loop .GT. 30) then
      threshold = Tol/2
    else if (loop .GT. 50) then
      threshold = Tol
    else if (loop .GT. 100) then
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
      
      call convertpolicy3(polprimeind1,polprimewgt1,labpol_int,lgrid_int(:,1))
      call convertpolicy3(polprimeind2,polprimewgt2,debpol_int,bgrid_int(1,:))
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
      
      marginals = 0.0
      if (sum(entering(2:Zsize)-entering(1:Zsize-1))>0)  then    
    	call coeff(coeffs,Zsize,2,l_y,v_entry)
	    result = -coeffs(1)/coeffs(2)
	    marginals = marginals_entering(Zsize,result,v_entry,cums)
      end if
      
      do kkk=1,Zsize
	    labentry(kkk) = lgrid_int(v_entry(kkk),1)*entering(kkk)
      end do
      
      labpol_ent = 1.0
      labpol_ent(1,:) = labentry
      size_active = 1.0-mzero+mzero*dot_product(s(:,1),entering) + marginals
      implied_consumption = Ybig - mzero*csi*dot_product(s(:,1),entering)

      if (Cons-implied_consumption .LT. 0.0) then
	    C_low = Cons
      elseif (Cons-implied_consumption .GT. 0.0) then
	    C_high = Cons
      end if

    end do   !! end of mkt clearing loop
    
   call transform_simp(nprimesimp,nprime,Zprob,polprimeind1,polprimewgt2,lgrid_int,nodes) 
   
   
   
   
   
   

end do  !! end of expectations loop


   
    
      









!end subroutine mapping_inverse

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
     integer :: iii
     
     objectif = kappa*(zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-wage*l_grid(jjj)-&
&		 b_grid(jjj) + b_grid*q_q) + beta*(1-kappa)*Q*expv0(:,kkk)
     
     
     do iii=1,vecsize**2
	if (    kappa*(zeta1(kkk)*(Ybig**(1/gamma))*l_grid(jjj)**(alpha-alpha/gamma)-&
&	wage*l_grid(jjj) - b_grid(jjj) + b_grid(iii)*q_q(iii))  <   0) then
	 objectif(iii) = -100.0
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end program Aggregator
