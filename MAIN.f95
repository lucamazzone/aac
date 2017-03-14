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
double precision  :: s(Zsize,1),s_alt(Zsize,1),mzero

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

double precision:: V_e(vecinterp,Zsize),l_y(Zsize),coeffs(2),result,marginals
integer :: entering(Zsize),v_entry(Zsize),polprimeind3(vecinterp,Zsize)
double precision :: labpol_ent(vecinterp,Zsize),labentry(Zsize),polprimewgt3(vecinterp,Zsize)

!!

double precision :: labdist(vecinterp,Zsize), dist(vecinterp,vecinterp,Zsize)

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

!!!!!!! experiments for valfun

 Cons= 0.5
 Nbig = 0.53
 Ybig = 0.62
 Nbig_1 = 0.5
 Ybig_1 = 0.6
 Cons_1 = 0.37
 zeta1 = zeta(:,1)
 wage = (Cons**eta)*(Nbig**chi)
 mzero = 0.15
 
 
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
!$OMP PARALLEL PRIVATE(kkk,q_q,jjj,obj) SHARED(qq,zeta1,Ybig,l_grid,b_grid)
!$OMP DO 
    do iii=1,nn
	kkk = (iii+vecsize**2-1)/(vecsize**2)
	q_q = reshape(qq(:,:,kkk),(/vecsize**2/))
	jjj = mod(iii-1,vecsize**2)+1
	obj = objectif(jjj,kkk,kappa,gamma,alpha,beta,zeta1,Ybig,wage,Q,q_q,l_grid,b_grid,expv0,vecsize,Zsize)
	vvalue(iii) = maxval(obj)
!	politics(iii) = maxloc(obj,1)
	policcc = maxloc(obj)
	politics(iii,1) = policcc(1)
    end do
!$OMP END DO
!	print*, 'operating thread n', omp_get_thread_num(), 'of', omp_get_num_threads()
!$OMP END PARALLEL
	vvalue0 = reshape(value0,(/nn/))
	epsilon = norm2(vvalue0-vvalue)
if (epsilon<0.0001)then
!print*, 'exiting!'
exit
else
    value0 = reshape(vvalue,(/vecsize**2,Zsize/))
!    print*, iter
end if
end do 


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
call convertpolicy2(polprimeind3,polprimewgt3,labpol_ent,lgrid_int(:,1))

call find_distribution(polprimeind1,polprimewgt1,polprimeind2,polprimewgt2,polprimeind3,polprimewgt3,&
& Zprob,mzero,entering,dist,labdist)

print*, labdist

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
integer :: zct,kct,ind
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

subroutine convertpolicy2(polprimeind,polprimewgt,policy,grid)
double precision, intent(in) :: policy(vecinterp,Zsize),grid(vecinterp)
double precision, intent(out) :: polprimewgt(vecinterp,Zsize)
integer, intent(out) :: polprimeind(vecinterp,Zsize)
integer :: kct, zct,ind
double precision :: wgt, primeval
do kkk=1,Zsize
    do iii  = 1,vecinterp
	primeval = policy(iii,kkk)
	ind = 1
	call hunt(grid,vecinterp,primeval,ind)
	if (ind<1) then
	    wgt = 0.0
	    ind = 1
	else if (ind >0 .AND. ind<vecinterp) then
	    wgt = (primeval-grid(ind))/(grid(ind+1)-grid(ind))
	else if (ind > vecinterp-1) then
	    wgt = 1.0
	    ind = vecinterp-1
	end if
	polprimeind(iii,kkk) = ind
	polprimewgt(iii,kkk) = wgt
    end do
end do



end subroutine convertpolicy2

subroutine find_distribution(polprimeind1,polprimewgt1,polprimeind2,polprimewgt2,polprimeind3,&
& polprimewgt3,Zprob,mzero,entering,dist,labdist)
implicit none
integer, intent(in) :: polprimeind1(vecinterp,vecinterp,Zsize),polprimeind2(vecinterp,vecinterp,Zsize)
integer, intent(in) :: polprimeind3(vecinterp,Zsize),entering(Zsize)
double precision, intent(in) :: polprimewgt1(vecinterp,vecinterp,Zsize),polprimewgt2(vecinterp,vecinterp,Zsize)
double precision, intent(in) :: polprimewgt3(vecinterp,Zsize),Zprob(Zsize,Zsize),mzero
double precision, intent(out) :: dist(vecinterp,vecinterp,Zsize),labdist(vecinterp,Zsize)

!!! other declarations

double precision :: distold(vecinterp,vecinterp,Zsize),onedimdistold(vecinterp,Zsize), summationold,disterror,disterror2
integer :: loop,ind1,ind2,zprime,ind3
double precision :: wgt1,wgt2,wgt3,summation,onedimdist(vecinterp,Zsize),dis_t(vecsize**2),onedimactivedist(vecinterp,Zsize)

! 1. Initialize Distribution

distold = 1.0
summationold = sum(reshape(distold,(/vecinterp*vecinterp*Zsize/)))

do kkk = 1,Zsize
    do jjj = 1,vecinterp
	do iii = 1,vecinterp
	    distold(iii,jjj,kkk) = distold(iii,jjj,kkk)/summationold
	end do
    end do
end do

do kkk = 1,Zsize
    do iii = 1,vecinterp
	onedimdistold(iii,kkk) = sum(distold(iii,:,kkk))
    end do
end do


dist = 1/(Zsize*vecinterp**2)

do kkk = 1,Zsize
    do iii = 1,vecinterp
	onedimdist(iii,kkk) = sum(dist(iii,:,kkk))
    end do
end do

! 2. Computing Distribution

loop = 0
disterror = 4.0

do while(disterror>0.01)

    loop = loop+1
    
    do kkk=1,Zsize
	do jjj=1,vecinterp
	    do iii=1,vecinterp
		ind1=polprimeind1(iii,jjj,kkk)
		wgt1=polprimewgt1(iii,jjj,kkk)
		ind2=polprimeind2(iii,jjj,kkk)
		wgt2=polprimewgt2(iii,jjj,kkk)
		    do zprime=1,Zsize
	dist(ind1,ind2,zprime) = dist(ind1,ind2,zprime)+distold(iii,jjj,zprime)*Zprob(kkk,zprime)*(1-wgt1)*(1-wgt2)
	dist(ind1+1,ind2,zprime)= dist(ind1+1,ind2,zprime)+distold(iii,jjj,zprime)*Zprob(kkk,zprime)*wgt1*(1-wgt2)
	dist(ind1,ind2+1,zprime)= dist(ind1,ind2+1,zprime)+distold(iii,jjj,zprime)*Zprob(kkk,zprime)*(1-wgt1)*wgt2
	dist(ind1+1,ind2+1,zprime)=dist(ind1+1,ind2+1,zprime)+distold(iii,jjj,zprime)*Zprob(kkk,zprime)*wgt1*wgt2
		    end do	
	    end do
	end do
    end do


    do kkk=1,Zsize
	do jjj=1,vecinterp
	    do iii=1,vecinterp
		if (dist(iii,jjj,kkk).LT.0.0) then
		dist(iii,jjj,kkk) = 0.0
		end if
	    end do
	end do
    end do
    
    do kkk = 1,Zsize
    summation = sum( reshape(dist(:,:,kkk),(/vecinterp*vecinterp,1/)))
	do jjj=1,vecinterp
	    do iii=1,vecinterp
		dist(iii,jjj,kkk) = dist(iii,jjj,kkk)/summation
	    end do
	end do
    end do

    summation = sum(reshape(dist,(/Zsize*vecinterp**2/)))
    dist = dist/summation
    disterror = norm2(  reshape(dist,(/Zsize*vecinterp**2/)) - reshape(distold,(/Zsize*vecinterp**2/)) )
    distold = dist

    if (loop .GT. 40) then
    exit
    end if

end do

do kkk=1,Zsize
    do iii=1,vecinterp
	onedimactivedist(iii,kkk)=sum(dist(iii,:,kkk))
    end do
end do


loop = 0
disterror2=4.0

do while(disterror2 .GT. 0.001)
loop = loop+1

    do kkk=1,Zsize
        do iii=1,vecinterp
    	    ind3=polprimeind3(iii,kkk)
    	    wgt3=polprimewgt3(iii,kkk)
    		do zprime=1,Zsize
    		    onedimdist(ind3,zprime) = onedimdist(ind3,zprime)+onedimdistold(iii,zprime)*Zprob(kkk,zprime)*(1-wgt3)
    		    onedimdist(ind3+1,zprime)= onedimdist(ind3+1,zprime)+onedimdistold(iii,zprime)*Zprob(kkk,zprime)*wgt3
    		end do
        end do
    end do
    
    
    do kkk=1,Zsize
	do iii=1,vecinterp
	    if (onedimdist(iii,kkk) .LT. 0.0) then
	    onedimdist(iii,kkk) = 0.0
	    end if
	end do
    end do
    
    do kkk=1,Zsize
	summation = sum(onedimdist(:,kkk))
	do iii=1,vecinterp
	    onedimdist(iii,kkk)= onedimdist(iii,kkk)/summation
	end do
    end do
    summation = sum(onedimdist)
    onedimdist = onedimdist/summation
    
    disterror2 = norm2(reshape(onedimdist,(/vecinterp*Zsize/))- reshape(onedimdistold,(/vecinterp*Zsize/)) )
    onedimdistold = onedimdist
    
    if (loop .GT. 30) then
    exit
    end if 

end do


! 3. Aggregating


do kkk=1,Zsize
    labdist(:,kkk) = onedimdist(:,kkk)*mzero*entering(kkk) + (1-mzero)*onedimactivedist(:,kkk)
end do


do kkk=1,Zsize
    summation = sum(labdist(:,kkk))
    do  iii=1,vecinterp
	labdist(iii,kkk) = labdist(iii,kkk)/summation
    end do
end do

summation = sum( reshape(labdist,(/vecinterp*Zsize/)) )
labdist = labdist/summation

end subroutine find_distribution
end program MAIN

 
