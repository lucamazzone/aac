  program MAIN 
  use params
  use library
  use solution_lib
  use omp_lib


  implicit none
  
!!!!!!!!!!!!!
integer, parameter :: snum = 2
integer :: zct,aggregate,curr_state,loop
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

double precision:: V_e(vecinterp,Zsize),l_y(Zsize),coeffs(2),result,marginals,size_active,implied_consumption,C_high,C_low
integer :: entering(Zsize),v_entry(Zsize),polprimeind3(vecinterp,Zsize)
double precision :: labpol_ent(vecinterp,Zsize),labentry(Zsize),polprimewgt3(vecinterp,Zsize),Nagg,Y_agg,Bagg,err_cons
double precision :: labdist(vecinterp,Zsize), dist(vecinterp,vecinterp,Zsize),nprimesimp(nsimp+1,Zsize),nprime(vecinterp,Zsize)
double precision :: N_1,Y_1,epsiloun,C_pred
double precision, allocatable :: momstoremat(:,:)


!! allocate where needed

allocate(momstoremat(Zsize,momnum))

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

 
 Nbig = 0.72
 Ybig = 0.7
 N_1 = 0.73
 Y_1 = 0.72
 C_pred = 0.7
 Cons_1 = 0.72
 zeta1 = zeta(:,1)
 mzero = 0.15
 

 loop = 0
 epsiloun = 1.0
 
do while( epsiloun > 0.01) 

C_low = 0.4
C_high = 1.1
loop = 0
 
do while( abs(C_high-C_low) .GT. 0.01 )
    loop = loop+1
    print*,loop
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
	    policcc = maxloc(obj)
	    politics(iii,1) = policcc(1)
	    end do
    !$OMP END DO
	    !print*, 'operating thread n', omp_get_thread_num(), 'of', omp_get_num_threads()
    !$OMP END PARALLEL
	    vvalue0 = reshape(value0,(/nn/))
	    epsilon = norm2(vvalue0-vvalue)
		if (epsilon<0.0001)then
		!print*, 'exiting!'
		exit
		else
		value0 = reshape(vvalue,(/vecsize**2,Zsize/))
		!print*, iter
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
    call transform_simp(nprimesimp,nprime,Zprob,polprimeind1,polprimewgt2,lgrid_int,nodes)
    Nagg = 0.0
    Y_agg = 0.0
    Bagg = 0.0
    call aggregate_var(Nagg,Y_agg,Bagg,nprime,dist,labdist,labpol_int,debpol_int,zeta1,alpha,gamma)
    size_active = 1.0-mzero+mzero*dot_product(s(:,1),entering) + marginals
    Nbig_1 = Nagg*size_active
    Ybig_1 = (size_active*Y_agg)**(gamma/(gamma-1))
    implied_consumption = Ybig - mzero*csi*dot_product(s(:,1),entering)

    if (Cons-implied_consumption .LT. 0.0) then
	C_low = Cons
    elseif (Cons-implied_consumption .GT. 0.0) then
	C_high = Cons
    end if

end do

epsiloun = abs(N_1 - Nbig_1) + abs(Y_1 - Ybig_1) + abs(C_pred - Cons)

N_1 = 0.5*N_1 + 0.5*Nbig_1
Y_1 = 0.5*Y_1 + 0.5*Ybig_1
C_pred = 0.5*C_pred + 0.5*Cons
Cons_1 = C_pred*Y_1/Ybig    

end do


do kkk=1,Zsize
    momstoremat(kkk,1) = dot_product(labdist(:,kkk),lgrid_int(:,1))/sum(labdist(:,kkk))
end do

do jjj=2,momnum
    do kkk=1,Zsize
	momstoremat(kkk,jjj) = dot_product(labdist(:,kkk),(lgrid_int(:,1)-momstoremat(kkk,1) )**jjj)/sum(labdist(:,kkk))
    end do
end do




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
integer :: ind1,zprime
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

subroutine pwl_basis_1d ( nd, xd, ni, xi, bk )

!*****************************************************************************80
!
!! PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) BK(NI,ND), the basis functions at the 
!    interpolation points.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) bk(ni,nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)

  bk(1:ni,1:nd) = 0.0D+00

  if ( nd == 1 ) then
    bk(1:ni,1:nd) = 1.0D+00
    return
  end if

  do i = 1, ni

    do j = 1, nd

      if ( j == 1 .and. xi(i) <= xd(j) ) then

        t = ( xi(i) - xd(j) ) / ( xd(j+1) - xd(j) )
        bk(i,j) = 1.0D+00 - t

      else if ( j == nd .and. xd(j) <= xi(i) ) then

        t = ( xi(i) - xd(j-1) ) / ( xd(j) - xd(j-1) )
        bk(i,j) = t

      else if ( xd(j-1) < xi(i) .and. xi(i) <= xd(j) ) then

        t = ( xi(i) - xd(j-1) ) / ( xd(j) - xd(j-1) )
        bk(i,j) = t

      else if ( xd(j) <= xi(i) .and. xi(i) < xd(j+1) ) then

        t = ( xi(i) - xd(j) ) / ( xd(j+1) - xd(j) )
        bk(i,j) = 1.0D+00 - t

      end if

    end do

  end do

  return
end subroutine pwl_basis_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


subroutine aggregate_var(Nagg,Y_agg,Bagg,nprime,dist,labdist,labpol_int,debtpol_int,zeta,alpha,gamma)
implicit none
double precision, intent(in) :: nprime(vecinterp,Zsize),dist(vecinterp,vecinterp,Zsize),labdist(vecinterp,Zsize),gamma
double precision, intent(in):: labpol_int(vecinterp,vecinterp,Zsize),debtpol_int(vecinterp,vecinterp,Zsize),zeta(Zsize),alpha
double precision, intent(out) :: Nagg,Y_agg,Bagg


Nagg=0.0
Y_agg=0.0
Bagg=0.0

do kkk=1,Zsize
    do jjj=1,vecinterp
	do iii=1,vecinterp
	    Bagg = Bagg + dist(iii,jjj,kkk)*debtpol_int(iii,jjj,kkk)
	end do
    end do
end do


do kkk=1,Zsize
    do iii=1,vecinterp
	Nagg = Nagg + labdist(iii,kkk)*nprime(iii,kkk)
	Y_agg = Y_agg + zeta(kkk)*labdist(iii,kkk)*(nprime(iii,kkk))**(alpha*(gamma-1)/gamma)
    end do
end do



end subroutine aggregate_var


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function Fk(kval,zct)
implicit none

!this function evaluates the density-proportional function
!Fk = exp(rho_1 * (k-m^zct_1)+....+rho_momuse * ((k - m^zct_1)^momuse - m^zct_momuse) )

double precision :: kval
integer :: zct

double precision :: rhovec(momuse),momvec(momuse)
integer :: momct

!extract correct value of rho from global variable
rhovec = rhomat(zct,:)
momvec = mommat(zct,:)

Fk = rhovec(1)*(kval - momvec(1))

do momct=2,momuse
    Fk = Fk + rhovec(momct) * ( (kval-momvec(1))**dble(momct) - momvec(momct) )
end do !momct

Fk = exp(Fk)

end function Fk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Pint()
implicit none

!Pint is \sum_{z=1}^n_z \int_kmin^kmax Fk(k,z) dk, where Fk is above, i.e. Pint is the objective that
!is minimized in the process of finding distributions matching the indicated moments

double precision :: Pint

integer :: zct,kct
double precision :: kval,wgt,addval

Pint = 0.0
intvec(:) = 0.0
do zct=1,znum
    do kct=1,nsimp+1
        
        kval = simpnodes(kct)
        wgt = simpweights(kct)
        addval = wgt * Fk(kval,zct)
        
        Pint = Pint + addval
        intvec(zct) = intvec(zct) + addval
    end do !kct
end do !zct

end function Pint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function gradPint()
implicit none

!gradPint is the znum * momuse x 1 gradient of Pint() from above, wrt to {rho^z_1,....,rho^1_momuse }_{z=1}^znum
    
double precision :: gradPint(znum*momuse)    

integer :: momct,zct,kct,gradct
double precision :: rhovec(momuse),momvec(momuse),kval,wgt

gradPint(:) = 0.0

do zct=1,znum

!extract moments and rho's for zct    
rhovec = rhomat(zct,:)
momvec = mommat(zct,:)
    
do momct=1,momuse
    !track location in the gradient
    gradct = (zct-1)*momuse + momct
    
    !perform integration for each entry
    do kct=1,nsimp+1
        
        kval = simpnodes(kct)
        wgt = simpweights(kct)
        if (momct>1) then
            gradPint(gradct) = gradPint(gradct) + wgt * ( (kval - momvec(1))**dble(momct) - momvec(momct) ) * Fk(kval,zct)
        else if (momct==1) then
            gradPint(gradct) = gradPint(gradct) + wgt * (kval - momvec(1)) * Fk(kval,zct)
        end if
        
    end do !kct
    
    
end do !momct
end do !zct
    
end function gradPint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


double precision function Pintrhovec(rhovec)
implicit none

double precision, intent(in) :: rhovec(znum*momuse)

integer :: zct,momct,ct
ct=0
do zct=1,znum
do momct=1,momct
    ct=ct+1
    rhomat(zct,momct) = rhovec(ct)
end do !momct
end do !zct

Pintrhovec = Pint()
    
end function Pintrhovec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine findrhoBroyden()
implicit none

double precision :: newJinv(znum*momuse,znum*momuse),oldJinv(znum*momuse,znum*momuse),oldbigrho(znum*momuse),&
    newbigrho(znum*momuse),oldgrad(znum*momuse),newgrad(znum*momuse),diffrho(znum*momuse),diffgrad(znum*momuse),&
    numval(znum*momuse),denomval,multval(znum*momuse)
double precision :: rhoerror,graderror,funcerror,oldfunc,newfunc
integer :: zct,momct,gradct,iter,ct1,ct2
   
!make initial guesses for rho and evaluate gradient at those guesses

do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -1.0
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    oldbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = oldbigrho(gradct)
end do !momct
end do !zct
oldgrad = gradPint()
oldfunc = Pint()

do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -5.0
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    newbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = newbigrho(gradct)
end do !momct
end do !zct
newgrad = gradPint()  
newfunc = Pint()

!what are the implied initial differences in rho and the gradient?
diffrho = newbigrho - oldbigrho
diffgrad = newgrad - oldgrad

!make initial guesses for inv Hessian approx - identity matrix
oldJinv(:,:) = 0.0
do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    oldJinv(gradct,gradct) = 1.0
end do !momct
end do !zct

!actually do the Broyden iterations
do iter=1,maxbroydit
    
    !set some stuff to zero
    numval(:) = 0.0
    denomval = 0.0
    
    !determine numerator
    numval = matmul(oldJinv,diffgrad)
    numval = diffrho - numval

    !determine denominator
    multval = matmul(oldJinv,diffgrad)
    denomval = sum(diffrho*multval)
    
    !determinemultval
    multval(:) = 0.0
    multval = matmul(diffrho,oldJinv)
    
    !what is new guess for inv Hessian?
    do ct1=1,znum*momuse
    do ct2=1,znum*momuse
        newJinv(ct1,ct2) = oldJinv(ct1,ct2) + (1.0 / denomval) * numval(ct1)* multval(ct2)
    end do !ct1
    end do !ct2
    
    !what is new guess for rho?
    oldbigrho = newbigrho
    oldgrad = newgrad
    oldJinv = newJinv
    oldfunc = newfunc
    
    !take new step, at fixed stepsize length
    newbigrho = newbigrho - stepsize*matmul(newJinv,newgrad)
       
    !what is the new value for the gradient, i.e. evaluate gradient at newbigrho
    do zct=1,znum
    do momct=1,momuse
        gradct = (zct-1)*momuse + momct
        rhomat(zct,momct) = newbigrho(gradct)
    end do !momct
    end do !zct
    newgrad = gradPint()
    newfunc = Pint()

    rhoerror = maxval(abs(newbigrho - oldbigrho))/maxval(abs(oldbigrho))
    graderror = maxval(abs(newgrad))
    
    if (iter>5000) then
    funcerror = abs(1.0 - newfunc/oldfunc)
    else if (iter<=5000) then
    funcerror = 1.0
    end if

    if (mod(iter,1000)==1) then
        write(*,"(A,I5,A,F7.4,A,F7.4,A,F7.4,A,F7.4)") "it = ",iter,", rho = ",rhoerror,", grad = ",&
            graderror,", fun = ",funcerror,", val = ",Pint()
    end if
    
    if (iter>10) then
    if (graderror<broydengradtol.or.funcerror<broydenfunctol.or.rhoerror<broydenrhotol) exit
    end if
    
    diffgrad = newgrad - oldgrad
    diffrho = newbigrho - oldbigrho
    
end do !iter

write(*,"(A,I5,A,F7.4,A,F7.4,A,F7.4,A,F7.4)") "it = ",iter-1,", rho = ",rhoerror,", grad = ",&
            graderror,", fun = ",funcerror,", val = ",Pint()

end subroutine findrhoBroyden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine findrhoSR1wolfe(initflag)
implicit none

double precision :: newJinv(znum*momuse,znum*momuse),oldJinv(znum*momuse,znum*momuse),oldbigrho(znum*momuse),&
    newbigrho(znum*momuse),oldgrad(znum*momuse),newgrad(znum*momuse),diffrho(znum*momuse),diffgrad(znum*momuse),&
    numval(znum*momuse),denomval,nummat(znum*momuse,znum*momuse),stepdirec(znum*momuse)
double precision :: rhoerror,graderror,funcerror,oldfunc,newfunc,stepval
integer :: zct,momct,gradct,iter,ct1,ct2,stepct,initflag
logical :: wolfe1,wolfe2,nannewgrad,nannewfunc

!make initial guesses for rho and evaluate gradient at those guesses

if (initflag==0) then
do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -1.0
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

else if (initflag==1) then

    rhomat = rhoSS

end if

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    oldbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = oldbigrho(gradct)
end do !momct
end do !zct
oldgrad = gradPint()
oldfunc = Pint()

if (initflag==0) then
do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -1.3
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

else if (initflag==1) then

    rhomat = rhoSS*1.1

end if

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    newbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = newbigrho(gradct)
end do !momct
end do !zct
newgrad = gradPint()  
newfunc = Pint()

!what are the implied initial differences in rho and the gradient?
diffrho = newbigrho - oldbigrho
diffgrad = newgrad - oldgrad

!make initial guesses for inv Hessian approx - identity matrix
oldJinv(:,:) = 0.0
do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    oldJinv(gradct,gradct) = 1.0
end do !momct
end do !zct

!actually do the SR1 iterations
do iter=1,maxbroydit
    
    !set some stuff to zero
    numval(:) = 0.0
    nummat(:,:) = 0.0
    denomval = 0.0
    
    !determine numerator vector
    numval = matmul(oldJinv,diffgrad)
    numval = diffrho - numval
    
    !determine numerator matrix
    do ct1=1,znum*momuse
    do ct2=1,znum*momuse
        nummat(ct1,ct2) = numval(ct1) * numval(ct2)
    end do
    end do
    
    !determine denominator value
    do ct1=1,znum*momuse
        denomval = denomval + numval(ct1)*diffgrad(ct1)
    end do 
    
    !what is new guess for inv Hessian?
    newJinv = oldJinv + (1.0 / denomval) * nummat
    
    !what is new guess for rho?
    oldbigrho = newbigrho
    oldgrad = newgrad
    oldJinv = newJinv
    oldfunc = newfunc
    
    !choose a step length, according to Wolfe conditions
    stepdirec = matmul(newJinv,newgrad)
    stepval = stepinit
    do stepct=1,maxstepnum
        
        newbigrho = oldbigrho - stepval*stepdirec
           
        !what is the new value for the gradient, i.e. evaluate gradient at newbigrho
        do zct=1,znum
        do momct=1,momuse
            gradct = (zct-1)*momuse + momct
            rhomat(zct,momct) = newbigrho(gradct)
        end do !momct
        end do !zct
        newgrad = gradPint()
        newfunc = Pint()
        
        !check for nonexistence of new function or gradient value
        nannewfunc = isnan(newfunc)
        nannewgrad = .FALSE.
        do ct1=1,znum*momuse
            if (isnan(newgrad(ct1))) nannewgrad = .TRUE.
        end do 
        
        !check wolfe sufficient decrease condition
        wolfe1 = ( newfunc <= (oldfunc + c1 * stepval * sum(stepdirec*oldgrad) ) )
        
        !check wolfe slope condition
        wolfe2 = ( sum(stepdirec * newgrad) >= c2 * sum(stepdirec * oldgrad ) )
        
        !if existence and wolfe conditions hold, then done 
        if ((nannewfunc.eqv..FALSE.).and.(nannewgrad.eqv..FALSE.).and.&
            (wolfe1.eqv..TRUE.).and.(wolfe2.eqv..TRUE.)) exit
        
        !if made it past, then you failed, so update
        if (stepct<(maxstepnum-1)) then
            stepval = stepval * stepfrac
        else if (stepct==(maxstepnum-1)) then
            stepval = stepfallback
        end if
        

    end do !stepct
    
    rhoerror = maxval(abs(newbigrho - oldbigrho))/maxval(abs(oldbigrho))
    graderror = maxval(abs(newgrad))
    funcerror = abs(newfunc-oldfunc)/abs(oldfunc)

    
    !call sleep(10)
    if (iter>10) then
    if (graderror<broydengradtol) exit
    end if
    
    diffgrad = newgrad - oldgrad
    diffrho = newbigrho - oldbigrho
    
end do !iter

write(*,"(A,I5,A,F7.4,A,F7.4,A,I5,A,F7.4)") "it = ",iter-1,", stepval = ",stepval,", grad = ",&
            graderror,", stepct = ",stepct-1,", val = ",Pint()

!check for issues and return SS coefficients if there are optimization failures
if ((nannewfunc.eqv..TRUE.).or.(nannewgrad.eqv..TRUE.)) then
    write(*,*) "ISSUE WITH OPTIMIZATION, RETURNING SS VALUES"
    rhomat = rhoSS
    intvec = intSS
end if


end subroutine findrhoSR1wolfe

end program MAIN

 
