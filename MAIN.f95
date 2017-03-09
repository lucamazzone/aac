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

integer :: iter,maxiter,nn=Zsize*vecsize**2,politics(Zsize*vecsize**2,1),policcc(1)
double precision :: value0(vecsize**2,Zsize),expv0(vecsize**2,Zsize),vvalue(Zsize*vecsize**2),zeta1(Zsize)
double precision :: vvalue0(Zsize*vecsize**2),l_grid(vecsize**2),b_grid(vecsize**2),q_q(vecsize**2)
double precision :: obj(vecsize**2),epsilon

!!

double precision :: labpol(Zsize*vecsize**2), debpol(Zsize*vecsize**2),lab_pol(vecsize,vecsize,Zsize),deb_pol(vecsize,vecsize,Zsize)
double precision :: labpol_int(Zsize*vecinterp**2), debpol_int(Zsize*vecinterp**2),zi(vecinterp)


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


!!!!!!! experiments for valfun

 Cons= 0.6
 Nbig = 0.6
 Ybig = 0.7
 Nbig_1 = 0.6
 Ybig_1 = 0.6
 Cons_1 = 0.55
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


call pwl_interp_2d(vecsize,vecsize,lgrid(:,1),bgrid(1,:),lab_pol(:,:,1),vecinterp,lgrid_int(1,:),bgrid_int(:,60),zi)

print*, zi


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
	 objectif(iii) = -1500.0
         end if
     end do
     
     end function
     
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
double precision, intent(in) :: policy(vecsize,vecsize,Zsize),grid(vecsize)
double precision, intent(out) :: polprimewgt(vecsize,vecsize,Zsize)
integer, intent(out) :: polprimeind(vecsize,vecsize,Zsize)
integer :: zct,kct,ind
double precision :: wgt,primeval

!grid = lgrid_int(:,1);
!policy = labpol_int;
  
  
do zzz = 1,Zsize
	do jjj = 1,vecsize
    		do iii = 1,vecsize

primeval = policy(jjj,iii,zzz)
ind = 1
call hunt(grid,vecsize,primeval,ind) 

if (ind<1) then
    wgt = 0.0
    ind = 1
else if (ind > 0.0 .AND. ind < vecsize ) then
    wgt = (primeval - grid(ind))/(grid(ind+1)-grid(ind));
else if (ind > vecsize-1) then
    wgt = 1.0;
    ind = vecsize-1;
end
   
polprimeind(jjj,iii,zzz) = ind;
polprimewgt(jjj,iii,zzz) = wgt;

        end 
    end
end


    
 end subroutine convertpolicy3
    
    end program MAIN

 
