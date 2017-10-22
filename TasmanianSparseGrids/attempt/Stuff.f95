module Stuff

contains

subroutine f_kk(kval,zct,rhomat,momstoremat,Zsize,momnum,F_k)
implicit none
integer, intent(in) :: Zsize,momnum,zct
double precision, intent (in) :: kval, rhomat(Zsize,momnum), momstoremat(Zsize,momnum)
double precision, intent (out) :: F_k
integer :: momct
double precision :: rhovec(momnum), momvec(momnum)

rhovec = rhomat(zct,:)
momvec = momstoremat(zct,:)

F_k = rhovec(1)*(kval-momvec(1))

do momct = 2,momnum
	F_k = F_k + rhovec(momct) * ( (kval-momvec(1))**dble(momct) - momvec(momct) )
end do

F_k = exp(F_k)

end subroutine f_kk


end module Stuff
