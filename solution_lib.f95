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
  real ( kind = 8 ) r8_huge
  integer ( kind = 4 ) r8vec_bracket5
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
      + ( yd(j) - yd(j+1) ) * ( xi(i) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

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

end module solution_lib



