

subroutine funfun(points,vals)

double precision, intent(in) :: points(2)
double precision, intent(out) :: vals

vals  = exp(-points(1))/(1.0 + 100.0*exp(-10.0*points(2)) )

end subroutine funfun

