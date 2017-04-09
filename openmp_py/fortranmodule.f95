module fortranmodule
use omp_lib

contains
    subroutine test(x,y)
	real(8), dimension(:), intent(in)   :: x
	real(8), dimension(:), intent(out)  :: y
	integer :: i,n
	integer :: num_threads
	n = size(x,1)
	
	
	!$OMP PARALLEL DO PRIVATE(I) FIRSTPRIVATE(N) SHARED(X,Y)
	do i=1,n
	    if(i == 1) then
		num_threads = OMP_get_num_threads()
		print*, 'num_threads running:' , num_threads
	    end if 
	    y(i) = sin(x(i)) + cos(x(i)) + exp(x(i)) + log(x(i))
	end do
	!$OMP END PARALLEL DO
    end subroutine test

end module fortranmodule