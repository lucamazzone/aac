program markov
integer, parameter :: start = 1
integer, parameter ::  samples = 15
integer, parameter ::  t = 30
double precision, allocatable ::  chains(:,:), Sprob(:,:) 
integer, parameter :: I = 8585739



allocate(chains(t-1,samples),Sprob(2,2))

Sprob(1,1) =   0.87158036240791925     
Sprob(2,1) =   0.12841963759208072      
Sprob(1,2) =   0.12841963759208075      
Sprob(2,2) =   0.87158036240791925

call srand(I)

call markovchain(chains,Sprob,start,samples,t)

open(unit=10006,file='chainsmat.txt',ACTION="write",STATUS="new")
do kkk=1,t-1
write(10006,'(*(F14.7))')(real(chains(kkk,jjj) ),jjj=1,samples)
end do
close(10006)

contains

	subroutine markovchain(chains,Sprob,start,samples,t)

        integer, intent(in) :: start,samples,t
        double precision, intent(out) :: chains(t-1,samples)
        double precision, intent(in)  :: Sprob(2,2)
        double precision :: X(t-1,1), triangular(2,2),cum(2,2),ppi(2+1,1),ppiz(2,1)
        double precision :: s(2,1),chain(t-1,1),rand
        integer :: state(t-1,2),V(2,1),iii,jjjj,kkk
	integer,parameter :: seed = 6477876

       	call srand(seed)
	print*, rand(),rand()        
	do iii = 1,samples
        
        do jjj = 1,t-1
        X(jjj,1) = rand()
        end do
        
        s = 0
        s(start,1) = 1
        triangular = 1.0
        triangular(2,1) = 0
        !ppi(1,1) = 0.0
        
        V(1,1) = 1
        V(2,1) = 2
        
        cum = matmul(Sprob,triangular)
!        print*, 'cum', cum


        state(1,:) = s(:,1)
        do kkk = 1,t-1
        ppiz = matmul(Sprob,s)   ! cum
!	print*, ppiz, 'ppiz'
        ppi(:,1) = [0.d0,ppiz(1,1),1.d0]
	!print*, 'ppi', ppi
        
        	!print*, 'x val'  , X(kkk,1)  !, 'big pp', ppi(jjj+1,1), 'small pp', ppi(jjj,1)
                !if ((X(kkk,1) .LT. ppi(jjj+1,1)) .AND. (X(kkk,1) .GT. ppi(jjj,1) )) then
                !print*, X(kkk,1), ppi(2,1),ppi(3,1)
                if (X(kkk,1) .GT. ppi(2,1) .AND. X(kkk,1) .LT. ppi(3,1) ) then
                !print*, 'x val', X(kkk,1), 'big  pp', ppi(jjj+1,1), 'small pp', ppi(jjj,1)
                !state(kkk,:) = s(:,1)
                s(2,1) = 1
                s(1,1) = 0
                else
                s(1,1) = 1
                s(2,1) = 0
                end if

            state(kkk,:) = s(:,1)
		chains(kkk,iii) = dot_product(state(kkk,:),V(:,1))
		!print*, chains(kkk,iii)
        end do
	!print*, 'chain', chains(:,iii)
	!chain(1,1) = start
        !chain = matmul(state(1,:),V(:,1))
	!chain(1,1)    =	DOT_PRODUCT(state(1,:), V(:,1))
	!print*, chain(1,1)
	!print*, chain
        !chains(:,iii) = chain(:,1)
        end do
        

	!print*, chains
	end subroutine

end program 
