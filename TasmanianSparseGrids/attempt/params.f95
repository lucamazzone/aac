module params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this module just lists parameters !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

  ! model  parameters

  double precision, parameter :: alpha = 0.55 !labor productivity
  double precision, parameter :: beta = 0.985 ! hh discount factor => r = 2,67% 
  double precision, parameter :: eta = 1.0 ! hh IES  = 1
  double precision, parameter :: chi = 0.5  ! Fritsch labor elasticity
  double precision, parameter :: gamma = 7.7 ! elasticity of subs b/w goods 7.7
  double precision, parameter :: rhoz = 0.7 ! serial corr of idiosync shocks
  double precision, parameter :: rhosigma = 0.75 ! serial corr of unc shocks
  double precision, parameter :: kappa = 0.4 ! Jensen effect
  double precision, parameter :: phi = 0.08 ! std of unc shocks  0.13
  double precision, parameter :: musigma = 0.09 ! mean of unc process 0.18
  double precision, parameter :: csi = 1.0  ! entry costs

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
  
  !! for loops
  integer, parameter :: t = 250
  integer, parameter :: samples = 200
  integer, parameter :: start = 1
 ! integer :: iii,jjj,kkk

  



  
 

 end module
