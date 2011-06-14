! gfortran qnefspl.f90 bsplvd.f90 bsplvb.f90 test_qnefspl.f90 -llapack -ldfftpack
program test_quasi_neutral
  use sll_quasi_neutral_solver
  implicit none

  type(quasi_neutral_plan), pointer :: qndat
  integer, parameter     :: k=5  ! spline degree
  real(8), parameter     :: eps=1.d-14
  real(8), parameter     :: pi = 3.1415926535897931_8
  integer                :: nr, ntheta ! number of points 
  real(8)                :: dr, dtheta ! cell size
  real(8)                :: rmin, rmax, Lr   ! computational domain in r
  integer                :: nmode !, sec ! sec is not used
  real(8)                :: err, ri, rg, ttg, cr, sr, st, coef, fij
  real(8), dimension(:), allocatable :: mth, kth
  real(8), dimension(:,:), allocatable :: mar, kar
  real(8), dimension(:,:), allocatable :: Uex, U, C, F  ! Exact solution, solution spline coefficients and rhs of system
  integer :: i,j, ii, jj, jjj, ig, jg  ! loop indices
  real(8), dimension(k+1) :: biatr,biatth
  real :: etime, t1, t2, elapsed(2)  ! timing 

  print*, '----------------------------------------------------'
  print*, 'test spline QN solver. Spline degree = ',k
  print*, '----------------------------------------------------'

  ! Test intialisation
  ! Allocation
  nr = 10
  rmin = 0.2_8
  rmax = 1.0_8
  dr=(rmax-rmin)/(nr-1)
  ntheta = 16
  dtheta = 1.0_8/ntheta 
  
  allocate(mth(k+1))
  allocate(kth(k+1))
  allocate(mar(k+1,nr+k-1))
  allocate(kar(k+1,nr+k-1))
  ! enter correct matrices (nr=10, ntheta=16, rmin=0.2, rmax=1.0, dtheta=1/ntheta)
  if (k==1) then
     mth = (/0.041666666666666671_8,0.010416666666666664_8/)
     kth = (/32._8,-16._8/)
  elseif (k==2) then
     mth = (/0.034375_8, 0.013541666666666666667_8, 0.0005208333333333333333333_8/)
     kth = (/16.0_8, -5.33333333333333333_8, -2.6666666666666666_8/)
  elseif (k==3) then
     mth = (/0.029960317460317461_8, 0.014769345238095241_8, 0.0014880952380952384_8, 1.2400793650793656e-05_8/)
     kth = (/10.6666666666666666666_8,-2.0_8,-3.2_8,-0.1333333333333333333_8/)
  endif
   ! call intialization routine
  qndat => new_qn_plan( k, rmin, nr, ntheta, dr, dtheta )
  ! check errors
  err = maxval(abs(qndat%mth-mth))
  if (err<eps) then
     print*, 'degree=',k,', mth ok'
  else
     print*, 'test_qnefspl: problem with mth. Diff=',err 
  endif
  err = maxval(abs(qndat%kth-kth))
  if (err<eps) then
     print*, 'degree=',k,', kth ok'
  else
     print*, 'test_qnefspl: problem with kth. Diff=',err 
  endif

  
  ! test solver
  !------------
  nr = 128
  rmin = 0.2_8
  rmax = 1.0_8
  Lr = rmax-rmin
  dr= Lr/(nr-1)
  ntheta = 128
  dtheta = 2*pi/ntheta 

qndat => new_qn_plan( k, rmin, nr, ntheta, dr, dtheta )
!  call new_qn(qndat, k, rmin, nr, ntheta, dr, dtheta)
  ! Allocation
  allocate(F(nr+k-1,ntheta))
  allocate(C(nr+k-1,ntheta))
  allocate(U(nr,ntheta))
  allocate(Uex(nr,ntheta))

  ! initialize rhs
  F(:,:) = 0.0
  nmode = 6
  ri = rmin
             
  t1=etime(elapsed)
  do jg = 1, k+1
     ttg =  0.5_8*dtheta*(qndat%xgauss(jg)+1.0_8)
     call bsplvb(qndat%knotsth,k+1,1,ttg,k+1,biatth)
     do i = 0, Nr-2 ! loop on cells (Nr-1 cells for Nr points)
        ri = rmin+i*dr
        do ig= 1, k+1
           rg = ri + 0.5_8*dr*(qndat%xgauss(ig)+1.0_8)  ! rescale Gauss points to be in interval [ri,ri+dr]
           call bsplvb(qndat%knotsr,k+1,1,rg,k+1+i,biatr)  
           sr = sin(pi*(rg-rmin)/Lr)
           cr = cos(pi*(rg-rmin)/Lr)
           coef = 0.25_8*dtheta*dr*qndat%wgauss(ig)*qndat%wgauss(jg)*rg 
           do j = 0, Ntheta-1              
              st = cos(nmode*(ttg+j*dtheta))   
              fij = (((nmode/rg)**2+(pi/Lr)**2)*sr - pi/(rg*Lr)*cr)*st
              do ii = 0, k
                 do jj = 0, k   
                    jjj = mod(j+jj,Ntheta)+1
                    F(i+ii+1,jjj) = F(i+ii+1,jjj) + coef*biatr(ii+1)*biatth(jj+1)*fij 
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
!-------------
!!$ do jg = 1, k+1
!!$     ttg =  0.5_8*dtheta*(qndat%xgauss(jg)+1.0_8)
!!$     call bsplvb(qndat%knotsth,k+1,1,ttg,k+1,biatth)
!!$     do i = 0, Nr-2 ! loop on cells (Nr-1 cells for Nr points)
!!$        ri = rmin+i*dr
!!$        do ig= 1, k+1
!!$           rg = ri + 0.5_8*dr*(qndat%xgauss(ig)+1.0_8)  ! rescale Gauss points to be in interval [ri,ri+dr]
!!$           call bsplvb(qndat%knotsr,k+1,1,rg,k+1+i,biatr)  
!!$           sr = sin(pi*(rg-rmin)/Lr)
!!$           cr = cos(pi*(rg-rmin)/Lr)
!!$           do ii = 0, k
!!$              do j = 0, Ntheta-1              
!!$                 st = cos(nmode*(ttg+j*dtheta))        
!!$                 do jj = 0, k            
!!$                    F(i+ii+1,mod(j+jj,Ntheta)+1) = F(i+ii+1,mod(j+jj,Ntheta)+1) + &
!!$                         0.25_8*dtheta*dr*qndat%wgauss(ig)*qndat%wgauss(jg)*biatr(ii+1)*biatth(jj+1)*rg &
!!$                         *(((nmode/rg)**2+(pi/Lr)**2)*sr - pi/(rg*Lr)*cr)*st
!!$                 enddo
!!$              enddo
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo
!------------
  t2 = etime(elapsed)
  print*, 'rhs time', t2-t1
  t1 =t2

  write (7,*) F
  ! solve system
  call apply_quasi_neutral_solver_plan(qndat,F,C)
  t2 = etime(elapsed)
  print*, 'solve time', t2-t1
  t1=t2

  ! compute solution at grid points from spline coefficients
  call evalsplgrid(qndat,C,U)
  t2 = etime(elapsed)
  print*, 'eval time', t2-t1
  t1=t2

  ! compute exact solution and error
  do j=1,ntheta
    do i=1,nr
       Uex(i,j)= sin(pi*(i-1)*dr/Lr)*cos(nmode*(j-1)*dtheta)   
    enddo
 enddo

 print*, 'error=', maxval(abs(U-Uex))
  
 write (8,*) U


  print*, 'end of test'

end program
