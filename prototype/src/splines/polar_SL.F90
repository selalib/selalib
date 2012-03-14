program radial_1d_SL
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use cubic_nonuniform_splines
  use numeric_constants
  implicit none


  type(cubic_nonunif_spline_1D), pointer :: spl_hrmt 
  sll_int32 :: N,Nr,Ntheta,i,j,err,test_case,step,nb_step,visu_step
  sll_real64,dimension(:,:), pointer :: f
  sll_real64 :: rmin,rmax,r,dr,x,y,dt,theta,dtheta

  !parameters
  N=128
  Nr=N
  Ntheta=N
  rmin = 0.2_f64
  rmax = 0.8_f64
  test_case = 1
  dt = 0.01_f64
  nb_step = 5000
  visu_step = 100
  
  print *,'#N=',N
  print *,'#T=',real(nb_step,f64)*dt
  
  dr = (rmax-rmin)/real(Nr,f64)
  dtheta = 2._f64*sll_pi/real(Ntheta,f64)
  SLL_ALLOCATE(f(Nr+1,Ntheta+1), err)
  
  if(test_case==1)then
    do i=1,Nr+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
      !f(i) = 1._f64
        !f(i,j) = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)
        f(i,j) = exp(-3000._f64*(r-1.5_f64*rmin)**2)
      enddo  
    enddo
  endif
  
  
  
  
  print *,'#f0 boundary:',f(1,1),f(Nr+1,1)
  !save f
  open(unit=900,file='f0.dat')  
    do i=1,Nr+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        theta = real(j,f64)*dtheta
        write(900,*) r,theta,f(i,j)
      enddo
      write(900,*) ' '
    enddo
  close(900)
  
  
  i=1;j=Ntheta/2
  r=rmin+real(i-1,f64)*dr
  theta=real(j-1,f64)*dtheta

  open(unit=900,file='carac.dat')  
  open(unit=800,file='fn.dat')    
    do i=1,Nr+1
      do j=1,Ntheta+1
        r=rmin+real(i-1,f64)*dr
        theta=real(j-1,f64)*dtheta
        write(800,'(f0.3,a,f0.3)',advance='no') r,' ',theta
        do step=1,nb_step
          !x = r + r*cos(theta)*dt
          !y = theta+sin(theta)*dt
          x = r + r*(r-rmin)*(rmax-r)*cos(theta)*dt
          y = theta+(-r*(rmax+rmin-2._f64*r)+(r-rmin)*(rmax-r))*sin(theta)*dt
          r = x
          theta = y
          x = r*cos(theta)
          y = r*sin(theta)
          if(mod(step,visu_step)==0)then
            write(900,*) x,y
            write(800,'(a,f0.3)',advance='no') ' ',exp(-1000._f64*(r-0.5_f64*(rmax+rmin))**2)
          endif  
        enddo
        write(800,*) ' '
        !f(i,j) = exp(-1000._f64*(r-0.5_f64*(rmax+rmin))**2)
      enddo
      write(900,*) ' '
      write(800,*) ' '
    enddo  
  close(900)
  close(800)

  
  

end program