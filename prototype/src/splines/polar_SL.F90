program radial_1d_SL
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use cubic_nonuniform_splines
  use numeric_constants
  implicit none


  type(cubic_nonunif_spline_1D), pointer :: spl_hrmt 
  sll_int32 :: N,Nr,Ntheta,i,j,err,test_case,step,nb_step,visu_step,nb_diag
  sll_real64,dimension(:,:), pointer :: f
  sll_real64 :: rmin,rmax,r,dr,x,y,dt,theta,dtheta,val,tmp
  sll_real64,dimension(:,:), pointer :: diag
  sll_real64 :: a1,a2,rr,thetath,k1r,k2r,k3r,k4r,k1theta,k2theta,k3theta,k4theta
  
  !parameters
  N=256
  Nr=128
  Ntheta=128
  rmin = 0.2_f64
  rmax = 0.8_f64
  test_case = 1
  dt = 0.1_f64
  nb_step = 50
  visu_step = 10
  nb_diag = 2
  
  a1=1._f64
  a2=2._f64
  
  
  print *,'#N=',N
  print *,'#T=',real(nb_step,f64)*dt
  
  dr = (rmax-rmin)/real(Nr,f64)
  dtheta = 2._f64*sll_pi/real(Ntheta,f64)
  SLL_ALLOCATE(f(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(diag(nb_diag,0:nb_step), err)
  
  if(test_case==1)then
    do i=1,Nr+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        theta = real(j,f64)*dtheta
      !f(i) = 1._f64
        !f(i,j) = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)
        x = r*cos(theta)
        y = r*sin(theta)
        r = sqrt(x*x+y*y)
        !f(i,j) = exp(-1000._f64*(x-0.35)**2-1000._f64*(y-0.35)**2)
        f(i,j) = exp(-3000._f64*(r-1.5_f64*rmin)**2)
        !f(i,j) = exp(-30._f64*(theta-0.5*sll_pi)**2)
      enddo  
    enddo
  endif
  
  
  tmp=0._f64
  do i=1,Nr!+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta!+1
      tmp = tmp + f(i,j)*r
    enddo
  enddo
  tmp = tmp*dr*dtheta
  
  diag = 0._f64
  diag(1,0) = tmp


  tmp=0._f64
  do i=1,Nr!+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta!+1
      tmp = tmp + f(i,j)
    enddo
  enddo
  tmp = tmp*dr*dtheta
  diag(2,0) = tmp
  
  print *,'#init mass=',diag(1,0)
  
  
  print *,'#f0 boundary:',f(1,1),f(Nr+1,1)
  !save f
  open(unit=900,file='f0.dat')  
    do i=1,Nr+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        theta = real(j,f64)*dtheta
        x = r*cos(theta)
        y = r*sin(theta)
        !write(900,*) r,theta,f(i,j)
        write(900,*) x,y,f(i,j)
      enddo
      write(900,*) ' '
    enddo
  close(900)
  
  
  i=1;j=Ntheta/2
  r=rmin+real(i-1,f64)*dr
  theta=real(j-1,f64)*dtheta

  open(unit=900,file='carac.dat')  
  open(unit=800,file='fn.dat')    
    do i=1,Nr!+1
      do j=1,Ntheta!+1
        r=rmin+real(i-1,f64)*dr
        theta=real(j-1,f64)*dtheta
        x = r*cos(theta)
        y = r*sin(theta)
        write(800,'(f0.3,a,f0.3)',advance='no') r,' ',theta
        write(800,'(a,f0.3,a,f0.3)',advance='no') ' ',x,' ',y
        do step=1,nb_step
          !x = r + r*cos(theta)*dt
          !y = theta+sin(theta)*dt
          !Explicit Euler
          rr = r + r*(r-rmin)*(rmax-r)*cos(theta)*dt
          !y = theta+(-r*(rmax+rmin-2._f64*r)+(r-rmin)*(rmax-r))*sin(theta)*dt
          thetath = theta-(r*(rmax+rmin-2._f64*r)+2._f64*(r-rmin)*(rmax-r))*sin(theta)*dt
          !r = rr
          !theta = thetath
          
          !RK4
          rr = r
          thetath = theta
          k1r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k1theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k1r
          thetath = theta +  0.5_f64*dt*k1theta
          k2r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k2theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k2r
          thetath = theta +  0.5_f64*dt*k2theta
          k3r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k3theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + dt*k3r
          thetath = theta +  dt*k3theta
          k4r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k4theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          r = r+dt/6._f64*(k1r+2._f64*k2r+2._f64*k3r+k4r)
          theta = theta+dt/6._f64*(k1theta+2._f64*k2theta+2._f64*k3theta+k4theta)
          
          
          
          
          
          !x = x + a1*dt
          !y = y + a2*dt
          
          
          
          !r = sqrt(x**2+y**2)
          !theta = acos(x/r)
          
          !r = x
          !theta = y
          
          x = r*cos(theta)
          y = r*sin(theta)
          !val = exp(-3000._f64*(r-1.5_f64*rmin)**2)
          
          val = exp(-3000._f64*(r-1.5_f64*rmin)**2)
          !val = exp(-30._f64*(theta-0.5*sll_pi)**2)
          !val = exp(-1000._f64*(x-0.35)**2-1000._f64*(y-0.35)**2)
          diag(1,step) = diag(1,step)+ val *(rmin+real(i-1,f64)*dr)*dr*dtheta
          diag(2,step) = diag(2,step)+ val *dr*dtheta
          if(mod(step,visu_step)==0)then
            write(900,*) x,y
            write(800,'(a,f0.3)',advance='no') ' ',val
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
  
  open(unit=900,file='diag.dat')
    do step=1,nb_step
      write(900,*) real(step,f64)*dt,diag(1,step),diag(1,step)/diag(1,0)-1._f64,diag(2,step),diag(2,step)/diag(2,0)-1._f64
    enddo
  close(900)

  
  

end program