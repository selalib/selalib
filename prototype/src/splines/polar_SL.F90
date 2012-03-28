program radial_1d_SL
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_splines
  use cubic_nonuniform_splines
  use numeric_constants
  implicit none

  type(sll_spline_2D), pointer :: spl_natper
  sll_int32  :: N,Nr,Ntheta,i,j,err,test_case,step,nb_step,visu_step,nb_diag,meth,jac,field_case
  sll_real64 :: rmin,rmax,r,dr,x,y,dt,theta,dtheta,val,tmp,tmp1,tmp2
  sll_real64 :: a1,a2,rr,thetath,k1r,k2r,k3r,k4r,k1theta,k2theta,k3theta,k4theta
  sll_real64,dimension(:,:), pointer :: f,fh,fh_buf
  sll_real64,dimension(:,:), pointer   :: rfeet,thetafeet
  sll_real64,dimension(:,:), pointer :: diag
	
  !parameters
  N=64
  Nr=N !128
  Ntheta=N!128
  rmin = 0.2_f64
  rmax = 0.8_f64
  test_case = 1
  field_case = 1
  dt = 1._f64/4._f64
  nb_step = 4
  visu_step = 1
  nb_diag = 4
  meth = 2
  jac = 1
	
  a1=-1.0_f64
  a2=-2.0_f64
  
  print *,'#N=',N
  print *,'#T=',real(nb_step,f64)*dt
  
  dr = (rmax-rmin)/real(Nr,f64)
  dtheta = 2._f64*sll_pi/real(Ntheta,f64)
	
  SLL_ALLOCATE(f(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(fh(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(fh_buf(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(rfeet(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(thetafeet(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(diag(nb_diag,0:nb_step), err)
	
  spl_natper => new_spline_2D(Nr+1, Ntheta+1, &
    rmin, rmax, &
    0._f64, 2._f64*sll_pi, &
    HERMITE_SPLINE, PERIODIC_SPLINE,&
    const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  
  do i=1,Nr+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta+1
      theta = real(j-1,f64)*dtheta
      x = r*cos(theta)
      y = r*sin(theta)
      r = sqrt(x*x+y*y)
			
      ! test-function
      if (test_case==1) then
        f(i,j) = exp(-3000._f64*(r-1.5_f64*rmin)**2)
      endif
      if (test_case==2) then
        f(i,j) = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)
      endif
      if (test_case==3) then
        f(i,j) = cos(theta)
      endif
      if (test_case==4) then
        if(r>=0.4.and.r<=0.6)then
          f(i,j) = exp(-30._f64*(theta-0.5*sll_pi)**2)
        else
          f(i,j) = 0._f64
        endif  
      endif
      if (test_case==5) then
        f(i,j) = 1._f64
      endif
    enddo  
  enddo

  
  fh = f
  
  ! diagnostic with the mass multiplied by the jacobien and the classical mass
  tmp1=0._f64
  tmp2=0._f64
  do i=1,Nr!+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta!+1
      tmp1 = tmp1 + f(i,j)*r
      tmp2 = tmp2 + f(i,j)
    enddo
  enddo
  tmp1 = tmp1*dr*dtheta
  tmp2 = tmp2*dr*dtheta
  diag(1,0) = tmp1
  diag(3,0) = tmp1
  diag(2,0) = tmp2
  diag(4,0) = tmp2
  print *,'#init mass=',diag(1,0)
  !checking that the function is close to 0 at the r boundary 
  print *,'#f0 boundary:',f(1,1),f(Nr+1,1)
  !saving f
  open(unit=900,file='f0.dat')  
    do i=1,Nr+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        theta = real(j-1,f64)*dtheta
        x = r*cos(theta)
        y = r*sin(theta)
        !write(900,*) r,theta,f(i,j)
        write(900,*) x,y,f(i,j)
      enddo
      write(900,*) ' '
    enddo
  close(900)
  
!  i=1;j=Ntheta/2
!  r=rmin+real(i-1,f64)*dr
!  theta=real(j-1,f64)*dtheta

  ! Analytic solution
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
					!y = theta+(-r*(rmax+rmin-2._f64*r)+(r-rmin)*(rmax-r))*sin(theta)*dt
					
          !Explicit Euler
          !rr = r + r*(r-rmin)*(rmax-r)*cos(theta)*dt
          !thetath = theta-(r*(rmax+rmin-2._f64*r)+2._f64*(r-rmin)*(rmax-r))*sin(theta)*dt
					!r = rr
          !theta = thetath
          
          !RK4
          if(field_case==1)then
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
          endif
          if(field_case==2)then
            theta = theta+dt
          endif
          if(field_case==3)then
            x = x + a1*dt
            y = y + a2*dt
            r = sqrt(x**2+y**2)
            if(y>=0)then
              theta = acos(x/r)
            else
              theta = 2._f64*sll_pi-acos(x/r)
            endif
          endif
          
          if(theta>2._f64*sll_pi)then
            theta= theta-2._f64*sll_pi
          endif
          if(theta<0._f64)then
            theta= theta+2._f64*sll_pi
          endif
          
          ! advection coef constant 
					!x = x + a1*dt
          !y = y + a2*dt
          
          !r = sqrt(x**2+y**2)
          !theta = acos(x/r)
          
          !r = x
          !theta = y
					



          
          x = r*cos(theta)
          y = r*sin(theta)
          					
					if (test_case==1) then
						val = exp(-3000._f64*(r-1.5_f64*rmin)**2)
					endif
					if (test_case==2) then
						val = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)
					endif
					if (test_case==3) then
						val = cos(theta)
					endif
					if (test_case==4) then
					  if(r>=0.4.and.r<=0.6)then
                        val = exp(-30._f64*(theta-0.5*sll_pi)**2)
                      else
                        val = 0._f64
                      endif
                      
					endif
                    if (test_case==5) then
                      val = 1._f64
                    endif
			
          diag(1,step) = diag(1,step) + val*(rmin+real(i-1,f64)*dr)*dr*dtheta
          diag(2,step) = diag(2,step) + val*dr*dtheta
					
          if(mod(step,visu_step)==0)then
            write(900,*) x,y
            write(800,'(a,f0.3)',advance='no') ' ',val
          endif
          f(i,j) = val
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
    do step=0,nb_step
      write(900,*) real(step,f64)*dt,diag(1,step),diag(1,step)/diag(1,0)-1._f64,diag(2,step),diag(2,step)/diag(2,0)-1._f64
    enddo
  close(900)
  
  
!  diag(1:nb_diag,1:nb_step) = 0._f64
  
  ! BSL
	if (meth==1) then
				
	do step=1,nb_step
	if(jac==1)then
	  do i=1,Nr
	    r=rmin+real(i-1,f64)*dr
	    do j=1,Ntheta
	      fh(i,j) = fh(i,j)*r
	    enddo
	    fh(i,Ntheta+1)=fh(i,1)
	  enddo
	endif
    call compute_spline_2D(fh,spl_natper)
    do i=1,Nr!+1
      do j=1,Ntheta!+1
        r=rmin+real(i-1,f64)*dr
        theta=real(j-1,f64)*dtheta
        x = r*cos(theta)
        y = r*sin(theta)    
				    
        !RK4
        if(field_case==1)then
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
        endif
        
        if(field_case==2)then
          theta = theta+dt
		endif

        if(field_case==3)then
            x = x + a1*dt
            y = y + a2*dt

            r = sqrt(x**2+y**2)
            if(y>=0)then
              theta = acos(x/r)
            else
              theta = 2._f64*sll_pi-acos(x/r)
            endif
            !if(abs(r-sqrt(x**2+y**2))>1.e-12)then
            !  print *,'r',r,sqrt(x**2+y**2),i,j
            !  stop
            !endif
            !if(abs(theta-acos(x/r))>1.e-12)then
            !  print *,theta,acos(x/r),i,j
            !  stop
            !endif
		endif

		
        if(theta>2._f64*sll_pi)then
          theta= theta-2._f64*sll_pi
        endif
        if(theta<0._f64)then
          theta= theta+2._f64*sll_pi
        endif
        
        
        val = interpolate_value_2D(r,theta,spl_natper)
				
        !if(abs(val-f(i,j))>1.e-12)then
        !  print *,i,j,fh(i,j),val,fh(i,j)-val,f(i,j),cos(theta)
        !endif

        if(jac==1)then
          val=val/(rmin+real(i-1,f64)*dr)
        endif
				
        fh(i,j) = val
        
				
        diag(3,step) = diag(3,step)+ val *(rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(4,step) = diag(4,step)+ val *dr*dtheta
        
      enddo
    enddo  
  enddo
	! L^{\infty} error 
  val = 0._f64
  do i=1,Nr
    do j=1,Ntheta
      val = max(val,abs(f(i,j)-fh(i,j)))
    enddo
  enddo  
  print*,val
	open(unit=900,file='diag2.dat')
    do step=0,nb_step
      write(900,*) real(step,f64)*dt,diag(3,step),diag(3,step)/diag(1,0)-1._f64,diag(4,step),diag(4,step)/diag(4,0)-1._f64
    enddo
  close(900)
	
	!FSL
	else
	
	do step=1,nb_step
	  if(jac==1)then
	    do i=1,Nr
	      r=rmin+real(i-1,f64)*dr
	      do j=1,Ntheta
	        fh(i,j) = fh(i,j)*r
	      enddo
	    enddo  
	  endif

    call compute_spline_2D(fh,spl_natper)
    do i=1,Nr+1
      do j=1,Ntheta+1
        r=rmin+real(i-1,f64)*dr
        theta=real(j-1,f64)*dtheta
        x = r*cos(theta)
        y = r*sin(theta)    
				    
        !RK4
        dt = -dt
        
        if(field_case==1)then
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
        endif
        
        
        if(field_case==2)then
          theta = theta+dt
		endif

        if(field_case==3)then
            x = x + a1*dt
            y = y + a2*dt

            r = sqrt(x**2+y**2)
            if(y>=0)then
              theta = acos(x/r)
            else
              theta = 2._f64*sll_pi-acos(x/r)
            endif
            !if(abs(r-sqrt(x**2+y**2))>1.e-12)then
            !  print *,'r',r,sqrt(x**2+y**2),i,j
            !  stop
            !endif
            !if(abs(theta-acos(x/r))>1.e-12)then
            !  print *,theta,acos(x/r),i,j
            !  stop
            !endif
		endif

		
        if(theta>2._f64*sll_pi)then
          theta= theta-2._f64*sll_pi
        endif
        if(theta<0._f64)then
          theta= theta+2._f64*sll_pi
        endif
        
        rfeet(i,j) = r
        thetafeet(i,j) = theta
        dt = -dt
      enddo
    enddo
    
    !print*,'r/theta-feet',rmin+real(17-1,f64)*dr,rfeet(17,12),real(12-1,f64)*dtheta,thetafeet(17,12)
    
		!print*,'avant',fh(17,12)
    call deposit_value_2D_bis(rfeet,thetafeet,spl_natper,fh)
    !call deposit_value_2D(rfeet,thetafeet,spl_natper,fh)
    
    if(jac==1)then
	  do i=1,Nr
	    r=rmin+real(i-1,f64)*dr
	    do j=1,Ntheta
	      fh(i,j) = fh(i,j)/r
	    enddo
	  enddo  
	endif

    
    !print*,'apres',fh(17,12)
    !stop
		do i=1,Nr!+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta!+1
        theta = real(j-1,f64)*dtheta
				
        diag(3,step) = diag(3,step)+ fh(i,j)*(rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(4,step) = diag(4,step)+ fh(i,j)*dr*dtheta
      enddo
    enddo
  enddo
  
  
  

  
! L^{\infty} error 
  val = 0._f64
  do i=1,Nr
    do j=1,Ntheta
      val = max(val,abs(f(i,j)-fh(i,j)))
    enddo
  enddo  
  print*,val
	open(unit=900,file='diag2.dat')
    do step=0,nb_step
      write(900,*) real(step,f64)*dt,diag(3,step),diag(3,step)/diag(1,0)-1._f64,diag(4,step),diag(4,step)/diag(4,0)-1._f64
    enddo
  close(900)
	
	end if
	
	


  open(unit=900,file='fh.dat')  
    do i=1,Nr!+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta!+1
        theta = real(j-1,f64)*dtheta
        x = r*cos(theta)
        y = r*sin(theta)
				
        !write(900,*) r,theta,f(i,j)
        write(900,*) r,theta,x,y,fh(i,j),f(i,j)
      enddo
      write(900,*) ' '
    enddo
  close(900)
end program