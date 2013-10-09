module bgk_mesh_construction
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_constants
  use cubic_non_uniform_splines
  !use utils
  implicit none

contains 
  subroutine construct_bgk_mesh(nc_eta1,nc_eta2,mesh_case,&
   &x1n_array,x2n_array,x1c_array,x2c_array,jac_array,integration_points,&
   &geom_x,geom_eta,alpha_mesh,N_x1,N_x2)
    use sll_constants
    implicit none
    sll_int32,intent(in)::nc_eta1,nc_eta2,mesh_case,n_x1,n_x2
    sll_real64,intent(in)::geom_x(2,2),geom_eta(2,2),alpha_mesh
    sll_real64,dimension(:,:),pointer::x1n_array,x2n_array,x1c_array,x2c_array
    sll_real64,dimension(:,:),pointer::jac_array
    sll_real64,dimension(:,:,:),pointer::integration_points
    sll_int32  :: i1,i2,err,i
    sll_real64 :: x1_min,x1_max,x2_min,x2_max,delta_x1,delta_x2,x1
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1,eta1c,eta2,eta2c
    sll_real64 :: val,tmp
    sll_real64 :: slope_mesh1,slope_mesh2,wk,ll,dxx,slope_mesh3
    sll_int    ::Nzon,Nzon2
    sll_real64 , dimension(4)         :: ws
    sll_int , dimension(4)         :: Ns

    
    

    SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x1c_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x2c_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(integration_points(3,N_x1+1,N_x2+1),err)

    
    x1_min=geom_x(1,1)
    x1_max=geom_x(2,1)
    x2_min=geom_x(1,2)
    x2_max=geom_x(2,2)

    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)

        


    delta_x1 = (x1_max-x1_min)/real(nc_eta1,f64)
    delta_x2 = (x2_max-x2_min)/real(nc_eta2,f64)
    
    delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
    delta_eta2 = (eta1_max-eta1_min)/real(nc_eta2,f64)
    
    
    if(mesh_case==1)then
      do i2=1,nc_eta2+1
        do i1=1,nc_eta1+1
          x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
          x2n_array(i1,i2) = x2_min+real(i2-1,f64)*delta_x2
          x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
          x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
          jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)
          !jac_array(i1,i2) = 1._f64!(x1_max-x1_min)*(x2_max-x2_min)
        enddo
      enddo
    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

      do i2=1,nc_eta2+1
        do i1=1,nc_eta1+1
          !eta1 value of intersecting point (eta2,x1)=constant
          integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
          !x2 value of intersecting point (eta2,x1)=constant
          integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
        enddo
      enddo
    
    
    
    endif

  if(mesh_case==2)then
    eta2 = 0.0_f64 
    eta2c = eta2 + 0.5_f64*delta_eta2
    do i2= 1, nc_eta2 + 1
      eta1 = 0.0_f64
      eta1c = 0.5_f64*delta_eta1
      do i1 = 1, nc_eta1 + 1
        x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
        x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
        x1c_array(i1,i2) = eta1c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
        x2c_array(i1,i2) = eta2c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
        !x1n_array(i1,i2) = (x1n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
        !x2n_array(i1,i2) = (x2n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
        !x1c_array(i1,i2) = (x1c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
        !x2c_array(i1,i2) = (x2c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
        jac_array(i1,i2) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)) * &
         (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c)) - &
         alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c) * &
         alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)
        eta1 = eta1 + delta_eta1
        eta1c = eta1c + delta_eta1
        x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
        x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
        x1c_array(i1,i2) = x1_min+x1c_array(i1,i2)*(x1_max-x1_min)
        x2c_array(i1,i2) = x2_min+x2c_array(i1,i2)*(x2_max-x2_min)
        jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
      end do
      eta2 = eta2 + delta_eta2
      eta2c = eta2c + delta_eta2
    end do

    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !  x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !  PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    val = 0._f64
    do i2=1,nc_eta2
      do i1=1,nc_eta1
        x1 = (real(i1,f64)-0.5_f64)/real(nc_eta1,f64)
        eta2 = (real(i2,f64)-0.5_f64)/real(nc_eta2,f64)
        tmp = alpha_mesh*sin(2._f64*sll_pi*eta2)
        do i=1,100
          val = val-(val+tmp*sin(2._f64*sll_pi*val)-x1)/&
           (1._f64+2._f64*sll_pi*tmp*cos(2._f64*sll_pi*val))
        enddo
        if(abs(val+tmp*sin(2._f64*sll_pi*val)-x1)>1.e-14)then
          print *,i1,i2,val+tmp*sin(2._f64*sll_pi*val)-x1,val
          print *,'Problem of convergence of Newton'
          stop
        endif
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = val
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+(x1-val+eta2)*(x2_max-x2_min)        
      enddo
    enddo
  endif



  if(mesh_case==3)then
     eta2 = eta2_min 
     eta2c = eta2_min + 0.5_f64*delta_eta2
     do i2= 1, nc_eta2 + 1
        eta1 = eta1_min
        eta1c = 0.5_f64*delta_eta1
        do i1 = 1, nc_eta1 + 1
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)**2
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x1c_array(i1,i2) = eta1c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)**2
           x2c_array(i1,i2) = eta2c + alpha_mesh * sin(2*sll_pi*eta1c)* sin(2*sll_pi*eta2c)
           !x1n_array(i1,i2) = (x1n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2n_array(i1,i2) = (x2n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x1c_array(i1,i2) = (x1c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2c_array(i1,i2) = (x2c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !jac_array(i1,i2) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)) * &
           !  (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c)) - &
           !  alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c) * &
           !  alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)
           !val =   1.0_f64 + 2._f64*alpha_mesh *sll_pi*sin(2*sll_pi*(eta1c+eta2c))
           !if(abs(jac_array(i1,i2)-val)>1e-13)then
           !  print *,jac_array(i1,i2),val
           !  stop
           !endif
           jac_array(i1,i2) = 1._f64+2._f64*sll_pi*alpha_mesh*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)&
           +2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1c)&
           -2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)**2&
           -4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)*cos(2._f64*sll_pi*eta1c)&
           +4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)**3*cos(2._f64*sll_pi*eta1c)
           eta1 = eta1 + delta_eta1
           eta1c = eta1c + delta_eta1
           x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
           x1c_array(i1,i2) = x1_min+x1c_array(i1,i2)*(x1_max-x1_min)
           x2c_array(i1,i2) = x2_min+x2c_array(i1,i2)*(x2_max-x2_min)
           jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
        end do
        eta2 = eta2 + delta_eta2
        eta2c = eta2c + delta_eta2
      end do
  
        
  open(unit=900,file='xn_array.dat')  
    do i1=1,N_x1
      !x1 = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
      x1 = x1_min+(real(i1,f64)-1._f64)*delta_x1
      do i2=1,N_x2
        !write(900,*) x1,integration_points(2,i1,i2),x1c_array(i1,i2),x2c_array(i1,i2)
        write(900,*) x1n_array(i1,i2),x2n_array(i1,i2)
      enddo  
    enddo
  close(900)
       
      

      !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
      !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
      !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
      !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

      !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
      !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,COMPACT)
      !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
      !   PERIODIC, eta2_min, eta2_max, nc_eta2, COMPACT, geom)



      !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

      val = 0._f64
      do i2=1,nc_eta2
        do i1=1,nc_eta1
          x1 = (real(i1,f64)-0.5_f64)/real(nc_eta1,f64)
          eta2 = (real(i2,f64)-0.5_f64)/real(nc_eta2,f64)
          tmp = alpha_mesh*sin(2._f64*sll_pi*eta2)**2
          do i=1,100
            val = val-(val+tmp*sin(2._f64*sll_pi*val)-x1)/&
            (1._f64+2._f64*sll_pi*tmp*cos(2._f64*sll_pi*val))
          enddo
          if(abs(val+tmp*sin(2._f64*sll_pi*val)-x1)>1.e-14)then
            print *,i1,i2,val+tmp*sin(2._f64*sll_pi*val)-x1,val
            print *,'Problem of convergence of Newton'
            stop
          endif
          
          !eta1 value of intersecting point (eta2,x1)=constant
          integration_points(1,i1,i2) = val
          
          
          !x2 value of intersecting point (eta2,x1)=constant
          integration_points(2,i1,i2) = x2_min+((x1-val)/sin(2._f64*sll_pi*eta2)+eta2)*(x2_max-x2_min)        
        
        enddo
        
        
      enddo
      !do i1=1,nc_eta1
      !  print *,i1,integration_points(i1,1), integration_points(i1,nc_eta2)
      !enddo
      !stop


    endif
  

   if(mesh_case==10)then
   ! 3 parts 
      slope_mesh1=1._f64/2._f64
      wk=2.82_f64
      ll=0.5_f64
      !stop
      Ns(1)=int((-wk-ll-x2_min)/delta_x2)+1
      ws(1)=x2_min+(Ns(1)-1)*delta_x2
      !Ns(2)=int((-wk+ll-x2_min)/delta_x2)+1
      !ws(2)=x2_min+(Ns(2)-1)*delta_x2
      
     ! Ns(3)=int((wk-ll-x2_min)/delta_x2)
      !ws(3)=x2_min+(Ns(3)-1)*delta_x2
      Ns(4)=int((wk+ll-x2_min)/delta_x2)+1
      ws(4)=x2_max-(ws(1)-x2_min)!x2_min+(Ns(4)-1)*delta_x2
     !iiw2=int((-wk-1-x2_min)/delta_x2)
    ! print*,ws(1)-ws(2),ws(3)-ws(4),delta_x2
     !stop
     !stop
     Nzon= int((ws(4)-ws(1))/(slope_mesh1*delta_x2))
     if(mod(N_x2-Nzon,2)==0) then
      Nzon=Nzon
     else
      Nzon=Nzon+1
     endif 
     ! print*,"slop",slope_mesh1,Nzon,N_x2-2*Nzon
      !stop
      !if(
      !Nzon= int(1._f64/(slope_mesh1*delta_x2))
      !slope_mesh2=(x2_max-x2_min-(ws(4)-ws(1)))/((N_x2-Nzon)*delta_x2)
        slope_mesh2=(ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64)
        slope_mesh1=(ws(4)-ws(1))/(Nzon*delta_x2)
       !slope_mesh2=(x2_max-x2_min-2)/((N_x2-2*Nzon)*delta_x2)
       !print*,(ws(4)-ws(1))/slope_mesh1+2*(ws(1)-x2_min)/slope_mesh2,(x2_max-x2_min)!, (ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64),Nzon,N_x2-Nzon
       !stop
       !(ws(2)-ws(1))-(ws(4)-ws(3)),
       !stop
      
      !do i1=1,nc_eta1+1  
       !do i2=1,nc_eta2+1  
         !x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       !enddo
      !enddo
     
 
     do i1=1,nc_eta1+1   

       dxx=(ws(1)-x2_min)/((N_x2-Nzon)/2._f64)
      
       do i2=1, (N_x2-Nzon)/2+1

       x2n_array(i1,i2) =x2_min+real(i2-1,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*dxx
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)*(ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64)
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/(slope_mesh2)
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)/slope_mesh2
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =x2_min+(real(i2,f64)-0.5_f64)*dxx
       
      enddo
     enddo

     do i1=1,nc_eta1+1   

       dxx=(ws(4)-ws(1))/(Nzon)
      
       do i2=(N_x2-Nzon)/2+1, (N_x2-Nzon)/2+Nzon +1

       x2n_array(i1,i2) =ws(1)+real(i2-1-(N_x2-Nzon)/2,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(1)+real(i2-(N_x2-Nzon)/2-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(4)-ws(1))/(slope_mesh1)
        !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)/slope_mesh1
                       
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(1)+real(i2-(N_x2-Nzon)/2-0.5_f64,f64)*dxx
      enddo
     enddo

     do i1=1,nc_eta1+1   

       dxx=(x2_max-ws(4))/((N_x2-Nzon)/2._f64)
      
       do i2= (N_x2-Nzon)/2+Nzon +1, N_x2+1

       x2n_array(i1,i2) =ws(4)+real(i2-1-((N_x2-Nzon)/2+Nzon),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2)  =ws(4)+real(i2-((N_x2-Nzon)/2+Nzon)-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/(slope_mesh2)
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)/slope_mesh2
        
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =ws(4)+real(i2-((N_x2-Nzon)/2+Nzon)-0.5_f64,f64)*dxx
      enddo
     enddo

        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
     


    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    !do i2=1,nc_eta2+1
      !do i1=1,nc_eta1+1
        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
      !enddo
    !enddo
    
   
     
    
  endif
! do i1=1,nc_eta1 +1  
!       do i2= 1, N_x2+1
!         print*, x1c_array(i1,i2),x2c_array(i1,i2),ws(1),ws(2),ws(3),ws(4)!x1c_array(i1,i2),x2c_array(i1,i2),integration_points(1,i1,i2),integration_points(2,i1,i2)  
!      enddo
!     enddo
!stop
  if(mesh_case==11)then

      slope_mesh1=1._f64/2._f64
      wk=2.82_f64
      ll=1._f64
       
      Ns(1)=int((-wk-ll-x2_min)/delta_x2)+1
      ws(1)=x2_min+(Ns(1)-1)*delta_x2
      
      ll=(int(ll/delta_x2)+1)*delta_x2
     ! print*,ll
     !stop
      !Ns(2)=int((-wk+ll-x2_min)/delta_x2)+1
      ws(2)=ws(1)+2*ll
      
     ! Ns(3)=int((wk-ll-x2_min)/delta_x2)
      !ws(3)=x2_min+(Ns(3)-1)*delta_x2
      Ns(4)=int((wk+ll-x2_min)/delta_x2)+1
      ws(4)=x2_max-(ws(1)-x2_min)!x2_min+(Ns(4)-1)*delta_x2
      ws(3)=ws(4)-2*ll
     !iiw2=int((-wk-1-x2_min)/delta_x2)
    ! print*,ws(1)-ws(2),ws(3)-ws(4),delta_x2
     !stop
     !stop
     Nzon= int((ws(2)-ws(1))/(slope_mesh1*delta_x2))

     !print*,Nzon,int((ws(4)-ws(3))/(slope_mesh1*delta_x2))
     !stop
     Nzon2=3!int((ws(3)-ws(2))/(delta_x2))

     if(mod(N_x2-2*Nzon-Nzon2,2)==0) then
      Nzon2=Nzon2
     else
      Nzon2=Nzon2+1
     endif 
     ! print*,"slop",slope_mesh1,Nzon,N_x2-2*Nzon
      !stop
      !if(
      !Nzon= int(1._f64/(slope_mesh1*delta_x2))
      !slope_mesh2=(x2_max-x2_min-(ws(4)-ws(1)))/((N_x2-Nzon)*delta_x2)
        slope_mesh2=(ws(1)-x2_min)/((N_x2-2*Nzon-Nzon2)*delta_x2/2._f64)

        slope_mesh1=(ws(2)-ws(1))/(Nzon*delta_x2)
        slope_mesh3=(ws(3)-ws(2))/(Nzon2*delta_x2)
       
       !slope_mesh2=(x2_max-x2_min-2)/((N_x2-2*Nzon)*delta_x2)
       !print*,slope_mesh1,slope_mesh2!, (ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64),Nzon,N_x2-Nzon
       !(ws(2)-ws(1))-(ws(4)-ws(3)),
       !stop
      
      !do i1=1,nc_eta1+1  
       !do i2=1,nc_eta2+1  
         !x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       !enddo
      !enddo
     
 
     do i1=1,nc_eta1+1   

       dxx=(ws(1)-x2_min)/((N_x2-2*Nzon-Nzon2)/2._f64)
      
       do i2=1, (N_x2-2*Nzon-Nzon2)/2+1

       x2n_array(i1,i2) =x2_min+real(i2-1,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*dxx
       
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/slope_mesh2
       
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =x2_min+(real(i2,f64)-0.5_f64)*dxx
       !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo
     
     

     do i1=1,nc_eta1+1   

       dxx=(ws(2)-ws(1))/(Nzon)
      
       do i2=(N_x2-2*Nzon-Nzon2)/2+1, (N_x2-2*Nzon-Nzon2)/2+Nzon +1

       x2n_array(i1,i2) =ws(1)+real(i2-1-(N_x2-2*Nzon-Nzon2)/2,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(1)+real(i2-(N_x2-2*Nzon-Nzon2)/2-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(2)-ws(1))/slope_mesh1
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(1)+real(i2-(N_x2-2*Nzon-Nzon2)/2-0.5_f64,f64)*dxx
       !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo

     do i1=1,nc_eta1+1   

       dxx=(ws(3)-ws(2))/(Nzon2)
      
       do i2= (N_x2-2*Nzon-Nzon2)/2+Nzon +1, (N_x2-2*Nzon-Nzon2)/2+Nzon +Nzon2 +1

       x2n_array(i1,i2) =ws(2)+real(i2-1-((N_x2-2*Nzon-Nzon2)/2+Nzon),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(2)+real(i2-((N_x2-2*Nzon-Nzon2)/2+Nzon)-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(3)-ws(2))/slope_mesh3
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(2)+real(i2-((N_x2-2*Nzon-Nzon2)/2+Nzon)-0.5_f64,f64)*dxx
       !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo 

       do i1=1,nc_eta1+1   

       !dxx=(ws(4)-ws(3))/(Nzon)
        dxx=(ws(2)-ws(1))/(Nzon)
       do i2=  (N_x2-2*Nzon-Nzon2)/2+Nzon +Nzon2 +1, (N_x2-2*Nzon-Nzon2)/2+2*Nzon +Nzon2 +1

       x2n_array(i1,i2) =ws(3)+real(i2-1-((N_x2-2*Nzon-Nzon2)/2+Nzon+Nzon2),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(3)+real(i2-((N_x2-2*Nzon-Nzon2)/2+Nzon+Nzon2)-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(2)-ws(1))/slope_mesh1
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(3)+real(i2-((N_x2-2*Nzon-Nzon2)/2+Nzon+Nzon2)-0.5_f64,f64)*dxx
      ! print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo 

     
     do i1=1,nc_eta1+1   

       !dxx=(x2_max-ws(4))/((N_x2-Nzon)/2._f64)
        dxx=(x2_max-ws(4))/((N_x2-2*Nzon-Nzon2)/2._f64)
       
       do i2= (N_x2-2*Nzon-Nzon2)/2+2*Nzon +Nzon2 +1, N_x2+1

       x2n_array(i1,i2) =ws(4)+real(i2-1-((N_x2-2*Nzon-Nzon2)/2+2*Nzon+Nzon2),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2)  =ws(4)+real(i2-((N_x2-2*Nzon-Nzon2)/2+2*Nzon+Nzon2)-0.5_f64,f64)*dxx
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)*slope_mesh2
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/slope_mesh2
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =ws(4)+real(i2-((N_x2-2*Nzon-Nzon2)/2+2*Nzon+Nzon2)-0.5_f64,f64)*dxx
        !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo
     !stop
        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
     


    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    !do i2=1,nc_eta2+1
      !do i1=1,nc_eta1+1
        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
      !enddo
    !enddo
    
   
     !print*,slope_mesh1,slope_mesh2,slope_mesh3,
   !print*, (2*(ws(1)-x2_min)/slope_mesh2+2*(ws(2)-ws(1))/slope_mesh1+(ws(3)-ws(2))/slope_mesh3)*(x1_max-x1_min), (x1_max-x1_min)*(x2_max-x2_min)!,2*sll_Pi
    
  endif
   

  
  
  
  
  
  
  end subroutine construct_bgk_mesh  
  
  subroutine compute_bgk_phi(L,Nx,mu,xi,tab_phi,tab_dphi)
    use sll_constants
    implicit none
    sll_real64,intent(in) ::L,mu,xi
    sll_int,intent(in) :: Nx
    sll_real64 :: xmin,xmax,x,dx,k2,phi(2)
    sll_int:: j
    sll_real64,dimension(:), pointer :: tab_phi,tab_dphi 
    xmin = 0._f64
    xmax = L/2._f64
    phi = 0._f64
    x   = xmin
    dx  = (xmax-xmin)/real(Nx,f64)
    k2  = mu*(1._f64-2._f64*xi+2._f64*phi(1))/(3._f64-2._f64*xi)*exp(-phi(1))
    tab_phi(1) = phi(1)
    tab_dphi(1)= phi(2)
    !write(10,* ) x, phi(1), phi(2),k2,v1,v2,1.d0/phi(2)
    do j=2, Nx+1
      x =xmin+ (j-1)*dx
      call Solve_Runge_K4(mu,xi,x, phi, dx, 2, Diff_Eq,0._f64,1._f64)
      k2=mu*(1._f64-2._f64*xi+2._f64*phi(1))/(3._f64-2._f64*xi)*exp(-phi(1))
      tab_phi(j) = phi(1)
      tab_dphi(j)= phi(2)
      ! write(10,* ) x, phi(1), phi(2),k2,v1,v2,1.d0/phi(2)
    enddo!j 
  
    !open(unit=10,file='phi.dat')
    !  write(10,*) Nx,L
    !  do j=1,Nx+1
    !    write(10,*) j,tab_phi(j),tab_dphi(j)
    !  enddo
    !close(10)
    
  end subroutine compute_bgk_phi
  
  subroutine construct_mesh_bgk(tab_phi,h_positions,tab_phi_positions,h_theta_positions,&
  Nx,Nen,Nen2,h_max,L,N_theta)
    sll_real64,dimension(:), pointer :: tab_phi,h_positions 
    sll_real64,dimension(:,:), pointer :: tab_phi_positions 
    sll_real64,dimension(:,:,:), pointer :: h_theta_positions 
    sll_int,intent(in)::Nx,Nen,Nen2,N_theta
    sll_int::jj,i,j,N_h,k
    sll_real64,intent(in) :: h_max,L
    sll_real64 :: delta_v,v,h,dx,xmin,xmax,x_factor,dtheta
    sll_real64 :: y_factor,total_length,tmp_loc,alpha,length
    
    
    N_h = Nen+Nen2
    
    xmin = 0._f64
    xmax= L/2._f64
    
    dx  = (xmax-xmin)/real(Nx,f64)
    
    !check compatibility between Nx and Nen
    jj= Nx/Nen
    if(jj==0)then
      print *,'bad compatibility between Nx=',Nx,' and Nen=',Nen
      print *,'Nx/Nen should be a non zero integer'
      stop
    endif
    if(jj*Nen/=Nx)then
      print *,'bad compatibility between Nx=',Nx,' and Nen=',Nen
      print *,'Nx/Nen should be a non zero integer'
      stop
    endif
     
    ! first define the h_positions for the grid
    jj=(Nx/Nen)
    do j=1,Nen+1
      h_positions(j)=tab_phi(1+(j-1)*jj)
    enddo


  
    delta_v = sqrt(2._f64*(h_positions(Nen+1)))/real(Nen,f64)
    do j=2,Nen
      v = real(j-1,f64)*delta_v
      !h_positions(j)=0.5_f64*v*v
    enddo
  
  
    do j=1,Nen2
      h_positions(Nen+1+j) = (sqrt(h_positions(Nen+1+j-1))+(sqrt(h_max)-sqrt(h_positions(Nen+1)))/real(Nen2,8))**2
    enddo

  
    print *,'#h_min error=',h_positions(1)
    print *,'#h_max error=',h_positions(N_h+1)-h_max
    h_positions(1)=0._f64
    h_positions(N_h+1)=h_max
  
  
    ! compute length for a given h
    do i=2,N_h+1
      h=h_positions(i)
      j=1
      do while((h-tab_phi(j)>=0._f64).and.(j<=Nx))
        tab_phi_positions(2,j) = sqrt(2._f64*(h-tab_phi(j)))
        tab_phi_positions(1,j) = xmin+real(j-1,f64)*dx
        j=j+1
      enddo
      if((j==Nx+1).and.(h-tab_phi(j)>=0._f64))then
        tab_phi_positions(1,j) = xmin+real(j-1,f64)*dx
        tab_phi_positions(2,j) = sqrt(2._f64*(h-tab_phi(j)))    
        j=j+1
      endif    
      jj=j-1
      !compute length
      y_factor = tab_phi_positions(1,jj)**2
      x_factor = tab_phi_positions(2,1)**2
      total_length = 0._f64
      do j=1,jj-1
        tmp_loc = y_factor*(tab_phi_positions(2,j+1)-tab_phi_positions(2,j))**2
        tmp_loc = tmp_loc+x_factor*(tab_phi_positions(1,j+1)-tab_phi_positions(1,j))**2
        total_length = total_length+sqrt(tmp_loc)    
      enddo
  
  
      if(total_length<=0._f64)then
        print *,'length is zero',total_length
        stop
      endif
      print *,'#length=',total_length
      !dtheta(i) =total_length/real(N_theta,f64)
      dtheta =total_length/real(N_theta,f64)
     !index_theta_positions(1,i) = 1
      length=0._f64
      j=1
      do k=2,N_theta
        do while(length<=real(k-1,8)*dtheta)
          tmp_loc = y_factor*(tab_phi_positions(2,j+1)-tab_phi_positions(2,j))**2
          tmp_loc = tmp_loc+x_factor*(tab_phi_positions(1,j+1)-tab_phi_positions(1,j))**2
          tmp_loc= sqrt(tmp_loc)
          length = length+tmp_loc
          j=j+1
        enddo  
        length = length-tmp_loc
        j=j-1
        !index_theta_positions(k,1) = j
        if((length>real(k-1,f64)*dtheta).or.(real(k-1,f64)*dtheta>=length+tmp_loc))then
          print *,'Problem of localization',j,length,real(k-1,f64)*dtheta,length+tmp_loc
          stop
        endif
        alpha= (real(k-1,f64)*dtheta-length)/tmp_loc
        if((alpha>=1._f64).or.(alpha<0._f64))then
          print *,'bad value of alpha=',alpha
        endif
      
        if(j>Nx)then
          print *,'j=',j,' Nx=',Nx
          print *,'length=',length,'Ntheta*dtheta=',real(N_theta,f64)*dtheta
          stop
        endif
        h_theta_positions(1,k,i) =  alpha*tab_phi_positions(1,j+1)+(1._f64-alpha)*tab_phi_positions(1,j)
        h_theta_positions(2,k,i) =  alpha*tab_phi_positions(2,j+1)+(1._f64-alpha)*tab_phi_positions(2,j)
        !print *,j,tmp,real(k-1,8)*dtheta(i),tmp+tmp_loc
      enddo
      !k=1
      h_theta_positions(1,1,i) =  tab_phi_positions(1,1)
      h_theta_positions(2,1,i) =  tab_phi_positions(2,1)
      !k=N_theta+1
      h_theta_positions(1,N_theta+1,i) =  tab_phi_positions(1,jj)
      h_theta_positions(2,N_theta+1,i) =  tab_phi_positions(2,jj)
    
    enddo

    open(unit=10,file='meshi.dat')
    !write(10,*) Nx,L
    do i= 1, N_h+1
      do k=1,N_theta+1
        write(10,* ) i,k, h_theta_positions(1,k,i), h_theta_positions(2,k,i)
      enddo
    enddo
    close(10)

    
    
  end subroutine construct_mesh_bgk
  
    

  

!*****************************************************************************************
!  Subroutine: To solve  the differential Equation with Runge Kutta 4 method             *
!*****************************************************************************************
  subroutine Solve_Runge_K4(mu,xi,X,Phi,dx,Neq,Diff_Eq,Y,C)
    use sll_constants
    implicit none
#ifdef __INTEL_COMPILER
    External Diff_Eq 
#endif
    sll_int   :: Neq
    sll_real64   :: x, phi(Neq), dx,mu,xi,y,c
    sll_real64     :: Step1(Neq),Step2(Neq),Step3(Neq), Step4(Neq), phi_temp(Neq)
    sll_int   :: i

    call Diff_Eq(mu,xi,x,phi,Step1,Neq,y,c)
    do i=1, Neq
      phi_temp(i) = phi(i) + Step1(i)*(dx/2._f64)
    enddo
    call Diff_Eq(mu,xi,x+dx/2._f64,phi_temp,Step2,Neq,y,c)
    do i=1, Neq
      phi_temp(i) = phi(i) + Step2(i)*(dx/2._f64)
    enddo
    call Diff_Eq(mu,xi,x+dx/2._f64,phi_temp,Step3,Neq,y,c)
    do i=1, Neq
      phi_temp(i) = phi(i) + Step3(i)*dx
    enddo
    call Diff_Eq(mu,xi,x+dx,phi_temp,Step4,Neq,y,c)
    do i=1, Neq
      phi(i)=phi(i)+(dx/6._f64)*(Step1(i)+2._f64*Step2(i)+2._f64*Step3(i) + Step4(i))
    enddo
    return
  end subroutine Solve_Runge_K4
!*****************************************************************************************
!  Subroutine: To define  the differential Equation                                      *
!*****************************************************************************************
  subroutine Diff_Eq(mu,xi,X,Phi,Dphi, Neq,Y,C)
    use sll_constants
    implicit none
    sll_real64   :: mu,xi,y,c
    sll_int  :: Neq
    sll_real64   :: x, phi(Neq), dphi(Neq)
    ! Differential Equations
    if((mu.ne.0).and.(xi.ne.0)) then
      dphi(1) =phi(2)
      dphi(2) = -mu*(3._f64-2._f64*xi+2._f64*phi(1))/(3._f64-2._f64*xi)*exp(-phi(1))+1._f64
    else
      dphi(1) =-C*phi(1)/Y
      dphi(2) =0._f64
    endif
    return
  end subroutine Diff_Eq




end module bgk_mesh_construction
