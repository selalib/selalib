module module_cg_curvi_function

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "selalib.h"

   use module_cg_curvi_structure
   use sll_boundary_condition_descriptors
    !--*Poisson*----
  
contains

!!!****************************************************************************************

!!***************************************************************************************
subroutine init_distribution_curvilinear(N_eta1,N_eta2, fcase,f,mesh_case,&
  &x1_array,x2_array,x1c,x2c,sigma_x1,sigma_x2)

implicit none


  sll_real64, dimension (:,:), intent(inout) :: f
  sll_real64, dimension (:,:), pointer       :: x1_array,x2_array
  sll_int32,  intent(in) :: N_eta1,N_eta2
  sll_int32,  intent(in) :: fcase,mesh_case
  sll_real64, intent(in) :: x1c,x2c,sigma_x1,sigma_x2
  sll_real64 :: x,y
  sll_int32  :: i, j



! test-function
    
 select case (fcase)
  case(1)

     do i=1,N_eta1+1
        do j=1,N_eta2+1
        f(i,j)=1.0_f64
        enddo
     enddo


  case(2)

     do i=1,N_eta1+1
           do j=1,N_eta2+1
              f(i,j) = cos(x2_array(i,j))
           end do
     end do
     

  case(3)

     do i=1,N_eta1+1
        do j=1,N_eta2+1
          x=(x1_array(i,j)-x1c)/sigma_x1
          f(i,j) = exp(-0.5*x**2) 
        end do
     end do

  case(4) 
   
     do i=1,N_eta1+1
        do j=1,N_eta2+1
          x=(x1_array(i,j)-x1c)/sigma_x1
          y=(x2_array(i,j)-x2c)/sigma_x2
          f(i,j) = exp(-0.5*(x**2+y**2))+1
        end do
     end do

  case(5)
     
     do i=1,N_eta1+1
       do j=1,N_eta2+1
        f(i,j) = 0._f64
        if ((i==(N_eta1+1)/2).and.(j==(N_eta2+1)/2)) then
          f(i,j) = 1._f64
        endif
       enddo
     enddo

   case(6)
      do i=1,N_eta1+1
       do j=1,N_eta2+1
        f(i,j) = 0._f64
        if ((i==(N_eta1+1)/2).and.(j==(N_eta2+1)/2)) then
          f(i,j) = 1._f64
        endif
       enddo
     enddo
    

  case default 

     print*,"f is not defined"
   
     stop

  end select 
  
end subroutine init_distribution_curvilinear





!!***********************************************************************

subroutine construct_mesh_transF(nc_eta1,nc_eta2,mesh_case,&
   &x1n_array,x2n_array,jac_array,&
   &geom_eta,alpha_mesh,geom_x,jac_matrix)

    implicit none
    sll_real64,dimension(:,:),   pointer :: x1n_array,x2n_array  
    sll_real64,dimension(:,:,:), pointer :: jac_matrix
    sll_real64,dimension(:,:),   pointer :: jac_array
    sll_int32, intent(in)  :: nc_eta1,nc_eta2,mesh_case
    sll_real64,intent(in)  :: geom_eta(2,2),alpha_mesh
    sll_real64,intent(out) :: geom_x(2,2)
    sll_int32  :: i1,i2,err
    sll_real64 :: x1_min,x1_max,x2_min,x2_max,delta_x1,delta_x2
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1,eta2
    sll_real64 :: eta1_n,eta2_n
    !sll_real64, intent(inout) :: x1c,x2c
    
    

    ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1))
    ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1))
    ALLOCATE(jac_matrix(4,nc_eta1+1, nc_eta2+1)); jac_matrix=0._f64
    ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1)); jac_array=0.0_f64

   

    
    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)

    
    delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
    delta_eta2 = (eta2_max-eta2_min)/real(nc_eta2,f64)
    
    ! cartesian mesh
    if(mesh_case==1)then
      x1_min=eta1_min
      x1_max=eta1_max
      x2_min=eta2_min
      x2_max=eta2_max
  
      
      delta_x1 = (x1_max-x1_min)/real(nc_eta1,f64)
      delta_x2 = (x2_max-x2_min)/real(nc_eta2,f64)
      do i2=1,nc_eta2+1
        do i1=1,nc_eta1+1
          x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
          x2n_array(i1,i2) = x2_min+real(i2-1,f64)*delta_x2
          jac_array(i1,i2) = 1.0_f64 !(x1_max-x1_min)*(x2_max-x2_min)

          jac_matrix(1,i1,i2)=1._f64
          jac_matrix(2,i1,i2)=0._f64
          jac_matrix(3,i1,i2)=0._f64
          jac_matrix(4,i1,i2)=1._f64
        enddo
      enddo 
    endif
  
    ! polar mesh
      if (mesh_case==2) then 
         x1_min=-eta1_max
         x1_max=eta1_max
         x2_min=-eta2_max
         x2_max=eta2_max
        do i1= 1, nc_eta1 + 1
            eta1 = eta1_min + real(i1-1,f64)*delta_eta1
           do i2= 1, nc_eta2 + 1
            eta2 = eta2_min + real(i2-1,f64)*delta_eta2
            x1n_array(i1,i2) = eta1*cos(eta2)
            x2n_array(i1,i2) = eta1*sin(eta2)
            jac_array(i1,i2) = eta1
            
            !--*Jacobien matrix for (r,theta)*---
            jac_matrix(1,i1,i2)=cos(eta2)         !dx/deta1
            jac_matrix(2,i1,i2)=-eta1*sin(eta2)   !dx/deta2
            jac_matrix(3,i1,i2)=sin(eta2)         !dy/deta1
            jac_matrix(4,i1,i2)=eta1*cos(eta2)    !dy/deta2

             !--*Jacobien matrix for (theta,r)*---
            jac_matrix(1,i1,i2)=-eta1*sin(eta2)   !dx/deta1
            jac_matrix(2,i1,i2)=cos(eta2)         !dx/deta2
            jac_matrix(3,i1,i2)=eta1*cos(eta2)    !dy/deta1
            jac_matrix(4,i1,i2)=sin(eta2)         !dy/deta2
            
            
           end do
        end do

      endif
      
      ! polar-like mesh
      if (mesh_case==3) then
        do i1= 1, nc_eta1 + 1
           eta1 = eta1_min + real(i1-1,f64)*delta_eta1
           do i2= 1, nc_eta2 + 1
              eta2 = eta2_min + real(i2-1,f64)*delta_eta2
              x1n_array(i1,i2) = eta1*eta1*cos(eta2)
              x2n_array(i1,i2) = eta1*eta1*sin(eta2)
              jac_array(i1,i2) = 2._f64*eta1*eta1*eta1

            jac_matrix(1,i1,i2)=2*eta1*cos(eta2)         !dx/deta1
            jac_matrix(2,i1,i2)=-eta1*eta1*sin(eta2)     !dx/deta2
            jac_matrix(3,i1,i2)=2*eta1*sin(eta2)         !dy/deta1
            jac_matrix(4,i1,i2)=eta1*eta1*cos(eta2)      !dy/deta2
           end do
        end do
      endif

     ! Collela mesh
    if(mesh_case==4)then
       x1_min = eta1_min
       x1_max = eta1_max
       x2_min = eta2_min
       x2_max = eta2_max
      do i1= 1, nc_eta1 + 1
         eta1 = eta1_min + real(i1-1,f64)*delta_eta1
         eta1_n=real(i1-1,f64)/real(nc_eta1,f64)
         do i2 = 1, nc_eta2 + 1
           eta2 = eta2_min + real(i2-1,f64)*delta_eta2
           eta2_n=real(i2-1,f64)/real(nc_eta2,f64)
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1_n) * sin(2*sll_pi*eta2_n)
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1_n) * sin(2*sll_pi*eta2_n)
           
           jac_array(i1,i2) = (1.0_f64+ alpha_mesh *(2._f64 *sll_pi/(eta1_max-eta1_min)) * &
           &cos (2*sll_pi*eta1_n) * sin (2*sll_pi*eta2_n)) * &
           & (1.0_f64 + alpha_mesh *2._f64 * sll_pi/(eta2_max-eta2_min)* & 
           & sin(2*sll_pi*eta1_n)*cos(2*sll_pi*eta2_n)) - &
           & alpha_mesh *2._f64 *sll_pi/(eta2_max-eta2_min) * sin (2*sll_pi*eta1_n) * cos (2*sll_pi*eta2_n) * &
           & alpha_mesh *2._f64 * sll_pi/(eta1_max-eta1_min) * cos (2*sll_pi*eta1_n) * sin (2*sll_pi*eta2_n)
        
          !x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
          !x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
          !jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)

            jac_matrix(1,i1,i2)= 1. + alpha_mesh *(2._f64 *sll_pi/(eta1_max-eta1_min)) * &
                                 &cos (2*sll_pi*eta1_n) * sin (2*sll_pi*eta2_n)       !dx/deta1
            jac_matrix(2,i1,i2)= alpha_mesh *2._f64 *sll_pi/(eta2_max-eta2_min)*sin(2*sll_pi*eta1_n)* & 
                                 & cos (2*sll_pi*eta2_n)                              !dx/deta2
            jac_matrix(3,i1,i2)= alpha_mesh *2._f64 * sll_pi/(eta1_max-eta1_min) * cos (2*sll_pi*eta1_n)* & 
                                 & sin (2*sll_pi*eta2_n)                              !dy/deta1
            jac_matrix(4,i1,i2)=1. + alpha_mesh *2._f64 * sll_pi/(eta2_max-eta2_min)* & 
                                 & sin(2*sll_pi*eta1_n)*cos(2*sll_pi*eta2_n)          !dy/deta2

            
            !jac_matrix(1,i1,i2)= jac_matrix(1,i1,i2)*(x1_max-x1_min)  !dx/deta1
            !jac_matrix(2,i1,i2)= jac_matrix(1,i1,i2)*(x1_max-x1_min)  !dx/deta2
            !jac_matrix(3,i1,i2)= jac_matrix(1,i1,i2)*(x2_max-x2_min)  !dy/deta1
            !jac_matrix(4,i1,i2)= jac_matrix(1,i1,i2)*(x2_max-x2_min)  !dy/deta2
          
      end do
    end do

  endif

  if(mesh_case==5)then
     x1_min = eta1_min
     x1_max = eta1_max
     x2_min = eta2_min
     x2_max = eta2_max
     do i1= 1, nc_eta1 + 1
       eta1 = eta1_min + real(i1-1,f64)*delta_eta1
       eta1_n=real(i1-1,f64)/real(nc_eta1,f64)
        do i2 = 1, nc_eta2 + 1
           eta2 = eta2_min + real(i2-1,f64)*delta_eta2
           eta2_n=real(i2-1,f64)/real(nc_eta2,f64)
           x1n_array(i1,i2) = eta1_n + alpha_mesh * sin(2*sll_pi*eta1_n) * sin(2*sll_pi*eta2_n)**2
           x2n_array(i1,i2) = eta2_n + alpha_mesh * sin(2*sll_pi*eta1_n) * sin(2*sll_pi*eta2_n)
          ! a refaire
          ! jac_array(i1,i2) = 1._f64+2._f64*sll_pi*alpha_mesh*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)&
          ! +2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1)&
          ! -2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)**2&
          ! -4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)*cos(2._f64*sll_pi*eta1)&
          !+4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)**3*cos(2._f64*sll_pi*eta1)
           x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
           !jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
        end do
      end do
  
    endif
        
  open(unit=900,file='xn_array.dat')  
    do i1=1,nc_eta1
      do i2=1,nc_eta2
        write(900,*) x1n_array(i1,i2),x2n_array(i1,i2)
      enddo  
    enddo
  close(900)

    geom_x(1,1)=x1_min
    geom_x(2,1)=x1_max
    geom_x(1,2)=x2_min
    geom_x(2,2)=x2_max

end subroutine construct_mesh_transF  
  







!!!************************************************************************
subroutine phi_analytique(phi_exact,plan,phi_case,x1n_array,x2n_array,a1,a2,x1c_r,x2c_r,jac_matrix)

    implicit none

    sll_real64, dimension(:,:), intent(inout), pointer :: phi_exact
    sll_real64, dimension(:,:), intent(in), pointer:: x1n_array,x2n_array 
    sll_real64,dimension(:,:,:),pointer,intent(in) :: jac_matrix
    type(sll_plan_adv_curvilinear), intent(inout), pointer :: plan
    sll_int32 :: i,j,N_eta1,N_eta2
    sll_int32, intent(in) :: phi_case
    sll_real64,intent(in) :: x1c_r,x2c_r
    sll_real64 :: a1,a2,eta1_min,eta2_min,delta_eta1,delta_eta2,eta1,eta2
    
   
  
    
    N_eta1 = plan%N_eta1
    N_eta2 = plan%N_eta2
    eta1_min=plan%eta1_min
    eta2_min=plan%eta2_min
    delta_eta1=plan%delta_eta1
    delta_eta2=plan%delta_eta2
    
    !-dphi/deta1=-(dx/deta1*dphi/dx+dy/deta1*dphi/dy)
    ! dphi/deta2=(dx/deta2*dphi/dx+dy/deta2*dphi/dy)

    select case(phi_case)
    case(1) ! translation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            phi_exact(i,j)=a1*x1n_array(i,j)+a2*x2n_array(i,j)
            plan%field(1,i,j)= -(a1*jac_matrix(1,i,j)+a2*jac_matrix(3,i,j))      
            plan%field(2,i,j)= a1*jac_matrix(2,i,j)+a2*jac_matrix(4,i,j) 
          end do
       end do

    case(2) !rotation
       do i=1,N_eta1+1 
          do j=1,N_eta2+1
            phi_exact(i,j)=((x1n_array(i,j)-x1c_r)**2+(x2n_array(i,j)-x2c_r)**2)*0.5_f64
            plan%field(1,i,j)= -((x1n_array(i,j)-x1c_r)*jac_matrix(1,i,j)+ (x2n_array(i,j)-x2c_r)*jac_matrix(3,i,j))
            plan%field(2,i,j)= (x1n_array(i,j)-x1c_r)*jac_matrix(2,i,j)+(x2n_array(i,j)-x2c_r)*jac_matrix(4,i,j)   
           !write(500,*) eta1_min+(i-1)*delta_eta1,phi_exact(i,j),(eta1_min+(i-1)*delta_eta1)*(eta1_min+(i-1)*delta_eta1)
          end do  
       end do
    case(3) !anisotropic rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            phi_exact(i,j)=0.5_f64*(a1*(x1n_array(i,j)-x1c_r)**2+a2*(x2n_array(i,j)-x2c_r)**2)
            plan%field(1,i,j)= -(a1*(x1n_array(i,j)-x1c_r)*jac_matrix(1,i,j)+ a2*(x2n_array(i,j)-x2c_r)*jac_matrix(3,i,j))
            plan%field(2,i,j)= a1*(x1n_array(i,j)-x1c_r)*jac_matrix(2,i,j)+a2*(x2n_array(i,j)-x2c_r)*jac_matrix(4,i,j)
          end do
       end do
   case default
    print*,'#no phi define'
    print*,'#phi_case =1'
   end select
    
  end subroutine phi_analytique

!!!************************************************************************
subroutine carac_analytique(phi_case,N_eta1,N_eta2,x1n_array,x2n_array,a1,a2,x1c,x2c,&
  x1_tab,x2_tab,t)

    implicit none

    !sll_real64, dimension(:,:), intent(inout), pointer :: phi_exact
    sll_real64, dimension(:,:), intent(in), pointer:: x1n_array,x2n_array 
    sll_real64, dimension(:,:), intent(out), pointer:: x1_tab,x2_tab 
    !type(sll_plan_adv_curvilinear), intent(inout), pointer :: plan
    !sll_real64,dimension(:,:),pointer,intent(in) :: jac_array
    sll_int32 :: i,j
    sll_int32, intent(in) :: phi_case,N_eta1,N_eta2
    sll_real64,intent(inout) :: x1c,x2c
    sll_real64,intent(in)  :: t
    sll_real64 :: a1,a2,t1
    
  
  
 t1=t  

    select case(phi_case)
    case(1) ! translation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            x1_tab(i,j) = x1n_array(i,j) -a2*t1
            x2_tab(i,j) = x2n_array(i,j) +a1*t1
          end do
       end do
    case(2) !rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            x1_tab(i,j) = x1c+cos(t1)*(x1n_array(i,j)-x1c)-sin(t1)*(x2n_array(i,j)-x2c)
            x2_tab(i,j) = x2c+cos(t1)*(x2n_array(i,j)-x2c)+sin(t1)*(x1n_array(i,j)-x1c)
          end do
       end do
    case(3) !anisotropic rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            x1_tab(i,j) = a1*x1c+a1*cos(t1)*(x1n_array(i,j)-x1c)-a2*sin(t1)*(x2n_array(i,j)-x2c)
            x2_tab(i,j) = a2*x2c+a2*cos(t1)*(x2n_array(i,j)-x2c)+a1*sin(t1)*(x1n_array(i,j)-x1c)            
          end do
       end do
   case default
    print*,'#no phi define'
    print*,'#phi_case =1'
   end select
   
   
    
  end subroutine carac_analytique


!!!*******************************************************************************************
!===================================
!  construction of plan_curvilinear_op
!===================================

  function new_curvilinear_op(geom_eta,N_eta1,N_eta2,grad_case,bc1_type,bc2_type) result(plan_sl)

    type(plan_curvilinear_op), pointer :: plan_sl
    sll_real64, intent(in)             :: geom_eta(2,2)
    sll_int32,  intent(in)             :: N_eta1, N_eta2,bc1_type,bc2_type
    sll_int32,  intent(in),   optional :: grad_case
    sll_int32  :: err
    sll_real64 :: eta1_min,eta2_min,eta1_max,eta2_max
    sll_real64 :: delta_eta1,delta_eta2
   
    ALLOCATE(plan_sl)

    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)


    if (.not. present(grad_case)) then
       plan_sl%grad_case=2
    else
       plan_sl%grad_case=grad_case
    end if

    
    if (bc1_type==sll_DIRICHLET) then 
        plan_sl%spl_phi => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
                                 & bc1_type,bc2_type,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    elseif (bc2_type==sll_DIRICHLET) then 
        plan_sl%spl_phi => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
                                 & bc1_type,bc2_type,const_slope_x2_min = 0._f64,const_slope_x2_max = 0._f64)
    else
        plan_sl%spl_phi => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
                                 & bc1_type,bc2_type)
    endif

 
  end function new_curvilinear_op

!===============================
!  deletion of plan_curvilinear_op
!===============================

  !>subroutine delete_plan_curvilinear_op(this)
  !>deletion of plan_curvilinear_op object
  subroutine delete_plan_curvilinear_op(this)

    implicit none

    type(plan_curvilinear_op), pointer :: this
    sll_int32 :: err
   
    call delete_spline_2d(this%spl_phi)
    DEALLOCATE(this)
    this=>null()

  end subroutine delete_plan_curvilinear_op

!===========================================
!  computation of gradient and divergence
!===========================================

  subroutine compute_grad_field(plan,phi,grad_phi, &
       & N_eta1, N_eta2,geom_eta,bc1_type,bc2_type)

    implicit none

    type(plan_curvilinear_op),    pointer       :: plan
    sll_real64, dimension(:,:),   intent(inout) :: phi
    sll_real64, dimension(:,:,:), intent(out)   :: grad_phi
    sll_int32,  intent(in),       optional      :: bc1_type,bc2_type
    sll_real64, intent(in)                      :: geom_eta(2,2)
    sll_int32,  intent(in)                      :: N_eta1, N_eta2
    
    sll_real64 :: eta1,eta2
    sll_real64 :: delta_eta1,delta_eta2
    sll_real64 :: eta1_min,eta2_min,eta1_max,eta2_max
    sll_int32 :: i,j

  
    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)

    delta_eta1= (eta1_max-eta1_min)/real(N_eta1,f64)
    delta_eta2= (eta2_max-eta2_min)/real(N_eta2,f64)

    
    if (plan%grad_case==1) then
       ! center formula for eta1 end eta2
       ! decenter on boundaries
       do i=2,N_eta1
          do j=2,N_eta2
             grad_phi(1,i,j)=(phi(i+1,j)-phi(i-1,j))/(2*delta_eta1)
             grad_phi(2,i,j)=(phi(i,j+1)-phi(i,j-1))/(2*delta_eta2)
          end do
       end do
       do j=1,N_eta2
          grad_phi(1,1,j)=(phi(2,j)-phi(1,j))/delta_eta1
          grad_phi(1,N_eta1+1,j)=(phi(N_eta1+1,j)-phi(N_eta1,j))/delta_eta1

          grad_phi(2,1,j)=(phi(1,j+1)-phi(1,j))/delta_eta2
          grad_phi(2,N_eta1+1,j)=(phi(N_eta1+1,j+1)-phi(N_eta1+1,j))/delta_eta2
       end do

       do i=1,N_eta1
          grad_phi(1,i,1)=(phi(i,2)-phi(i,1))/delta_eta2
          grad_phi(1,i,N_eta2+1)=(phi(i,N_eta1+1)-phi(i,N_eta1))/delta_eta2

          grad_phi(2,1,j)=(phi(1,j+1)-phi(1,j))/delta_eta2
          grad_phi(2,N_eta1+1,j)=(phi(i+1,N_eta2+1)-phi(i,N_eta2+1))/delta_eta2
       end do

    
    else if (plan%grad_case==2) then
       ! using splines for eta1 and eta2

       call compute_spline_2D(phi,plan%spl_phi)

       do j=1,N_eta2+1
          eta2=eta2_min+real(j-1,f64)*delta_eta2
          do i=1,N_eta1+1
             eta1=eta1_min+real(i-1,f64)*delta_eta1
             grad_phi(1,i,j)=interpolate_x1_derivative_2D(eta1,eta2,plan%spl_phi)
             grad_phi(2,i,j)=interpolate_x2_derivative_2D(eta1,eta2,plan%spl_phi)
          end do
       end do

    else
     
       print*,'no choosen way to compute grad'
    end if
     
    
     do j=1,N_eta2+1
        do i=1,N_eta1+1
           grad_phi(1,i,j)=-grad_phi(1,i,j)
        end do
     end do
     

    if (present(bc1_type)) then
        if(bc1_type==sll_PERIODIC) grad_phi(:,N_eta1+1,:)=grad_phi(:,1,:)
    endif
 
    if (present(bc2_type)) then
       if(bc2_type==sll_PERIODIC) grad_phi(:,:,N_eta2+1)=grad_phi(:,:,1)
    endif
  end subroutine compute_grad_field


!!!************************************************************************
subroutine correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,eta1,eta2)

   implicit none

    sll_real64, intent(in)    :: eta1_min,eta1_max,eta2_min,eta2_max
    sll_real64, intent(inout) :: eta1,eta2  
    sll_int32,  intent(in)    :: bc1_type,bc2_type
     
  
! --- Corrections on the BC ---
        if (bc1_type.eq.sll_DIRICHLET) then
          eta1 = min(max(eta1,eta1_min),eta1_max)
        endif
        if (bc2_type.eq.sll_DIRICHLET) then
          eta2 = min(max(eta2,eta2_min),eta2_max)
        endif
        if (bc1_type==sll_PERIODIC) then
          do while (eta1>eta1_max)
            eta1 = eta1-(eta1_max-eta1_min)
          enddo
          do while (eta1<eta1_min)
            eta1 = eta1+(eta1_max-eta1_min)
          enddo
        endif
        if (bc2_type==sll_PERIODIC) then
          do while (eta2>eta2_max)
            eta2 = eta2-(eta2_max-eta2_min)
          enddo
          do while (eta2<eta2_min)
            eta2 = eta2+(eta2_max-eta2_min)
          enddo
        endif
  end subroutine correction_BC    
!!!************************************************************************
!*****
subroutine plot_f1(f,x1,x2,iplot,N_eta1,N_eta2,delta_eta1,delta_eta2,eta1_min,eta2_min)!

  use sll_xdmf
  use sll_hdf5_io_serial
  sll_int32 :: file_id
  sll_int32 :: error
  sll_real64, dimension (:,:), allocatable :: f
  sll_real64, dimension(:,:), pointer :: x1 !allocatable :: x1
  sll_real64, dimension(:,:), pointer :: x2 !allocatable :: x2
  sll_int32, intent(in) :: N_eta1,N_eta2
  sll_int32 :: i, j
  sll_int32, intent(in) :: iplot
  character(len=4)      :: cplot
  sll_int32             :: nnodes_x1, nnodes_x2
  sll_real64 :: eta1,eta2,delta_eta1,delta_eta2,eta1_min,eta2_min

  nnodes_x1 = N_eta1+1
  nnodes_x2 = N_eta2+1

  if (iplot == 1) then
     call sll_hdf5_file_create("curvilinear_mesh-x1.h5",file_id,error)
     call sll_hdf5_write_array(file_id,x1,"/x1",error)
     call sll_hdf5_file_close(file_id, error)
     call sll_hdf5_file_create("curvilinear_mesh-x2.h5",file_id,error)
     call sll_hdf5_write_array(file_id,x2,"/x2",error)
     call sll_hdf5_file_close(file_id, error)
  end if

  call int2string(iplot,cplot)
  call sll_xdmf_open("f"//cplot//".xmf","curvilinear_mesh",nnodes_x1,nnodes_x2,file_id,error)
  call sll_xdmf_write_array("f"//cplot,f,"values",error,file_id,"Node")
  call sll_xdmf_close(file_id,error)!
!
 end subroutine plot_f1
!!!****************************************************************************************
subroutine init_distribution_cart(eta1_min,eta1_max,eta2_min,eta2_max,N_eta1,N_eta2, &  
           & fcase,f,mesh_case,mode,alpha,eta1_r,eta1_rm,sigma_x1,sigma_x2, &  
           & landau_alpha, landau_mode)

implicit none


  sll_real64, dimension (:,:), intent(inout) :: f
  sll_int32,  intent(in) :: N_eta1,N_eta2,mode
  sll_int32,  intent(in) :: fcase,mesh_case
  sll_real64, intent(in) :: eta1_min, eta1_max, eta2_min, eta2_max
  sll_real64 :: delta_eta1,delta_eta2
  sll_real64 :: eta1,eta2,eta1c,eta2c
  sll_int32  :: i, j
  sll_real64, optional :: eta1_r,eta1_rm,alpha
  sll_real64, optional :: landau_alpha, landau_mode
  sll_real64, optional :: sigma_x1,sigma_x2



  delta_eta1= (eta1_max-eta1_min)/real(N_eta1,f64)
  delta_eta2= (eta2_max-eta2_min)/real(N_eta2,f64)

  eta1c = 0.25_f64*eta1_max+0.75_f64*eta1_min
  eta2c = 0.25_f64*eta2_max+0.75_f64*eta2_min
!
!
!! test-function
!    
!
  select case (fcase)
   case(0)
     do j=1,N_eta2+1
        do i=1,N_eta1+1
             f(i,j)=1.0_f64
        enddo
     enddo

   case(1)

     print*,'#r_1 = ', eta1_r
     print*,'#r_2 = ',eta1_rm
     print*,'#alpha, mode = ',alpha,mode
     do j=1,N_eta2+1
       eta2=eta2_min+real(j-1,f64)*delta_eta2
        do i=1,N_eta1+1
       eta1=eta1_min+real(i-1,f64)*delta_eta1
        f(i,j)=0.0_f64
           if((eta1>eta1_r).and.(eta1<eta1_rm)) then
             f(i,j)=1.0_f64+ alpha*cos(mode*eta2)
           endif
        enddo
     enddo
    
  case(2)
     do j=1,N_eta2+1
           do i=1,N_eta1+1
              eta1=eta1_min+real(i-1,f64)*delta_eta1
              f(i,j) = cos(eta1)
           end do
     end do

!
   case(3) !khp testcase 
   
     do j=1,N_eta2+1
          eta2=eta2_min+real(j-1,f64)*delta_eta2
        do i=1,N_eta1+1
          eta1=eta1_min+real(i-1,f64)*delta_eta1
          f(i,j) = sin(eta2)+landau_alpha*cos(landau_mode*eta1)
        end do
     end do

   case(4) 
     
     do i=1,N_eta1+1
         eta1=eta1_min+real(i-1,f64)*delta_eta1
        do j=1,N_eta2+1
         eta2=eta2_min+real(j-1,f64)*delta_eta2
                 f(i,j) = exp(-((eta1-eta1c)/sigma_x1)**2)*exp(-((eta2-eta2c)/sigma_x2)**2)
        end do
     end do

   case(5)
     
     do i=1,N_eta1+1
       do j=1,N_eta2+1
        f(i,j) = 0._f64
        if ((i==(N_eta1+1)/2).and.(j==(N_eta2+1)/2)) then
          f(i,j) = 1._f64
        endif
       enddo
     enddo
!   
!  !case(6) 
!
!    ! open(25,file=f_file,action="read")
!    ! read(25,*)f
!    ! close(25)
!
   case default 

     print*,"f is not defined"
     stop
!
   end select 
!
end subroutine init_distribution_cart


!*****************************************************************************************
subroutine lire_appel(phi,f,N_eta1,N_eta2)
implicit none

  sll_real64, dimension (:,:), allocatable, intent(inout) :: f
  sll_real64, dimension (:,:), allocatable, intent(inout) :: phi
  sll_real64, dimension (:), allocatable :: vec_phi
  sll_int32 :: i, j,err
  sll_int32, intent(in) :: N_eta1,N_eta2

ALLOCATE(vec_phi((N_eta1+1)*(N_eta2+1)))
open (unit=28,file="file_rho.dat",action="write",status="old")
      do i=1,N_eta1+1
          do j=1,N_eta2+1 
             write(28,*) f(i,j)
          enddo
      enddo
   close(28)
   call system('./quasi_neutre')
   open(unit=29,file='values_U.dat',action="read",status="old")
    do i = 1,(N_eta1+1)*(N_eta2+1)
         read (unit=29,fmt=*) vec_phi(i)
    enddo 
    close(unit=29)
    !vec_phi=0.
    phi(1:(N_eta1+1),1:(N_eta2+1))=reshape(vec_phi,(/(N_eta1+1),(N_eta2+1)/))

DEALLOCATE(vec_phi) 
end subroutine lire_appel
!***************************************************************************
!subroutine Poisson_solver_curvilinear(eta1_min,eta1_max,eta2_min,eta2_max, &
!                                    & phi,init,BC_P, option,apr_U,point1,point2, &
!                                    & f,tab_b11,tab_b12,tab_b21,tab_b22, tab_c, &
!                                    & tab_dF11,tab_dF12,tab_dF21,tab_dF22,values_U, &
!                                    & spline_degree,N_eta1,N_eta2)
!
!sll_real64, dimension (:,:), allocatable :: phi
!sll_real64, dimension (:,:), allocatable :: f
!  integer   :: BC_P ! boundary conditions 0 if periodic-Drichlet and 1 if Dirichlet-Dirichlet  
!  real(8)   :: period
!  real(8), dimension(:), pointer :: apr_U
!  integer                        :: option ! 0 for analytique and 1 for interpolation 
!  type(initial) :: init
!  real(8), dimension(:,:), allocatable :: values_U
!  real(8), dimension(:), allocatable  :: point1
!  real(8), dimension(:), allocatable  :: point2
!  real(8), dimension(:,:), pointer,OPTIONAL  :: tab_b11,tab_b12,tab_b21,tab_b22
!  real(8), dimension(:,:), pointer,OPTIONAL  :: tab_c,tab_dF11,tab_dF12,tab_dF21,tab_dF22
!integer:: N_eta1,N_eta2,nbpt1,nbpt2,num_pts_eta1,num_pts_eta2 
!  integer :: spline_degree
!  real(8)   :: eta1_min, eta2_min, eta1_max, eta2_max
!  real(8)   :: conv_order_L1,conv_order_L2,conv_order_Linf
!
! num_pts_eta1 = N_eta2 +1
! num_pts_eta2 = N_eta1 +1
! nbpt1 = num_pts_eta1
! nbpt2 = num_pts_eta2
! period = eta2_max - eta2_min
! phi=0._f64
!
!
! call general_qn_solver(init,BC_P, option,apr_U,point1,point2, &
!       transpose(f), tab_b11,tab_b12,tab_b21,tab_b22,&
!        tab_dF11,tab_dF12,tab_dF21,tab_dF22,tab_c)
!   
!
! print *,'#general_qn_solver done'       
! values_U = 0.0_8
! if(BC_P==0) then
!  print*,'pass get values perdir',size(point1)
!        call get_values_function_perdir(init,period, apr_U, &  
!                                     & point1,nbpt1-1,point2,nbpt2,values_U )
! endif
!
! if(BC_P==1) then
!        call get_values_function_dirdir(init,apr_U, &
!                                     & point1,nbpt1,point2,nbpt2,values_U )
! endif
!
! if(BC_P==2) then
!         call get_values_function_perper(init,eta2_max-eta2_min,eta1_max-eta1_min, apr_U, &  
!                                     & point1,nbpt1,point2,nbpt2,values_U )
! endif
!
!
! phi=transpose(values_U)
!print*,'phi poisson ' , maxval(abs(f)),maxval(abs(phi))
!
! !---------------------------------------------------------------------------
!    !!!!!!!!!!!!!!!! CONVERGENCE TEST !!!!!!!!!!!!!!!!
!    !---------------------------------------------------------------------------
!    !call convergence_order(apr_U,init,BC_P,conv_order_L1,conv_order_L2,conv_order_Linf)
!
!   ! print*," ordre de convergence L1 =:)", conv_order_L1
!   ! print*," ordre de convergence L2 =:)", conv_order_L2
!   ! print*," ordre de convergence Linf =:)", conv_order_Linf
!    !open(unit=800,file='conv_order.dat',position="append")
!    !write(800,*) num_pts_eta1,conv_order_L1,conv_order_L2,conv_order_Linf  
!    !close(800)  
!  call free_tab()
!end subroutine Poisson_solver_curvilinear
!!**********************************************************************************

!!!************************************************************************


end module module_cg_curvi_function




!!!****************************************************************************************



!!***************************************************************************************
