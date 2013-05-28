module module_cg_curvi_function

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

   use module_cg_curvi_structure
   use sll_constants
   use sll_cubic_splines
   use sll_fft
   use sll_poisson_2d_polar
contains

!!!****************************************************************************************

!!***************************************************************************************
subroutine init_distribution_curvilinear(N_eta1,N_eta2, fcase,f,mesh_case,&
  &x1_array,x2_array,x1c,x2c,sigma_x1,sigma_x2)

implicit none


  sll_real64, dimension (:,:), allocatable, intent(inout) :: f
  sll_real64, dimension (:,:), pointer, intent(in) :: x1_array,x2_array
  sll_int32 :: i, j
  sll_int32, intent(in) :: N_eta1,N_eta2
  sll_int32, intent(in):: fcase,mesh_case
  sll_real64,intent(in) :: x1c,x2c,sigma_x1,sigma_x2
  sll_real64 :: x,y



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
    sll_int32,intent(in) :: nc_eta1,nc_eta2,mesh_case
    sll_real64,intent(in) :: geom_eta(2,2),alpha_mesh
    sll_real64,intent(out) :: geom_x(2,2)
    sll_real64,dimension(:,:),pointer,intent(out) :: x1n_array,x2n_array  
    sll_real64,dimension(:,:,:),pointer,intent(out) :: jac_matrix
    sll_real64,dimension(:,:),pointer,intent(out) :: jac_array
    sll_int32  :: i1,i2,err
    sll_real64 :: x1_min,x1_max,x2_min,x2_max,delta_x1,delta_x2
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1,eta2
    sll_real64 :: eta1_n,eta2_n
    !sll_real64, intent(inout) :: x1c,x2c
    
    

    SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(jac_matrix(4,nc_eta1+1, nc_eta2+1), err); jac_matrix=0._f64
    SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), err); jac_array=0.0_f64

   

    
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
         x2_min=-eta1_max
         x2_max=eta1_max
        do i1= 1, nc_eta1 + 1
            eta1 = eta1_min + real(i1-1,f64)*delta_eta1
           do i2= 1, nc_eta2 + 1
            eta2 = eta2_min + real(i2-1,f64)*delta_eta2
            x1n_array(i1,i2) = eta1*cos(eta2)
            x2n_array(i1,i2) = eta1*sin(eta2)
            jac_array(i1,i2) = eta1

            jac_matrix(1,i1,i2)=cos(eta2)         !dx/deta1
            jac_matrix(2,i1,i2)=-eta1*sin(eta2)   !dx/deta2
            jac_matrix(3,i1,i2)=sin(eta2)         !dy/deta1
            jac_matrix(4,i1,i2)=eta1*cos(eta2)    !dy/deta2
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
    sll_real64 ::delta_eta1,delta_eta2
    sll_real64, intent(in) :: geom_eta(2,2)
    sll_real64 :: eta1_min,eta2_min,eta1_max,eta2_max
    sll_int32, intent(in) :: N_eta1, N_eta2,bc1_type,bc2_type
    sll_int32, intent(in), optional :: grad_case
    sll_int32 :: err
   
    SLL_ALLOCATE(plan_sl,err)

    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)


    if (.not. present(grad_case)) then
       plan_sl%grad_case=2
    else
       plan_sl%grad_case=grad_case
    end if

    
    if (bc1_type==HERMITE_SPLINE) then 
        plan_sl%spl_phi => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
                                 & bc1_type,bc2_type,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
    elseif (bc2_type==HERMITE_SPLINE) then 
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
    SLL_DEALLOCATE(this,err)
    this=>null()

  end subroutine delete_plan_curvilinear_op

!===========================================
!  computation of gradient and divergence
!===========================================

  subroutine compute_grad_field(plan,phi,grad_phi, &
       & N_eta1, N_eta2,geom_eta,bc1_type,bc2_type,jac_matrix,x1n_array,x2n_array)

    implicit none

    type(plan_curvilinear_op), pointer              :: plan
    sll_real64, dimension(:,:), intent(inout) :: phi
    sll_real64, dimension(:,:,:), intent(out) :: grad_phi
    sll_real64 :: eta1,eta2
    sll_real64 ::delta_eta1,delta_eta2
    sll_real64, intent(in) :: geom_eta(2,2)
    sll_real64 :: eta1_min,eta2_min,eta1_max,eta2_max
    sll_int32, intent(in) :: N_eta1, N_eta2
    sll_int32, intent(in), optional :: bc1_type,bc2_type
    sll_int32 :: i,j


    sll_real64, dimension(:,:), intent(in), pointer,optional :: x1n_array,x2n_array 
   sll_real64,dimension(:,:,:),pointer,intent(in), optional :: jac_matrix
  
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

   ! else if (plan%grad_case==3) then
   !    ! using splines for x1 and x2

   !    call compute_spline_2D(phi,plan%spl_phi)

   !    do j=1,N_eta2+1
   !       do i=1,N_eta1+1
   !      grad_phi(1,i,j)=interpolate_x1_derivative_2D(x1n_array(i,j),x2n_array(i,j),plan%spl_phi)
   !      grad_phi(2,i,j)=interpolate_x2_derivative_2D(x1n_array(i,j),x2n_array(i,j),plan%spl_phi)
   !      grad_phi(1,i,j)=jac_matrix(1,i,j)*grad_phi(1,i,j)+jac_matrix(3,i,j)*grad_phi(2,i,j)
   !     grad_phi(2,i,j)=jac_matrix(2,i,j)*grad_phi(1,i,j)+jac_matrix(4,i,j)*grad_phi(2,i,j)       
   !       end do
   !    end do

    else
     
       print*,'no choosen way to compute grad'
    end if
     
    
     do j=1,N_eta2+1
        do i=1,N_eta1+1
           grad_phi(1,i,j)=-grad_phi(1,i,j)
        end do
     end do
     

    if (present(bc1_type)) then
        if(bc1_type==PERIODIC_SPLINE) grad_phi(:,N_eta1+1,:)=grad_phi(:,1,:)
    endif
 
    if (present(bc2_type)) then
       if(bc2_type==PERIODIC_SPLINE) grad_phi(:,:,N_eta2+1)=grad_phi(:,:,1)
    endif
  end subroutine compute_grad_field


!!!************************************************************************
subroutine correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,eta1,eta2)

   implicit none

    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,eta1,eta2  
    sll_int32  :: bc1_type,bc2_type
     
  
! --- Corrections on the BC ---
        if (bc1_type.eq.HERMITE_SPLINE) then
          eta1 = min(max(eta1,eta1_min),eta1_max)
        endif
        if (bc2_type.eq.HERMITE_SPLINE) then
          eta2 = min(max(eta2,eta2_min),eta2_max)
        endif
        if (bc1_type==PERIODIC_SPLINE) then
          do while (eta1>eta1_max)
            eta1 = eta1-(eta1_max-eta1_min)
          enddo
          do while (eta1<eta1_min)
            eta1 = eta1+(eta1_max-eta1_min)
          enddo
        endif
        if (bc2_type==PERIODIC_SPLINE) then
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
subroutine plot_f1(f,x1,x2,iplot,N_eta1,N_eta2,delta_eta1,delta_eta2,eta1_min,eta2_min)

  use sll_xdmf
  use sll_hdf5_io
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
  call sll_xdmf_close(file_id,error)

 end subroutine plot_f1
!!!****************************************************************************************
subroutine init_distribution_cartesian(eta1_min,eta1_max,eta2_min,eta2_max,N_eta1,N_eta2, &  
           & delta_eta1,delta_eta2,fcase,f,mesh_case,mode,alpha,eta1_r,eta1_rm,sigma_x1,sigma_x2)

implicit none


  sll_real64, dimension (:,:), allocatable, intent(inout) :: f
  sll_int32 :: i, j,mode
  sll_int32, intent(in) :: N_eta1,N_eta2
  sll_int32, intent(in):: fcase,mesh_case
  sll_real64, intent(in) :: eta1_min, eta1_max, eta2_min, eta2_max,delta_eta1,delta_eta2,sigma_x1,sigma_x2
  sll_real64 :: eta1,eta2,eta1c,eta2c
  sll_real64 :: eta1_r,eta1_rm,alpha

!!  eta1c = 0.5_f64*(eta1_max+eta1_min)
!!  eta2c = 0.5_f64*(eta2_max+eta2_min)
!
  eta1c = 0.25_f64*eta1_max+0.75_f64*eta1_min
  eta2c = 0.25_f64*eta2_max+0.75_f64*eta2_min
!
!
!! test-function
!    
!
  select case (fcase)
   case(1)

     print*,'eta1_r', eta1_r
     do j=1,N_eta2+1
       eta2=eta2_min+real(j-1,f64)*delta_eta2
        do i=1,N_eta1+1
       eta1=eta1_min+real(i-1,f64)*delta_eta1
        f(i,j)=0.0_f64
           if((eta1>eta1_r).and.(eta1<eta1_rm)) then
             f(i,j)=1.0_f64+ alpha*cos(mode*eta2)
           endif
         write(119,*) eta1,eta2,f(i,j)
        enddo
     enddo

!
  case(2)
     do j=1,N_eta2+1
           do i=1,N_eta1+1
              eta2=eta2_min+real(i-1,f64)*delta_eta2
              f(i,j) = cos(eta2)
           end do
     end do

!
  ! case(3)
!
    ! do i=1,N_eta1+1
    !     eta1=eta1_min+real(i-1,f64)*delta_eta1
    !    do j=1,N_eta2+1
    !      f(i,j) = exp(-100._f64*(eta1-eta1c)**2) 
    !    end do
    ! end do

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
end subroutine init_distribution_cartesian

!!!************************************************************************
subroutine diagnostic_1(f,plan_sl,phi_ref,int_r,bc,rmin,rmax,dr,dtheta,nr &
       & ,ntheta,nb_step,fcase,scheme,carac,grad,mode,&
       & l10,l20,e0,dt,alpha,r1,r2)

implicit none

 type(sll_SL_curvilinear), pointer :: plan_sl
sll_real64, dimension (:,:)  , allocatable :: phi_ref
sll_real64, dimension (:,:)  , allocatable :: f
sll_real64, dimension (:)  , allocatable :: int_r
sll_int32 :: bc(2)
sll_real64 :: r,theta,x,y
sll_real64, intent(in) :: r1,r2,rmax,rmin,dr,dtheta,alpha,dt
sll_real64 :: k1,k2,k3,c1,c2,c3 
sll_real64 :: k1_mode,k2_mode,k3_mode,c1_mode,c2_mode,c3_mode
sll_int32 ::  NEUMANN,NEUMANN_MODE0,DIRICHLET
sll_real64 :: w0,w, l1, l2, e, re, im,tmp,temps
sll_real64, intent(out) :: l10,l20,e0
sll_int32 :: mode,nb_step,i,j
sll_real64  :: temps_mode,err_loc
sll_int32, intent(in) :: fcase, scheme, carac, grad,nr, ntheta

 !mode=1.5
DIRICHLET=1
NEUMANN=2
NEUMANN_MODE0=3

  if(bc(1)==DIRICHLET)then
    k1 = (r1**2-r2**2+2.0_f64*r1**2*log(rmax/r1) + &
      2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
    k2 = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
      2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
    k3 = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
      2.0_f64*r2**2*log(r2/rmin))/(4.0_f64*log(rmin/rmax))
    c1 = (2.0_f64*r1**2*log(rmax/r1)+2.0_f64*r2**2*log(r2/rmax) + &
      r1**2-r2**2)*log(rmin)/(-4.0_f64*log(rmin/rmax))
    c2 = (2.0_f64*r2**2*log(rmin)*log(r2/rmax) + &
      2.0_f64*r1**2*log(rmax)*log(rmin/r1)+r1**2*log(rmax) - &
      r2**2*log(rmin))/(-4.0_f64*log(rmin/rmax))
    c3 = (r1**2-r2**2+2.0_f64*r2**2*log(r2/rmin) + &
      2.0_f64*r1**2*log(rmin/r1))*log(rmax) / &
      (-4.0_f64*log(rmin/rmax))
  endif
  if((bc(1)==NEUMANN).or.(bc(1)==NEUMANN_MODE0))then
    k1 = 0._f64
    k2 = r1**2/2._f64
    k3 = r1**2/2._f64-r2**2/2._f64
    c1 = r1**2*log(r1)/2._f64+r2**2/4._f64-r2**2*log(r2)/2._f64
    c1 = c1-r1**2*log(rmax)/2._f64+r2**2*log(rmax)/2._f64 - &
      r1**2/4._f64
    c2 = r2**2/4._f64-r2**2*log(r2)/2._f64-r1**2 * &
      log(rmax)/2._f64+r2**2*log(rmax)/2._f64
    c3 = -log(rmax)*(r1**2-r2**2)/2._f64
  endif

  if (mode==0) then 
    if(bc(1)==DIRICHLET)then
      k1_mode = (r1**2-r2**2+2.0_f64*r1**2*log(rmax/r1) + &
        2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
      k2_mode = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
        2.0_f64*r2**2*log(r2/rmax))/(4.0_f64*log(rmin/rmax))
      k3_mode = (r1**2-r2**2+2.0_f64*r1**2*log(rmin/r1) + &
        2.0_f64*r2**2*log(r2/rmin))/(4.0_f64*log(rmin/rmax))
      c1_mode = (2.0_f64*r1**2*log(rmax/r1) + &
        2.0_f64*r2**2*log(r2/rmax)+r1**2-r2**2)*log(rmin) / &
        (-4.0_f64*log(rmin/rmax))
      c2_mode = (2.0_f64*r2**2*log(rmin)*log(r2/rmax) + &
        2.0_f64*r1**2*log(rmax)*log(rmin/r1) + &
        r1**2*log(rmax)-r2**2*log(rmin))/(-4.0_f64*log(rmin/rmax))
      c3_mode = (r1**2-r2**2+2.0_f64*r2**2*log(r2/rmin) + &
        2.0_f64*r1**2*log(rmin/r1))*log(rmax) / &
        (-4.0_f64*log(rmin/rmax))
    endif
    if((bc(1)==NEUMANN).or.(bc(1)==NEUMANN_MODE0))then
      k1_mode = 0._f64
      k2_mode = r1**2/2._f64
      k3_mode = r1**2/2._f64-r2**2/2._f64
      c1_mode = r1**2*log(r1)/2._f64+r2**2/4._f64-r2**2 * &
        log(r2)/2._f64
      c1_mode = c1_mode-r1**2*log(rmax)/2._f64+r2**2 * &
        log(rmax)/2._f64-r1**2/4._f64
      c2_mode = r2**2/4._f64-r2**2*log(r2)/2._f64-r1**2 * &
        log(rmax)/2._f64+r2**2*log(rmax)/2._f64
      c3_mode = -log(rmax)*(r1**2-r2**2)/2._f64
    endif
  endif

  if(mode==1)then
    if((bc(1)==NEUMANN))then
      k1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)/(6._f64*(rmin**2+rmax**2))
      k2_mode = (3._f64*r2*rmax**2-r2**3 + &
        3*rmin**2*r1+r1**3)/(6._f64*(rmin**2+rmax**2))
      k3_mode = (-3._f64*r2*rmin**2-r2**3 + &
        3*rmin**2*r1+r1**3)/(6._f64*(rmin**2+rmax**2))
      c1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)*rmin**2/(6._f64*(rmin**2+rmax**2))
      c2_mode = -(3._f64*r1*rmin**2*rmax**2 + &
        r2**3*rmin**2-3*rmin**2*rmax**2*r2+rmax**2*r1**3) / &
        (6._f64 *(rmin**2+rmax**2))
      c3_mode = -(-3._f64*r2*rmin**2-r2**3 + &
        3*rmin**2*r1+r1**3)*rmax**2/(6._f64*(rmin**2+rmax**2))
    endif

    if((bc(1)==DIRICHLET).or.(bc(1)==NEUMANN_MODE0))then
      k1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)/(6._f64*(-rmin**2+rmax**2))
      k2_mode = (3._f64*r2*rmax**2-r2**3 - &
        3*rmin**2*r1+r1**3)/(6._f64*(-rmin**2+rmax**2))
      k3_mode = (3._f64*r2*rmin**2-r2**3 - &
        3*rmin**2*r1+r1**3)/(6._f64*(-rmin**2+rmax**2))
      c1_mode = (-3._f64*r1*rmax**2-r2**3 + &
        3*rmax**2*r2+r1**3)*rmin**2/(6._f64*(rmin**2-rmax**2))
      c2_mode = (-3._f64*r1*rmin**2*rmax**2 - &
        r2**3*rmin**2+3*rmin**2*rmax**2*r2+rmax**2*r1**3) / &
        (6._f64 *(rmin**2-rmax**2))
      c3_mode = (-3._f64*r1*rmin**2-r2**3 + &
        3*rmin**2*r2+r1**3)*rmax**2/(6._f64*(rmin**2-rmax**2))
    endif
  endif

  if(mode==3)then
    if((bc(1)==DIRICHLET).or.(bc(1)==NEUMANN_MODE0))then
      k1_mode = (r1*r2*(r2**5-r1**5) - &
        5._f64*rmax**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
      k2_mode = (r1*r2*(r2**5-r1**5) - &
        5._f64*(rmin**6*r2-rmax**6*r1))/(30._f64*r2*r1 * &
        (rmin**6-rmax**6))    
      k3_mode = (-r1*r2*(r1**5-r2**5) - &
        5._f64*rmin**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
      c1_mode = (-r1*r2*(r2**5-r1**5) + &
        5._f64*rmax**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
      c2_mode = (-r1*r2*(rmin**6*r2**5-rmax**6*r1**5) + &
        5._f64*(rmin*rmax)**6*(r2-r1)) / &
        (30._f64*r2*r1*(rmin**6-rmax**6))
      c3_mode = rmax**6*(-r1*r2*(r2**5-r1**5) + &
        5._f64*rmin**6*(r2-r1))/(30._f64*r2*r1*(rmin**6-rmax**6))
    endif
  endif

  if(mode==7)then
    if((bc(1)==DIRICHLET).or.(bc(1)==NEUMANN_MODE0))then
      k1_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r1**5*rmax**14 + &
        9._f64*r2**5*rmax**14
      k1_mode = k1/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      k2_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r2**5*rmin**14 - &
        9._f64*r1**5*rmax**14
      k2_mode = k2/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      k3_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r2**5*rmin**14 + &
        9._f64*r1**5*rmin**14
      k3_mode = k3/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      c1_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5-9._f64*r1**5*rmax**14 + &
        9._f64*r2**5*rmax**14
      c1_mode = rmin**14*c1 / &
        (630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      c2_mode = 5._f64*rmin**14*r1**5*r2**14 + &
        9._f64*rmax**14*rmin**14*r1**5
      c2_mode = c2_mode-9._f64*rmax**14*r2**5*rmin**14 + &
        5._f64*rmax**14*r1**14*r2**5
      c2_mode = -c2_mode/(630._f64*r1**5*r2**5*(rmin**14+rmax**14))
      c3_mode = -5._f64*r1**5*r2**14 + &
        5._f64*r1**14*r2**5+9._f64*r1**5*rmin**14 - &
        9._f64*r2**5*rmin**14
      c3_mode = -rmax**14*c3_mode / &
        (630._f64*r1**5*r2**5*(rmin**14+rmax**14))
    endif
  endif

  tmp = 0._f64
  l1  = 0.0_f64
  l2  = 0.0_f64

  open (unit=20,file='CGinit.dat')
  do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    if(mode==0)then
      if (r<r1) then
        temps_mode = k1_mode*log(r)+c1_mode
      else if (r>r2) then
        temps_mode = k3_mode*log(r)+c3_mode
      else
        temps_mode = k2_mode*log(r)+c2_mode-r**2/4.0_f64
      end if
    end if

    if((mode>=3).or.(mode==1))then
      if (r<r1) then
        temps_mode = k1_mode*r**mode+c1_mode/r**(mode)
      else if (r>r2) then
        temps_mode = k3_mode*r**mode+c3_mode/r**mode
      else
        temps_mode = k2_mode*r**mode+c2_mode/r**mode + &
          r**2/(mode**2-4._f64)
      end if
    end if


    if (r<r1) then
      temps = k1*log(r)+c1
    else if (r>r2) then
      temps = k3*log(r)+c3
    else
      temps = k2*log(r)+c2-r**2/4.0_f64
    end if

    temps_mode = temps_mode*alpha
    do j = 1,ntheta+1
      theta   = real(j-1,f64)*dtheta
      x       = r*cos(theta)
      y       = r*sin(theta)
      err_loc = abs(plan_sl%phi(i,j) - &
        temps_mode*cos(mode*theta)-temps)
      phi_ref(i,j) = temps_mode*cos(mode*theta)+temps

      write(20,*) r, theta, x, y, plan_sl%phi(i,j), &
        temps+temps_mode*cos(mode*theta)
      tmp = max(tmp,err_loc)
      if (i==1 .or. i==nr+1) then
        l1 = l1+err_loc*r/2.0_f64
        l2 = l2+err_loc**2*r/2.0_f64
      else
        l1 = l1+err_loc*r
        l2 = l2+err_loc**2*r
      end if
    end do
    write(20,*)' '
  end do
  close(20)
  l1 = l1*dr*dtheta
  l2 = sqrt(l2*dr*dtheta)
  print*,"#error for phi in initialization", &
    nr,dr,tmp/(1._f64+abs(alpha)), &
    l1/(1._f64+abs(alpha)),l2/(1._f64+abs(alpha))

  open(unit=23,file='diagMF.dat')
  write(23,*)'#fcase',fcase,'scheme',scheme, &
    'mode',mode,'grad',grad,'carac',carac
  write(23,*)'#nr',nr,'ntheta',ntheta,'alpha',alpha
  write(23,*)'#nb_step = ',nb_step,'  dt = ',dt
  write(23,*)'#   t   //   w   //   l1 rel  //   l2  rel //   e   //   re   //   im'

  do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)/r
  end do

  w0    = 0.0_f64
  l10   = 0.0_f64
  l20   = 0.0_f64
  e0    = 0.0_f64
  int_r = 0.0_f64
  do j = 1,ntheta
    w0  = w0+(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
    l10 = l10+abs(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
    l20 = l20+(f(1,j)/2.0_f64)**2*rmin+(f(nr+1,j)/2.0_f64)**2*rmax
    e0  = e0+rmin*(plan_sl%adv%field(1,1,j))**2/2.0_f64 + &
      rmax*(plan_sl%adv%field(1,nr+1,j))**2/2.0_f64 + &
      rmin*(plan_sl%adv%field(2,1,j))**2/2.0_f64 + &
      rmax*(plan_sl%adv%field(2,nr+1,j))**2/2.0_f64
    int_r(j) = (f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
    do i = 2,nr
      r   = rmin+real(i-1,f64)*dr
      w0  = w0+r*f(i,j)
      l10 = l10+r*abs(f(i,j))
      l20 = l20+r*f(i,j)**2
      e0  = e0+r*(plan_sl%adv%field(1,i,j)**2 + &
        plan_sl%adv%field(2,i,j)**2)
      int_r(j) = int_r(j)+f(i,j)*r
    end do
  end do

  w0    = w0*dr*dtheta
  l10   = l10*dr*dtheta
  l20   = sqrt(l20*dr*dtheta)
  e0    = e0*dr*dtheta/2.0_f64
  int_r = int_r*dr
  call fft_apply_plan(plan_sl%poisson%pfwd,int_r,int_r)
  write(23,*)'#t=0',w0,l10,l20,e0
  write(23,*)0.0_f64,w0,1.0_f64,1.0_f64,0.0_f64,e0, &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
    real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7)), &
    aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7))

 do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)*r
 end do
end subroutine diagnostic_1

!!******************
subroutine diagnostic_2(f,plan_sl,phi_ref,int_r,rmin,rmax,dr,dtheta,nr &
       & ,ntheta,step,l10,l20,e0,mode,dt,alpha)
implicit none
type(sll_SL_curvilinear), pointer :: plan_sl
sll_real64, dimension (:,:)  , allocatable :: phi_ref
sll_real64, dimension (:,:)  , allocatable :: f
sll_real64, dimension (:)  , allocatable :: int_r
sll_real64 :: r,theta
sll_real64 :: rmax,rmin
sll_real64 :: w, l10, l1, l20, l2, e, e0, re
sll_int32 :: mode,nb_step
sll_real64 :: dr,dtheta,alpha,dt
sll_int32  :: nr, ntheta,i,j,step


  do i = 1,nr+1
      r = rmin+real(i-1,f64)*dr
      plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)/r
    end do
    w     = 0.0_f64
    l1    = 0.0_f64
    l2    = 0.0_f64
    e     = 0.0_f64
    int_r = 0.0_f64
    do j = 1,ntheta
      w  = w+(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
      l1 = l1+abs(f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
      l2 = l2+(f(1,j)/2.0_f64)**2*rmin+(f(nr+1,j)/2.0_f64)**2*rmax
      e  = e+rmin*(plan_sl%adv%field(1,1,j))**2/2.0_f64 + &
        rmax*(plan_sl%adv%field(1,nr+1,j))**2/2.0_f64 + &
        rmin*(plan_sl%adv%field(2,1,j))**2/2.0_f64 + &
        rmax*(plan_sl%adv%field(2,nr+1,j))**2/2.0_f64
      int_r(j) = (f(1,j)*rmin+f(nr+1,j)*rmax)/2.0_f64
      int_r(j) = (plan_sl%phi(1,j)*rmin + &
        plan_sl%phi(nr+1,j)*rmax)/2.0_f64
      do i = 2,nr
        r  = rmin+real(i-1,f64)*dr
        w  = w+r*f(i,j)
        l1 = l1+r*abs(f(i,j))
        l2 = l2+r*f(i,j)**2
        e  = e+r*(plan_sl%adv%field(1,i,j)**2 + &
          plan_sl%adv%field(2,i,j)**2)
        int_r(j) = int_r(j)+plan_sl%phi(i,j)*r
      end do
    end do
    w     = w*dr*dtheta
    l1    = l1*dr*dtheta
    l2    = sqrt(l2*dr*dtheta)
    e     = e*dr*dtheta/2.0_f64
    int_r = int_r*dr
    call fft_apply_plan(plan_sl%poisson%pfwd,int_r,int_r)
    write(23,*) dt*real(step,f64),w,l1/l10,l2/l20,e-e0,e, &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &   !$7
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,0)), &  !$8
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &   !$9
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,1)), &  !$10
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &   !$11
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,2)), &  !$12
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &   !$13
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,3)), &  !$14
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &   !$15
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,4)), &  !$16
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &   !$17
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,7)), &  !$18
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), & 
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-1)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-2)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-3)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-4)), &
      real(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7)), &
      aimag(fft_get_mode(plan_sl%poisson%pfwd,int_r,ntheta-7))

do i = 1,nr+1
    r = rmin+real(i-1,f64)*dr
    plan_sl%adv%field(2,i,:) = plan_sl%adv%field(2,i,:)*r
  end do
end subroutine diagnostic_2
!*****************************************************************************************
subroutine lire_appel(phi,f,N_eta1,N_eta2)
implicit none

  sll_real64, dimension (:,:), allocatable, intent(inout) :: f
  sll_real64, dimension (:,:), allocatable, intent(inout) :: phi
  sll_real64, dimension (:), allocatable :: vec_phi
  sll_int32 :: i, j,err
  sll_int32, intent(in) :: N_eta1,N_eta2

SLL_ALLOCATE(vec_phi((N_eta1+1)*(N_eta2+1)),err)
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

SLL_DEALLOCATE_ARRAY(vec_phi,err) 
end subroutine lire_appel
!***************************************************************************


end module module_cg_curvi_function




!!!****************************************************************************************



!!***************************************************************************************
