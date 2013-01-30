module module_cg_curvi_function

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

   use module_cg_curvi_structure
   use numeric_constants
contains

!!!****************************************************************************************

subroutine init_distribution_cartesian(eta1_min,eta1_max,eta2_min,eta2_max,N_eta1,N_eta2, &  
           & delta_eta1,delta_eta2,fcase,f,mesh_case)

implicit none


  sll_real64, dimension (:,:), allocatable, intent(inout) :: f
  sll_int32 :: i, j
  sll_int32, intent(in) :: N_eta1,N_eta2
  sll_int32, intent(in):: fcase,mesh_case
  sll_real64, intent(in) :: eta1_min, eta1_max, eta2_min, eta2_max,delta_eta1,delta_eta2
  sll_real64 :: eta1,eta2,eta1c,eta2c


!  eta1c = 0.5_f64*(eta1_max+eta1_min)
!  eta2c = 0.5_f64*(eta2_max+eta2_min)

  eta1c = 0.25_f64*eta1_max+0.75_f64*eta1_min
  eta2c = 0.25_f64*eta2_max+0.75_f64*eta2_min


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
              eta2=eta2_min+real(j-1,f64)*delta_eta2
              f(i,j) = cos(eta2)
           end do
     end do
     

   case(3)

     do i=1,N_eta1+1
         eta1=eta1_min+real(i-1,f64)*delta_eta1
        do j=1,N_eta2+1
          f(i,j) = exp(-100._f64*(eta1-eta1c)**2) 
        end do
     end do

   case(4) 
     
     do i=1,N_eta1+1
         eta1=eta1_min+real(i-1,f64)*delta_eta1
        do j=1,N_eta2+1
         eta2=eta2_min+real(j-1,f64)*delta_eta2
              if ((mesh_case==1).or.(mesh_case==4)) then
                 f(i,j) = exp(-5_f64*(eta1-eta1c)**2)*exp(-5_f64*(eta2-eta2c)**2)
              endif
              if ((mesh_case==2).or.(mesh_case==3)) then
                f(i,j) = exp(-100._f64*(eta1-eta1c)**2)*exp(-30._f64*(eta2-eta2c)**2)
              endif
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
   
  !case(6) 

    ! open(25,file=f_file,action="read")
    ! read(25,*)f
    ! close(25)

   case default 

     print*,"f is not defined"
     print*,'see variable "fcase" in file selalib/prototype/src/simulation/CG_polar.F90'
     print*,'can not go any further'
     print*,'exiting...'
     stop

   end select 

end subroutine init_distribution_cartesian

!!***************************************************************************************
subroutine init_distribution_curvilinear(N_eta1,N_eta2, fcase,f,mesh_case,x1_tab,x2_tab,x1c,x2c)

implicit none


  sll_real64, dimension (:,:), allocatable, intent(inout) :: f
  sll_real64, dimension (:,:), pointer, intent(in) :: x1_tab,x2_tab
  sll_int32 :: i, j
  sll_int32, intent(in) :: N_eta1,N_eta2
  sll_int32, intent(in):: fcase,mesh_case
  sll_real64,intent(in) :: x1c,x2c


  !eta1c = 0.5_f64*(eta1_max+eta1_min)
  !eta2c = 0.5_f64*(eta2_max+eta2_min)

  


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
              f(i,j) = cos(x2_tab(i,j))
           end do
     end do
     

  case(3)

     do i=1,N_eta1+1
        do j=1,N_eta2+1
          f(i,j) = exp(-100._f64*(x1_tab(i,j)-x1c)**2) 
        end do
     end do

  case(4) 
     
     do i=1,N_eta1+1
        do j=1,N_eta2+1
              if ((mesh_case==1).or.(mesh_case==4)) then
                 f(i,j) = exp(-5_f64*(x1_tab(i,j)-x1c)**2)*exp(-5_f64*(x2_tab(i,j)-x2c)**2)
              endif
              if ((mesh_case==2).or.(mesh_case==3)) then
                f(i,j) = exp(-5._f64*(x1_tab(i,j)-x1c)**2)*exp(-5._f64*(x2_tab(i,j)-x2c)**2)
              endif
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
   

  case default 

     print*,"f is not defined"
   
     stop

  end select 
  
end subroutine init_distribution_curvilinear





!!***********************************************************************

subroutine construct_mesh_transF(nc_eta1,nc_eta2,mesh_case,&
   &x1n_array,x2n_array,jac_array,&
   &geom_eta,alpha_mesh,x1c,x2c,geom_x)

    implicit none
    sll_int32,intent(in) :: nc_eta1,nc_eta2,mesh_case
    sll_real64,intent(in) :: geom_eta(2,2),alpha_mesh
    sll_real64,intent(out) :: geom_x(2,2)
    sll_real64,dimension(:,:),pointer,intent(out) :: x1n_array,x2n_array  
    sll_real64,dimension(:,:),pointer,intent(out) :: jac_array
    sll_int32  :: i1,i2,err
    sll_real64 :: x1_min,x1_max,x2_min,x2_max,delta_x1,delta_x2
    sll_real64,intent(out) :: x1c,x2c
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1,eta2
    sll_real64 :: eta1_n,eta2_n
    
    

    SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), err); jac_array=0.0_f64

   

    
    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)

    
    delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
    delta_eta2 = (eta1_max-eta1_min)/real(nc_eta2,f64)
    
    ! cartesian mesh
    if(mesh_case==1)then
      x1_min=eta1_min
      x1_max=eta1_max
      x2_min=eta2_min
      x2_max=eta2_max
      x1c=(x1_min+x1_max)/2
      x2c=(x2_min+x2_max)/2
      
      delta_x1 = (x1_max-x1_min)/real(nc_eta1,f64)
      delta_x2 = (x2_max-x2_min)/real(nc_eta2,f64)
      do i2=1,nc_eta2+1
        do i1=1,nc_eta1+1
          x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
          x2n_array(i1,i2) = x2_min+real(i2-1,f64)*delta_x2
          jac_array(i1,i2) = 1.0_f64 !(x1_max-x1_min)*(x2_max-x2_min)
        enddo
      enddo 
    endif
  
    ! polar mesh
      if (mesh_case==2) then 
         x1_min=-eta1_max
         x1_max=eta1_max
         x2_min=-eta1_max
         x2_max=eta1_max
        x1c=((eta1_min+eta1_max)/2)*cos((eta2_min+eta2_max)/2)
        x2c=((eta1_min+eta1_max)/2)*sin((eta2_min+eta2_max)/2)
        do i1= 1, nc_eta1 + 1
            eta1 = eta1_min + real(i1-1,f64)*delta_eta1
           do i2= 1, nc_eta2 + 1
            eta2 = eta2_min + real(i2-1,f64)*delta_eta2
            x1n_array(i1,i2) = eta1*cos(eta2)
            x2n_array(i1,i2) = eta1*sin(eta2)
            jac_array(i1,i2) = eta1
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
           end do
        end do
      endif

     ! Collela mesh
    if(mesh_case==4)then
       x1_min = eta1_min
       x1_max = eta1_max
       x2_min = eta2_min
       x2_max = eta2_max
       x1c=(x1_min+x1_max)/2
       x2c=(x2_min+x2_max)/2
      do i1= 1, nc_eta1 + 1
         eta1 = eta1_min + real(i1-1,f64)*delta_eta1
         eta1_n=real(i1-1,f64)/real(nc_eta1,f64)
         do i2 = 1, nc_eta2 + 1
           eta2 = eta2_min + real(i2-1,f64)*delta_eta2
           eta2_n=real(i2-1,f64)/real(nc_eta2,f64)
           x1n_array(i1,i2) = eta1_n + alpha_mesh * sin(2*sll_pi*eta1_n) * sin(2*sll_pi*eta2_n)
           x2n_array(i1,i2) = eta2_n + alpha_mesh * sin(2*sll_pi*eta1_n) * sin(2*sll_pi*eta2_n)
           
           jac_array(i1,i2) = (1.0_f64 + alpha_mesh *(2._f64 *sll_pi/(eta1_max-eta1_min)) * &
           &cos (2*sll_pi*eta1_n) * sin (2*sll_pi*eta2_n)) * &
           & (1.0_f64 + alpha_mesh *2._f64 * sll_pi/(eta2_max-eta2_min) * sin (2*sll_pi*eta1_n) * cos (2*sll_pi*eta2_n)) - &
           & alpha_mesh *2._f64 *sll_pi/(eta2_max-eta2_min) * sin (2*sll_pi*eta1_n) * cos (2*sll_pi*eta2_n) * &
           & alpha_mesh *2._f64 * sll_pi/(eta1_max-eta1_min) * cos (2*sll_pi*eta1_n) * sin (2*sll_pi*eta2_n)
        
          x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
          x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
          !jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
      end do
    end do

  endif

  if(mesh_case==5)then
     x1_min = eta1_min! + alpha_mesh * sin(2*sll_pi*eta1_min) * sin(2*sll_pi*eta2_min)**2
     x1_max = eta1_max! + alpha_mesh * sin(2*sll_pi*eta1_max) * sin(2*sll_pi*eta2_max)**2
     x2_min = eta2_min! + alpha_mesh * sin(2*sll_pi*eta1_min) * sin(2*sll_pi*eta2_min)
     x2_max = eta2_max! + alpha_mesh * sin(2*sll_pi*eta1_max) * sin(2*sll_pi*eta2_max)
     x1c=(x1_min+x1_max)/2
     x2c=(x2_min+x2_max)/2
     do i1= 1, nc_eta1 + 1
       eta1 = eta1_min + real(i1-1,f64)*delta_eta1
        do i2 = 1, nc_eta2 + 1
           eta2 = eta2_min + real(i2-1,f64)*delta_eta2
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)**2
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           jac_array(i1,i2) = 1._f64+2._f64*sll_pi*alpha_mesh*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)&
           +2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1)&
           -2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)**2&
           -4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)*cos(2._f64*sll_pi*eta1)&
           +4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)**3*cos(2._f64*sll_pi*eta1)
           eta2 = eta2 + delta_eta2
           !x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           !x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
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
subroutine phi_analytique(phi_exact,plan,phi_case,x1n_array,x2n_array,a1,a2,x1c,x2c)

    implicit none

    sll_real64, dimension(:,:), intent(inout), pointer :: phi_exact
    sll_real64, dimension(:,:), intent(in), pointer:: x1n_array,x2n_array 
    type(sll_plan_adv_curvilinear), intent(inout), pointer :: plan
    sll_int32 :: i,j,N_eta1,N_eta2
    sll_int32, intent(in) :: phi_case
    sll_real64,intent(in) :: x1c,x2c
    sll_real64 :: a1,a2
   
  
    
    N_eta1 = plan%N_eta1
    N_eta2 = plan%N_eta2
     

    select case(phi_case)
    case(1) ! translation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            phi_exact(i,j)=a1*x1n_array(i,j)+a2*x2n_array(i,j)
            plan%field(1,i,j)= -a1
            plan%field(2,i,j)= -a2
          end do
       end do
    case(2) !rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            phi_exact(i,j)=((x1n_array(i,j)-x1c)**2+(x2n_array(i,j)-x2c)**2)*0.5
            plan%field(1,i,j)= -(x1n_array(i,j)-x1c)
            plan%field(2,i,j)= -(x2n_array(i,j)-x2c)
          end do
       end do
    case(3) !anisotropic rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            phi_exact(i,j)=0.5_f64*(a1*(x1n_array(i,j)-x1c)**2+a2*(x2n_array(i,j)-x2c)**2)
            plan%field(1,i,j)= -a1*(x1n_array(i,j)-x1c)
            plan%field(2,i,j)= -a2*(x2n_array(i,j)-x2c)
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
    sll_int32 :: i,j
    sll_int32, intent(in) :: phi_case,N_eta1,N_eta2
    sll_real64,intent(inout) :: x1c,x2c
    sll_real64,intent(in)  :: t
    sll_real64 :: a1,a2
    
   
  
    

    select case(phi_case)
    case(1) ! translation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            !phi_exact(i,j)=a1*x1n_array(i,j)+a2*x2n_array(i,j)
            !plan%field(1,i,j)= -a1
            !plan%field(2,i,j)= -a2
            x1_tab(i,j) = x1n_array(i,j) -a2*t 
            x2_tab(i,j) = x2n_array(i,j) +a1*t 
          end do
       end do
    case(2) !rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            !x1c=0._f64
            !x2c=0._f64
            !phi_exact(i,j)=((x1n_array(i,j)-x1c)**2+(x2n_array(i,j)-x2c)**2)*0.5
            !plan%field(1,i,j)= -(x1n_array(i,j)-x1c)
            !plan%field(2,i,j)= -(x2n_array(i,j)-x2c)
            x1_tab(i,j) = x1c+cos(t)*(x1n_array(i,j)-x1c)-sin(t)*(x2n_array(i,j)-x2c)
            x2_tab(i,j) = x2c+cos(t)*(x2n_array(i,j)-x2c)+sin(t)*(x1n_array(i,j)-x1c)
          end do
       end do
    case(3) !anisotropic rotation
       do i=1,N_eta1+1
          do j=1,N_eta2+1
            !phi_exact(i,j)=0.5_f64*a1*(x1n_array(i,j)-x1c)**2+0.5_f64*a2*(x2n_array(i,j)-x2c)**2
            !plan%field(1,i,j)= -a1*(x1n_array(i,j)-x1c)
            !plan%field(2,i,j)= -a2*(x2n_array(i,j)-x2c)
            x1_tab(i,j) = a1*x1c+a1*cos(t)*(x1n_array(i,j)-x1c)-a2*sin(t)*(x2n_array(i,j)-x2c)
            x2_tab(i,j) = a2*x2c+a2*cos(t)*(x2n_array(i,j)-x2c)+a1*sin(t)*(x1n_array(i,j)-x1c)            
          end do
       end do
   case default
    print*,'#no phi define'
    print*,'#phi_case =1'
   end select
    
  end subroutine carac_analytique





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

subroutine plot_f(iplot,N_eta1,N_eta2,delta_eta1,delta_eta2,eta1_min,eta2_min)

  use sll_xdmf
  use sll_hdf5_io
  sll_int32 :: file_id
  sll_int32 :: error
  sll_real64, dimension(:,:), allocatable :: x1
  sll_real64, dimension(:,:), allocatable :: x2
  sll_int32, intent(in) :: N_eta1,N_eta2
  sll_int32 :: i, j
  sll_int32, intent(in) :: iplot
  character(len=4)      :: cplot
  sll_int32             :: nnodes_x1, nnodes_x2
  sll_real64 :: eta1,eta2,delta_eta1,delta_eta2,eta1_min,eta2_min

  nnodes_x1 = N_eta1+1
  nnodes_x2 = N_eta2+1

  if (iplot == 1) then

     SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
     SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
     do j=1,nnodes_x2
     do i=1,nnodes_x1
        eta1 = eta1_min+real(i-1,f32)*delta_eta1
        eta2 = eta2_min+real(j-1,f32)*delta_eta2
        x1(i,j) = eta1*cos(eta2)
        x2(i,j) = eta1*sin(eta2)
     end do
     end do
     call sll_hdf5_file_create("curvilinear_mesh-x1.h5",file_id,error)
     call sll_hdf5_write_array(file_id,x1,"/x1",error)
     call sll_hdf5_file_close(file_id, error)
     call sll_hdf5_file_create("curvilinear_mesh-x2.h5",file_id,error)
     call sll_hdf5_write_array(file_id,x2,"/x2",error)
     call sll_hdf5_file_close(file_id, error)
     deallocate(x1)
     deallocate(x2)

  end if

  call int2string(iplot,cplot)
  call sll_xdmf_open("f"//cplot//".xmf","curvilinear_mesh",nnodes_x1,nnodes_x2,file_id,error)
  !call sll_xdmf_write_array("f"//cplot,f,"values",error,file_id,"Node")
  call sll_xdmf_close(file_id,error)

 end subroutine plot_f

end module module_cg_curvi_function
