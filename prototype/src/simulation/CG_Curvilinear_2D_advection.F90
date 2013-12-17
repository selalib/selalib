module curvilinear_2D_advection
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "selalib.h"

 
  use module_cg_curvi_function
  use module_cg_curvi_structure
  !--*Poisson*----
  !use sll_general_coordinate_qn_solver
   use sll_mudpack_colella

  implicit none

  

contains

!===================================
!  creation of sll_plan_adv_curvilinear
!===================================
  

  function new_plan_adv_curvilinear(geom_eta,delta_eta1,&  
 &delta_eta2,dt,N_eta1,N_eta2,carac_case,bc1_type,bc2_type) result(plan_sl)

    type(sll_plan_adv_curvilinear), pointer :: plan_sl
    sll_real64, intent(in) :: delta_eta1,delta_eta2,dt
    sll_real64, intent(in) :: geom_eta(2,2)
    sll_int32,  intent(in) :: N_eta1,N_eta2
    sll_int32,  intent(in) :: carac_case,bc1_type,bc2_type
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max
    

    sll_int32 ::err

    ALLOCATE(plan_sl)
    ALLOCATE(plan_sl%field(2,N_eta1+1,N_eta2+1))
    
    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)
    
    plan_sl%field=0.0_f64
    plan_sl%eta1_min=eta1_min
    plan_sl%eta1_max=eta1_max
    plan_sl%eta2_min=eta2_min
    plan_sl%eta2_max=eta2_max
    plan_sl%delta_eta1=delta_eta1
    plan_sl%delta_eta2=delta_eta2
    plan_sl%dt=dt
    plan_sl%N_eta1=N_eta1
    plan_sl%N_eta2=N_eta2
    plan_sl%carac_case=carac_case
    plan_sl%bc1_type=bc1_type
    plan_sl%bc2_type=bc2_type

 if (bc1_type==sll_DIRICHLET) then 
  plan_sl%spl_f => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
         & bc1_type,bc2_type,const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
 elseif (bc2_type==sll_DIRICHLET) then 
  plan_sl%spl_f => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
         & bc1_type,bc2_type,const_slope_x2_min = 0._f64,const_slope_x2_max = 0._f64)
 else
 plan_sl%spl_f => new_spline_2D(N_eta1+1,N_eta2+1,eta1_min,eta1_max,eta2_min,eta2_max, &
         & bc1_type,bc2_type)
 endif
  end function new_plan_adv_curvilinear

!===================================
!  deletion of sll_plan_adv_curvilinear
!===================================

  ! delete_plan_adv_curvilinear(plan_sl)
  ! deletion of sll_plan_adv_curvilinear object
  subroutine delete_plan_adv_curvilinear(plan_sl)

    implicit none

    type(sll_plan_adv_curvilinear), intent(inout), pointer :: plan_sl

    sll_int32 :: err

    if (associated(plan_sl)) then
       DEALLOCATE(plan_sl%field)
       call delete_spline_2d(plan_sl%spl_f)
       DEALLOCATE(plan_sl)
       plan_sl=>null()
    end if

  end subroutine delete_plan_adv_curvilinear

!======================================================================================================
!  creation of sll_SL_curvilinear type
!======================================================================================================

  ! grad_case  : integer, see function new_curvilinear_op
  ! carac_case : integer, see function new_plan_adv_curvilnear
  !              to choose the scheme for advection
  !                1 : using explicit Euler method with linear interpolation
  !                2 : using explicit Euler method with spline interpolation
  !                5 : using fixed point method
  
  function new_SL(geom_eta,delta_eta1,delta_eta2,dt,& 
 & N_eta1,N_eta2,grad_case,carac_case,bc,bc1_type,bc2_type) result(plan_sl)

    type(sll_SL_curvilinear), pointer :: plan_sl
    sll_real64, intent(in) :: geom_eta(2,2)
    sll_int32,  intent(in) :: N_eta1,N_eta2 
    sll_int32,  intent(in) :: bc(2)
    sll_int32,  intent(in) :: grad_case,carac_case,bc1_type,bc2_type
    sll_real64 :: delta_eta1,delta_eta2,dt
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max
    !sll_real64,dimension(:,:),pointer::jac_array

    sll_int32 :: err

    ALLOCATE(plan_sl)
    ALLOCATE(plan_sl%phi(N_eta1+1,N_eta2+1))

    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)

     plan_sl%poisson => new_plan_poisson_polar(delta_eta1,eta1_min,N_eta1,N_eta2,bc)
     plan_sl%grad => new_curvilinear_op(geom_eta,N_eta1,N_eta2,grad_case,bc1_type,bc2_type)
     plan_sl%adv => new_plan_adv_curvilinear(geom_eta,delta_eta1,delta_eta2,dt,N_eta1,N_eta2,carac_case,bc1_type,bc2_type)

  end function new_SL

!=================================================================================================
!  deletion of sll_SL_pola_eta1
!=================================================================================================

  ! subroutine delete_SL_curvilinear(plan_sl)
  ! deletion of sll_SL_curvilinear object
  subroutine delete_SL_curvilinear(plan_sl)

    implicit none

    type(sll_SL_curvilinear), pointer :: plan_sl

    sll_int32 :: err

    if (associated(plan_sl)) then
       call delete_plan_adv_curvilinear(plan_sl%adv)
       call delete_plan_poisson_polar(plan_sl%poisson)
       call delete_plan_curvilinear_op(plan_sl%grad)
       DEALLOCATE(plan_sl)
    end if

  end subroutine delete_SL_curvilinear

!==================================================================================
!  advection
!==================================================================================

  ! subroutine advect_CG_polar(plan,in,out,jac_array)
  ! compute step for Center-Guide equation
  ! plan :  sll_plan_adv_pola_eta1 object
  ! in   :  distribution function at time t_n, size (N_eta1+1)*(N_eta2+1)
  ! out  :  distribution function at time t_(n+1), size (N_eta1+1)*(N_eta2+1)
  ! plan%carac_case : integer to choose the scheme for advection
  !                   1 : using explicit Euler method with linear interpolation
  !                   2 : using explicit Euler method with spline interpolation
  !                   5 : using fixed point method
  ! jac_array: the jacobian of the tranformation point F, jac_array = det(DF) 


  subroutine advect_CG_curvilinear(plan,fn,fnp1,jac_array,mesh_case)

    implicit none

    type(sll_plan_adv_curvilinear) , pointer      :: plan
    sll_real64, dimension(:,:) ,     intent(in)   :: fn
    sll_real64, dimension(:,:) ,     intent(out)  :: fnp1
    sll_real64, dimension(:,:) ,     pointer      :: jac_array
    sll_int32,  intent(in)                        :: mesh_case
    
    sll_real64 :: eta1_loc,eta2_loc,eta1,eta1n,eta2,eta20,eta2n
    sll_real64 :: tolr,a_eta1,a_eta2,eta10,eta2_min,eta2_max !,tolth
    sll_real64 :: dt, delta_eta1, delta_eta2, eta1_min, eta1_max
    sll_real64 :: phi_loc(2,2)
    sll_int32  :: N_eta1, N_eta2
    sll_int32  :: i,j,maxiter,iter,k_eta1,k_eta2
    sll_int32  :: bc1_type,bc2_type,ii,jj
   
  
   
   
    N_eta1=plan%N_eta1
    N_eta2=plan%N_eta2
    dt=plan%dt
    delta_eta1=plan%delta_eta1
    delta_eta2=plan%delta_eta2
    eta1_min=plan%eta1_min
    eta1_max=plan%eta1_max
    eta2_min=plan%eta2_min
    eta2_max=plan%eta2_max
    bc1_type=plan%bc1_type
    bc2_type=plan%bc2_type
 
    !construction of spline coefficients for f
    call compute_spline_2D(fn,plan%spl_f)

    fnp1=0._f64
    
    if (plan%carac_case==1) then
       !explicit Euler with linear interpolation  
      
       fnp1 = fn 
       do j=1,N_eta2+1
          do i=1,N_eta1+1
             eta1=eta1_min+real(i-1,f64)*delta_eta1
             eta2=eta2_min+real(j-1,f64)*delta_eta2
      
             eta2=eta2-dt*plan%field(1,i,j)/jac_array(i,j)
             eta1=eta1-dt*plan%field(2,i,j)/jac_array(i,j)
             call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,&
                &eta1,eta2)
                eta1_loc=(eta1-eta1_min)/(eta1_max-eta1_min)
                eta1_loc=eta1_loc*real(N_eta1,f64)
                k_eta1=floor(eta1_loc)+1
                eta1_loc=eta1_loc-real(k_eta1-1,f64)
                if(((k_eta1-1).gt.(N_eta1)).or.((k_eta1-1).lt.0))then
                  print *,"#bad value of k_eta1=",k_eta1,N_eta1,eta1_loc,eta1,i,j,iter
                endif
                if((k_eta1-1)==N_eta1)then
                  k_eta1=N_eta1
                  if (abs(eta1_loc)>1.e-13) print *,'#eta1_loc=',eta1_loc
                  eta1_loc=1._f64
                endif
                
                eta2_loc=(eta2-eta2_min)/(eta2_max-eta2_min)
                eta2_loc=eta2_loc*real(N_eta2,f64)
                k_eta2=floor(eta2_loc)+1
                eta2_loc=eta2_loc-real(k_eta2-1,f64)
                if(((k_eta2-1).gt.(N_eta2)).or.((k_eta2-1).lt.0))then
                  print *,"#bad value of k_eta2=",k_eta2,N_eta2
                endif
                if((k_eta2-1)==N_eta2)then
                  k_eta2=N_eta2
                  if (abs(eta2_loc)>1.e-13) print *,'#eta2_loc=',eta2_loc
                  eta2_loc=1._f64
                endif

             
             fnp1(i,j)=(1.0_f64-eta2_loc)*((1.0_f64-eta1_loc)*fn(k_eta1,k_eta2) &
             &+eta1_loc*fn(k_eta1+1,k_eta2)) +eta2_loc*((1.0_f64-eta1_loc)* &
             & fn(k_eta1,k_eta2+1)+eta1_loc*fn(k_eta1+1,k_eta2+1))
            
          end do

       end do
  end if

  if (plan%carac_case==2) then
      !explicit Euler with "Spline interpolation"
       do i=1,N_eta1+1
          do j=1,N_eta2+1

             eta10=eta1_min+real(i-1,f64)*delta_eta1
             eta20=eta2_min+real(j-1,f64)*delta_eta2
        
             eta2=eta20-(dt*plan%field(1,i,j)/jac_array(i,j))
             eta1=eta10-(dt*plan%field(2,i,j)/jac_array(i,j))

             call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,&
                &eta1,eta2)
             fnp1(i,j)=interpolate_value_2d(eta1,eta2,plan%spl_f)

          end do
       end do
       
  end if

    
    if (plan%carac_case==5) then
       !using fixed point method
    
       !initialization
       maxiter=25
       tolr=1e-12
       eta1n=0.0_f64
       eta2n=0.0_f64

       do j=1,N_eta2+1 !N_eta2+1
          eta20=eta2_min+real(j-1,f64)*delta_eta2
          eta2=eta20
          k_eta2=j
          do i=1,N_eta1+1 !N_eta1+1
             eta10=eta1_min+real(i-1,f64)*delta_eta1
             eta1=eta10
             eta1_loc=0.0_f64
             eta2_loc=0.0_f64
             k_eta1=i
             a_eta1=0.0_f64
             a_eta2=0.0_f64
             iter=0
 
           do while (((iter<maxiter) .and. (abs((eta1n-eta1))+abs((eta2n-eta2))>tolr)).or.(iter==0))    
      
                eta1_loc=(eta1-eta1_min)/(eta1_max-eta1_min)
                eta1_loc=eta1_loc*real(N_eta1,f64)
                k_eta1=floor(eta1_loc)+1
                eta1_loc=eta1_loc-real(k_eta1-1,f64)
                if(((k_eta1-1).gt.(N_eta1)).or.((k_eta1-1).lt.0))then
                  print *,"#bad value of k_eta1=",k_eta1,N_eta1,eta1_loc,eta1,i,j,iter
                endif
                if((k_eta1-1)==N_eta1)then
                  k_eta1= N_eta1
                  if (abs(eta1_loc)>1.e-13) print *,'#eta1_loc=',eta1_loc
                  eta1_loc=1._f64 
                endif  

             
                eta2_loc=(eta2-eta2_min)/(eta2_max-eta2_min)
                eta2_loc=eta2_loc*real(N_eta2,f64)
                k_eta2=floor(eta2_loc)+1
                eta2_loc=eta2_loc-real(k_eta2-1,f64)
                if(((k_eta2-1).gt.(N_eta2)).or.((k_eta2-1).lt.0))then
                  print *,"#bad value of k_eta2=",k_eta2,N_eta2
                endif
                if((k_eta2-1)==N_eta2)then
                  k_eta2=N_eta2
                  if (abs(eta2_loc)>1.e-13) print *,'#eta2_loc=',eta2_loc
                  eta2_loc=1._f64 
                endif
             
               do jj=0,1
                 do ii=0,1               
     phi_loc(1+ii,1+jj) = plan%field(2,k_eta1+ii,k_eta2+jj)/jac_array(k_eta1+ii,k_eta2+jj)
                 enddo
               enddo
               phi_loc(1,1) = (1.0_f64-eta1_loc)*phi_loc(1,1) + eta1_loc*phi_loc(2,1)
               phi_loc(1,2) = (1.0_f64-eta1_loc)*phi_loc(1,2) + eta1_loc*phi_loc(2,2)
               a_eta1 = (1.0_f64-eta2_loc)*phi_loc(1,1) + eta2_loc*phi_loc(1,2)
               a_eta1 = 0.5_f64*dt*a_eta1

               do jj=0,1
                 do ii=0,1               
     phi_loc(1+ii,1+jj) = plan%field(1,k_eta1+ii,k_eta2+jj)/jac_array(k_eta1+ii,k_eta2+jj)
                 enddo
               enddo
               phi_loc(1,1) = (1.0_f64-eta1_loc)*phi_loc(1,1) + eta1_loc*phi_loc(2,1)
               phi_loc(1,2) = (1.0_f64-eta1_loc)*phi_loc(1,2) + eta1_loc*phi_loc(2,2)
               a_eta2 = (1.0_f64-eta2_loc)*phi_loc(1,1) + eta2_loc*phi_loc(1,2)
               a_eta2 = 0.5_f64*dt*a_eta2

                eta1n=eta1
                eta2n=eta2
                eta1=eta10-a_eta1
                eta2=eta20-a_eta2              
                call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,& 
                                 & eta2_min,eta2_max,eta1,eta2)                           
                iter=iter+1
                
             end do

             if (iter==maxiter .and. abs((eta1n-eta1))+abs((eta2n-eta2))>tolr) then
                print*,'#no convergence in fixed point method',i,j !,step
             end if
             eta1=eta10-2.0_f64*a_eta1
             eta2=eta20-2.0_f64*a_eta2
            call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max, & 
                             & eta2_min,eta2_max,eta1,eta2)   
           fnp1(i,j)=interpolate_value_2d(eta1,eta2,plan%spl_f)
          end do

       end do
  end if




    if(bc2_type==sll_PERIODIC) fnp1(:,N_eta2+1)=fnp1(:,1) 
    if(bc1_type==sll_PERIODIC) fnp1(N_eta1+1,:)=fnp1(1,:)
 

  end subroutine advect_CG_curvilinear

!****************************************************************************
!****************************************************************************
  ! subroutine SL_classic(plan,in,out)
  ! computes the classic semi-Lagrangian scheme for Vlasov-Poisson equation
  ! plan : sll_SL_pola_eta1 object, contains plan for Poisso, gradient and advection
  ! in   : distribution function at time n, size (N_eta1+1)*(N_eta2+1)
  ! out  : distribution function at time n+1, size (N_eta1+1)*(N_eta2+1)

  subroutine SL_order_1(plan,inn,outt,jac_array,step,mesh_case,N_eta1,N_eta2,geom_eta)

    implicit none

    type(sll_SL_curvilinear)   , intent(inout), pointer :: plan
    type(mudpack_2d)                                    :: poisson_df
    sll_real64,  dimension(:,:), intent(inout)          :: inn
    sll_real64,  dimension(:,:), intent(out)            :: outt
    sll_real64,  dimension(:,:), pointer                :: jac_array
    sll_real64,  intent(in)                             :: geom_eta(2,2)
    sll_int32,   intent(in)                             :: N_eta1, N_eta2
    sll_int32,   intent(in)                             :: step,mesh_case
    sll_int                                             :: i,j
    sll_real64                                          :: eta1_min,eta1_max,eta2_min,eta2_max

   
    eta1_min = geom_eta(1,1)
    eta1_max = geom_eta(2,1)
    eta2_min = geom_eta(1,2)
    eta2_max = geom_eta(2,2)
      
     !---*Fourier method for poisson's equation*-----
     !call poisson_solve_polar(plan%poisson,inn,plan%phi)
     !-----------------****------------------------------------
    call solve_poisson_colella_mudpack(poisson_df,plan%phi, inn)
   
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field,N_eta1,N_eta2,geom_eta)
    
    call advect_CG_curvilinear(plan%adv,inn,outt,jac_array,mesh_case)
   
  end subroutine SL_order_1

!!*********************************************************************************
  ! subroutine SL_ordre_2(plan,in,out)
  ! computes the semi-Lagrangian scheme order 2
  ! plan : sll_SL_pola_eta1 object, contains plan for Poisso, gradient and advection
  ! in   : distribution function at time n, size (N_eta1+1)*(N_eta2+1)
  ! out  : distribution function at time n+1, size (N_eta1+1)*(N_eta2+1)
  subroutine SL_order_2(plan,inn,outt,jac_array,mesh_case,N_eta1,N_eta2,geom_eta)

    implicit none

    type(sll_SL_curvilinear),   pointer     ::  plan
    sll_real64, dimension(:,:), intent(in)  :: inn
    sll_real64, dimension(:,:), intent(out) :: outt
    sll_real64, dimension(:,:), pointer     ::jac_array
    sll_real64, intent(in)                  :: geom_eta(2,2)
    sll_int32,  intent(in)                  :: N_eta1, N_eta2
    sll_int32,  intent(in)                  :: mesh_case
    sll_real64                              :: dt

    dt=plan%adv%dt
    
    plan%adv%dt=dt/2.0_f64

    !poisson_solve
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field,N_eta1,N_eta2,geom_eta)
    
    
    call advect_CG_curvilinear(plan%adv,inn,outt,jac_array,mesh_case)
  
    
    !!we just obtained f^(n+1/2)
    !poisson
    call compute_grad_field(plan%grad,plan%phi,plan%adv%field,N_eta1,N_eta2,geom_eta)
    !!we just obtained E^(n+1/2)
    plan%adv%dt=dt
    call advect_CG_curvilinear(plan%adv,inn,outt,jac_array,mesh_case)
    
  end subroutine SL_order_2



!=================================================================================================
!  Plot results !=================================================================================================
! visucase : 0 gnuplot
!            1 visit
  subroutine print2d(dom,ftab,Nx,Ny,visucase,step,filename)
    sll_int32,intent(in)::Nx,Ny,visucase,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx,0:Ny),intent(in)::ftab
    character(len=*),intent(in)::filename

    if(visucase==0)then
       !gnuplot
       call printgp2d(dom,ftab,Nx,Ny,step,filename)
    endif
    if(visucase==1)then
       !vtk
       call printvtk2d(dom,ftab,Nx,Ny,step,filename)
    endif
  end subroutine print2d

!******************************************************
  subroutine print2dper(dom,ftab,Nx,Ny,visucase,step,filename)
    sll_int32,intent(in)::Nx,Ny,visucase,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    character(len=*),intent(in)::filename

    if(visucase==0)then
       !gnuplot
       call printgp2dper(dom,ftab,Nx,Ny,step,filename)
    endif
    if(visucase==1)then
       !vtk
       call printvtk2dper(dom,ftab,Nx,Ny,step,filename)
    endif
  end subroutine print2dper

!********************************************************
  subroutine printgp2dper(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    sll_int32::i,j
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.dat'

    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    open(unit=900,file=str)
    do j=0,Ny-1
       do i=0,Nx-1
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          write(900,*) z(0),z(1),ftab(i,j)
       enddo
       i=Nx
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       write(900,*) z(0),z(1),ftab(0,j)      
       write(900,*) ''      
    enddo
    j=Ny
    do i=0,Nx-1
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       write(900,*) z(0),z(1),ftab(i,0)
    enddo
    i=Nx
    z(0)=dom(0,0)+real(i,f64)*dz(0)
    z(1)=dom(0,1)+real(j,f64)*dz(1)
    write(900,*)z(0),z(1),ftab(0,0)
    write(900,*)''
    close(900)  
  end subroutine printgp2dper

!************************************************************
  subroutine printgp2d(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx,0:Ny),intent(in)::ftab
    sll_int32::i,j
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.dat'

    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    open(unit=900,file=str)
    do j=0,Ny
       do i=0,Nx
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          write(900,*) z(0),z(1),z(1)*cos(z(0)),z(1)*sin(z(0)),ftab(i,j)
       enddo
       write(900,*) ''      
    enddo
    close(900)  
  end subroutine printgp2d

!!*****************************************************************
  subroutine printvtk2dper(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    sll_int32::i,j
    sll_int32,intent(in):: step
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
    write(str,*) 'f'//trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    !open(unit=900,file='f.vtk')
    open(unit=900,file=str,form='formatted')
    write(900,'(A)')                  '# vtk DataFile Version 2.0'
    write(900,'(A)')                  'Exemple'
    write(900,'(A)')                  'ASCII'
    write(900,'(A)')                  'DATASET STRUCTURED_POINTS'
    write(900,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx+1,' ', Ny+1,' ', 1
    write(900,'(A,I0,A,I0,A,I0)') 'ORIGIN ', floor(dom(0,0)+0.1),' ' , floor(dom(0,1)+0.1),' ' , 0
    !write(900,'(A,F10.4,A,F10.4,A,F10.4)') 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*) 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*)
    write(900,'(A,I0)')           'POINT_DATA ',(Nx+1)*(Ny+1)
    write(900,'(A,I0)')           'SCALARS f float ',1
    write(900,'(A)')                  'LOOKUP_TABLE default'

    do j=0,Ny-1
       do i=0,Nx-1
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          !write(900,'(F0.8)') ftab(i,j)
          write(900,*) ftab(i,j)
       enddo
       i=Nx
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       !write(900,'(F0.8)') ftab(0,j)            
       write(900,*) ftab(0,j)            
    enddo
    j=Ny
    do i=0,Nx-1
       z(0)=dom(0,0)+real(i,f64)*dz(0)
       z(1)=dom(0,1)+real(j,f64)*dz(1)
       !write(900,'(F0.8)') ftab(i,0)
       write(900,*) ftab(i,0)
    enddo
    i=Nx
    z(0)=dom(0,0)+real(i,f64)*dz(0)
    z(1)=dom(0,1)+real(j,f64)*dz(1)
    !write(900,'(F0.8)') ftab(0,0)	  	       
    write(900,*) ftab(0,0)
    close(900)  
  end subroutine printvtk2dper

!!*************************************************************
  subroutine printvtk2d(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx,0:Ny),intent(in)::ftab
    sll_int32::i,j
    sll_int32,intent(in):: step
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character(len=80)::str,str2
    write(str2,*)step
    !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
    write(str,*) 'f'//trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    !open(unit=900,file='f.vtk')
  
    open(unit=900,file=str,form='formatted')
    write(900,'(A)')                  '# vtk DataFile Version 2.0'
    write(900,'(A)')                  'Exemple'
    write(900,'(A)')                  'ASCII'
    write(900,'(A)')                  'DATASET STRUCTURED_POINTS'
    write(900,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx+1,' ', Ny+1,' ', 1
    write(900,'(A,I0,A,I0,A,I0)') 'ORIGIN ', floor(dom(0,0)+0.1),' ' , floor(dom(0,1)+0.1),' ' , 0
    !write(900,'(A,F10.4,A,F10.4,A,F10.4)') 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*) 'SPACING ', dz(0),' ', dz(1),' ', 0.25 
    write(900,*)
    write(900,'(A,I0)')           'POINT_DATA ',(Nx+1)*(Ny+1)
    write(900,'(A,I0)')           'SCALARS f float ',1
    write(900,'(A)')                  'LOOKUP_TABLE default'
    
    do j=0,Ny
       do i=0,Nx
          z(0)=dom(0,0)+real(i,f64)*dz(0)
          z(1)=dom(0,1)+real(j,f64)*dz(1)
          !write(900,'(F0.8)') ftab(i,j)
          write(900,*) ftab(i,j)
       enddo
    enddo
    close(900)  
  end subroutine printvtk2d





end module curvilinear_2D_advection
