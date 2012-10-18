module poisson_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use sll_tridiagonal
  use numeric_constants
  implicit none
  !>type sll_plan_poisson_polar
  !>type for the Poisson solver in polar coordinate
  type sll_plan_poisson_polar
     sll_real64 :: dr, rmin
     sll_int32 :: nr, ntheta
     sll_int32 :: bc(2)
     type(sll_fft_plan), pointer :: pfwd,pinv
     sll_real64, dimension(:,:), allocatable :: f_fft
     sll_comp64, dimension(:), allocatable :: fk,phik
     !for the tridiagonal solver
     sll_real64, dimension(:), allocatable :: a,cts
     sll_int32, dimension(:), allocatable :: ipiv
  end type sll_plan_poisson_polar

  !flags for boundary conditions
  !>boundary conditions can be at TOP_ or BOT_ and take value NEUMANN or DIRICHLET
  !>ex : TOP_DIRICHLET or BOT_NEUMANN
  integer, parameter :: DIRICHLET=1
  integer, parameter :: NEUMANN=2
  integer, parameter :: NEUMANN_MODE0=3

contains

!========================================
!  creation of sll_plan_poisson_polar
!========================================

  !>new_plan_poisson_polar(dr,rmin,nr,ntheta)
  !>build a sll_plan_poisson_polar object for the Poisson solver in polar coordinate
  !>dr : size of space in direction r
  !>rmin : interior radius
  !>nr and ntheta : number of space in direction r and theta
  !>bc : boundary conditions, can be combined with +
  !>bc is optionnal and default is Dirichlet condition in rmin and rmax
  function new_plan_poisson_polar(dr,rmin,nr,ntheta,bc) result(this)

    implicit none

    sll_real64 :: dr, rmin
    sll_int32 :: nr, ntheta
    sll_int32, optional :: bc(2)
    type(sll_plan_poisson_polar), pointer :: this

    sll_int32 :: err
    sll_real64, dimension(:), allocatable :: buf

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(buf(ntheta),err)
    SLL_ALLOCATE(this%f_fft(nr+1,ntheta+1),err)
    SLL_ALLOCATE(this%fk(nr+1),err)
    SLL_ALLOCATE(this%phik(nr+1),err)
!!$    SLL_ALLOCATE(this%a(3*(nr+1)),err)
!!$    SLL_ALLOCATE(this%cts(7*(nr+1)),err)
!!$    SLL_ALLOCATE(this%ipiv(nr+1),err)
    SLL_ALLOCATE(this%a(3*(nr-1)),err)
    SLL_ALLOCATE(this%cts(7*(nr-1)),err)
    SLL_ALLOCATE(this%ipiv(nr-1),err)

    this%dr=dr
    this%rmin=rmin
    this%nr=nr
    this%ntheta=ntheta
    if (present(bc)) then
       this%bc=bc
    else
       this%bc(1)=-1
       this%bc(2)=-1
    end if

    this%pfwd => fft_new_plan(ntheta,buf,buf,FFT_FORWARD,FFT_NORMALIZE)
    this%pinv => fft_new_plan(ntheta,buf,buf,FFT_INVERSE)
    SLL_DEALLOCATE_ARRAY(buf,err)

  end function new_plan_poisson_polar

!======================================
! deletion of sll_plan_poisson_polar
!======================================

  !>delete_plan_poisson_polar(plan)
  !>delete a sll_plan_poisson_polar object
  subroutine delete_plan_poisson_polar(this)

    implicit none

    type(sll_plan_poisson_polar), intent(inout), pointer :: this
    sll_int32 :: err
    if (associated(this)) then
       call fft_delete_plan(this%pfwd)
       call fft_delete_plan(this%pinv)
       SLL_DEALLOCATE_ARRAY(this%f_fft,err)
       SLL_DEALLOCATE_ARRAY(this%fk,err)
       SLL_DEALLOCATE_ARRAY(this%phik,err)
       SLL_DEALLOCATE_ARRAY(this%a,err)
       SLL_DEALLOCATE_ARRAY(this%cts,err)
       SLL_DEALLOCATE_ARRAY(this%ipiv,err)
       SLL_DEALLOCATE(this,err)
    end if

  end subroutine delete_plan_poisson_polar

!===================
!  Poisson solver
!===================

  !>subroutine poisson_solve_polar(plan,f,phi)
  !>poisson solver for polar system : -\Delta (phi)=f
  !>plan : sll_plan_poisson_polar, contains data for the solver
  !>f : distribution function, size (nr+1)*(ntheta+1), input
  !>phi : unknown field, size (nr+1)*(ntheta+1), output
  !>initialization must be done outside the solver
  subroutine poisson_solve_polar(plan,f,phi)

    implicit none

    type(sll_plan_poisson_polar), intent(inout), pointer :: plan
    sll_real64, dimension(plan%nr+1,plan%ntheta+1), intent(in) :: f
    sll_real64, dimension(plan%nr+1,plan%ntheta+1), intent(out) :: phi

    sll_real64 :: rmin,dr
    sll_int32 :: nr, ntheta,bc(2)

    sll_real64 :: r
    sll_int32::i,k,ind_k
    sll_real64:: kval

    nr=plan%nr
    ntheta=plan%ntheta
    rmin=plan%rmin
    dr=plan%dr

    bc = plan%bc
    plan%f_fft=f


    do i=1,nr+1
       call fft_apply_plan(plan%pfwd,plan%f_fft(i,1:ntheta),plan%f_fft(i,1:ntheta))
    end do

   ! poisson solver
    do k=0,ntheta-1!ntheta/2
       
       ind_k=k
       !do i=1,nr+1
       if( ind_k .gt. ntheta/2 ) then
                ind_k = ind_k - ntheta
       end if
       kval=real(ind_k,f64)
       !kval=1.5

       do i=2,nr
          r=rmin+real(i-1,f64)*dr
          plan%a(3*(i-1))=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
          plan%a(3*(i-1)-1)=2.0_f64/dr**2+(kval/r)**2
          plan%a(3*(i-1)-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)


!          plan%a(3*(i-1))=-plan%a(3*(i-1))
!          plan%a(3*(i-1)-1)=-plan%a(3*(i-1)-1)
!          plan%a(3*(i-1)-2)=-plan%a(3*(i-1)-2)

!!$          plan%a(3*i)=-1.0_f64/dr**2-1.0_f64/(2.0_f64*dr*r)
!!$          plan%a(3*i-1)=2.0_f64/dr**2+(real(ind_k,f64)/r)**2
!!$          plan%a(3*i-2)=-1.0_f64/dr**2+1.0_f64/(2.0_f64*dr*r)
          plan%fk(i)=fft_get_mode(plan%pfwd,plan%f_fft(i,1:ntheta),k)!ind_k)          
       enddo
       
       !if(k==0)then
       !  print *,'second membre en input'
       !  print *,plan%fk
       !endif
       
       !print *,k,sum(abs(plan%fk(1:nr+1)))
       plan%phik=0.0_f64
       !plan%a(1)=0.0_f64
       !plan%a(3*(nr-1))=0.0_f64

        
        !boundary condition at rmin
        if(bc(1)==DIRICHLET)then !Dirichlet
          plan%a(1)=0.0_f64
        endif
        if(bc(1)==NEUMANN)then
          plan%a(2)=plan%a(2)+plan%a(1) !Neumann
          plan%a(1)=0._f64
        endif
        if(bc(1)==NEUMANN_MODE0)then 
          if(k==0)then!Neumann for mode zero
            plan%a(2)=plan%a(2)+plan%a(1)
            plan%a(1)=0._f64
          else !Dirichlet for other modes
            plan%a(1)=0._f64
          endif  
        endif

        !boundary condition at rmax
        if(bc(2)==DIRICHLET)then !Dirichlet
          plan%a(3*(nr-1))=0.0_f64
        endif
        if(bc(2)==NEUMANN)then
          plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1)) !Neumann
          plan%a(3*(nr-1))=0.0_f64
        endif
        if(bc(2)==NEUMANN_MODE0)then 
          if(k==0)then!Neumann for mode zero
            plan%a(3*(nr-1)-1)=plan%a(3*(nr-1)-1)+plan%a(3*(nr-1))
            plan%a(3*(nr-1))=0.0_f64
          else !Dirichlet for other modes
            plan%a(3*(nr-1))=0.0_f64
          endif  
        endif
         
          
       call setup_cyclic_tridiag(plan%a,nr-1,plan%cts,plan%ipiv)
       call solve_cyclic_tridiag(plan%cts,plan%ipiv,plan%fk(2:nr),nr-1,plan%phik(2:nr))



        !boundary condition at rmin
        if(bc(1)==1)then !Dirichlet
          plan%phik(1)=0.0_f64
        endif
        if(bc(1)==2)then
          plan%phik(1)=plan%phik(2) !Neumann
        endif
        if(bc(1)==3)then 
          if(k==0)then!Neumann for mode zero
            plan%phik(1)=plan%phik(2)
          else !Dirichlet for other modes
            plan%phik(1)=0.0_f64
          endif  
        endif

        !boundary condition at rmax
        if(bc(2)==1)then !Dirichlet
          plan%phik(nr+1)=0.0_f64
        endif
        if(bc(2)==2)then
          plan%phik(nr+1)=plan%phik(nr) !Neumann
        endif
        if(bc(2)==3)then 
          if(k==0)then!Neumann for mode zero
            plan%phik(nr+1)=plan%phik(nr)
          else !Dirichlet for other modes
            plan%phik(nr+1)=0.0_f64
          endif  
        endif






       do i=1,nr+1
          call fft_set_mode(plan%pinv,phi(i,1:ntheta),plan%phik(i),k)!ind_k)
       end do
       !print *,k,'s',bc,sum(abs(plan%phik(1:nr+1)))

       !if(k==0)then
         !print *,'output'
         !print *,plan%phik
       !endif


    end do

    ! FFT INVERSE
    do i=1,nr+1
       call fft_apply_plan(plan%pinv,phi(i,1:ntheta),phi(i,1:ntheta))
    end do

    phi(:,ntheta+1)=phi(:,1)

  end subroutine poisson_solve_polar

end module poisson_polar
