!> @brief 
!> Selalib 2D poisson solver for cartesian coordinates, Dirichlet BC's
!   
!> @authors                    
!> Nhung PHAM 
!> Edwin CHACON-GOLCHER
!                                  
!> @details
!> Add in here all sequential poisson solvers that use Dirichlet boundary
!> conditions.
!**************************************************************************

module sll_poisson_2d_dirichlet_cartesian

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
!  use sll_fft
  use sll_constants
   implicit none

  !> Structure to store data from Poisson solver. This
  !> solver is parallel on structured cartesian mesh. 
  type poisson_2d_dirichlet_cartesian
     sll_int32                           :: ncx    !< number of cells  
     sll_int32                           :: ncy    !< number of cells

     sll_real64                          :: Lx     !< domain length 
     sll_real64                          :: Ly     !< domain length
     sll_real64, dimension(:), pointer   :: vkgd
     sll_real64, dimension(:), pointer   :: vkgi
     sll_real64, dimension(:), pointer   :: vkgs
     sll_int32,dimension(:),pointer      :: prof
     sll_int32,dimension(:),pointer      :: indi
     sll_int32,dimension(:),pointer      :: indj
     sll_int32,dimension(:),pointer      :: kld
     sll_real64, dimension(:),pointer    :: val
     sll_real64, dimension(:),pointer    :: f
     sll_real64, dimension(:),pointer    :: vx
     sll_int32                           :: nsky
  end type poisson_2d_dirichlet_cartesian

  interface sll_delete
     module procedure delete_poisson_2d_dirichlet_cartesian
  end interface sll_delete

contains

  function new_poisson_2d_dirichlet_cartesian_plan( &
    ncx, &            
    ncy, &            
    Lx, &    
    Ly ) result(plan)

    type (poisson_2d_dirichlet_cartesian), pointer :: plan
    sll_int32                        :: ncx          !< number of cells in x
    sll_int32                        :: ncy          !< number of cells in y
    sll_real64                       :: Lx           !< length x
    sll_real64                       :: Ly           !< length y
    !sll_real64, dimension(:,:)          :: vkgs
    sll_int32                        :: nn,n,ch
    sll_int32 :: ii,jj,i,j,p,k,l
    sll_int32     :: ierr
    sll_real64    ::v
    SLL_ALLOCATE(plan, ierr)


    plan%ncx = ncx
    plan%ncy = ncy
    plan%Lx  = Lx
    plan%Ly  = Ly
    nn=ncx
    n=nn*nn
    ch=nn*nn+4*nn*(nn-1)
    SLL_ALLOCATE( plan%indi(ch), ierr )
    SLL_ALLOCATE( plan%indj(ch), ierr )
    SLL_ALLOCATE( plan%val(ch), ierr )
    p=0
    do k=1,n
       p=p+1
       plan%indi(p)=k
       plan%indj(p)=k
       plan%val(p)=4
    end do
    do ii=1,nn-1
       do jj=1,nn
          p=p+1
          i=ii+(jj-1)*nn
          j=ii+1+(jj-1)*nn
          plan%val(p)=-1
          plan%indi(p)=i
          plan%indj(p)=j
          p=p+1
          plan%val(p)=-1
          plan%indi(p)=j
          plan%indj(p)=i
       end do
    end do
    do ii=1,nn
       do jj=1,nn-1
          i=ii+(jj-1)*nn
          j=ii+jj*nn
          p=p+1
          plan%val(p)=-1
          plan%indi(p)=i
          plan%indj(p)=j
          p=p+1
          plan%val(p)=-1
          plan%indi(p)=j
          plan%indj(p)=i
       end do
    end do
  SLL_ALLOCATE(plan%prof(n),ierr)
    SLL_ALLOCATE(plan%kld(n+1),ierr)
    plan%prof=0
    do k=2,ch
       i=plan%indi(k)
       j=plan%indj(k)
       plan%prof(j)=max(plan%prof(j),j-i)
       plan%prof(i)=max(plan%prof(i),i-j)
    end do
    plan%kld(1)=1
    do i=2,n+1
       plan%kld(i)=plan%kld(i-1)+plan%prof(i-1)
    end do
    plan%nsky=plan%kld(n+1)-1
    SLL_ALLOCATE(plan%vkgs(plan%nsky),ierr)
    SLL_ALLOCATE(plan%vkgi(plan%nsky),ierr)
    SLL_ALLOCATE(plan%vkgd(n),ierr)
    plan%vkgs=0.d0
    plan%vkgi=0.d0
    !definir vkgs,vkgi,vkgd
    do l=1,ch
       i=plan%indi(l)
       j=plan%indj(l)
       v=plan%val(l)
       if(i.eq.j)then
          plan%vkgd(i)=v
       end if
       if (i<j) then
          k=plan%kld(j+1)-j+i
          plan%vkgs(k)=v
       end if
       if (i>j) then
          k=plan%kld(i+1)-i+j
          plan%vkgi(k)=v
       end if
    end do
    SLL_ALLOCATE(plan%f(n),ierr)
    SLL_ALLOCATE(plan%vx(n),ierr)

  end function new_poisson_2d_dirichlet_cartesian_plan




  !> Note that the equation that is solved is: \f$ \Delta \phi = \rho \f$
  !> Thus the user is responsible for giving the proper sign to the source term.
  subroutine solve_poisson_2d_dirichlet_cartesian(plan, rho, phi)
    type (poisson_2d_dirichlet_cartesian), pointer :: plan !< self object
    sll_real64, dimension(:,:)        :: rho      !< charge density
    sll_real64, dimension(:,:)        :: phi      !< electric potential
    sll_int32                         :: ncx      !< global size
    sll_int32                         :: i
    sll_int32                         :: ii, jj
    sll_real64                        :: void
    sll_int32                         :: nn,n,ier
    sll_int32                         :: ifac,isol,mp,nsym
    !sll_real64, dimension(:,:)          :: vkgs
    !sll_int32, dimen                  :: indi,indj
    ifac=1  ! we compute the LU decomposition
    isol=1  ! we solve the linear system
    nsym=1  ! we do not take into account the symetry of M
    mp=6   ! write the log on screen
    ncx=size(rho,1)
    nn=ncx
    n=nn*nn
    !rho(1,1)=rho(1,1)+1
    do ii=1,nn
       do jj=1,nn
          i=ii+(jj-1)*nn
          plan%f(i)=rho(ii,jj)
       end do
    end do
    call sol(plan%vkgs,plan%vkgd,plan%vkgi,plan%f,plan%kld,plan%vx,&
         n,mp,ifac,isol,nsym,void,ier,plan%nsky)
    do ii=1,nn
       do jj=1,nn
          i=ii+(jj-1)*nn
          phi(ii,jj)=plan%vx(i)
       end do
    end do
  end subroutine solve_poisson_2d_dirichlet_cartesian


!> Delete the Poisson solver object
  subroutine delete_poisson_2d_dirichlet_cartesian(plan)
    type (poisson_2d_dirichlet_cartesian), pointer :: plan
    sll_int32                                           :: ierr

    if( .not. associated(plan) ) then
       print *, 'ERROR, delete_poisson_2d_dirichlet_cartesian_plan(): ', &
            'passed plan is not associated.'
       STOP
    end if

    SLL_DEALLOCATE(plan, ierr)
  end subroutine delete_poisson_2d_dirichlet_cartesian



end module sll_poisson_2d_dirichlet_cartesian
