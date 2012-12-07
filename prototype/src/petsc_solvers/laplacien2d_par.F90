module f90module
#include <finclude/petscdmdef.h>
use petscdmdef
type userctx
   type(DM)     :: da
   PetscInt     :: xs,xe,xm,gxs,gxe,gxm
   PetscInt     :: ys,ye,ym,gys,gye,gym
   PetscInt     :: mx,my
   PetscMPIInt  :: rank
end type userctx

end module f90module

program main
#include <finclude/petscdef.h>
use petsc
use f90module
implicit none

MPI_Comm                :: comm
type(Vec)               :: X,F,localF
type(Mat)               :: J
type(KSP)               :: ksp
PetscInt                :: Nx,Ny,N,mx,my,ifive,ithree
PetscBool               :: flg

double precision        :: rtol
PetscInt                :: max_nonlin_its,one
PetscInt                :: lin_its
PetscInt                :: m
PetscScalar             :: mone
PetscErrorCode          :: ierr
type(userctx)           :: user
type(userctx), pointer  :: puser

mone           = -1.d0
rtol           = 1.d-8
max_nonlin_its = 10
one            = 1
ifive          = 5
ithree         = 3

call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
comm = PETSC_COMM_WORLD
call MPI_Comm_rank(PETSC_COMM_WORLD,user%rank,ierr)


mx = 4
my = 4
call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-mx',user%mx,flg,ierr)
call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-my',user%my,flg,ierr)
N = mx*my

call KSPCreate(comm,ksp,ierr)

Nx = PETSC_DECIDE
Ny = PETSC_DECIDE
call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-Nx',Nx,flg,ierr)
call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-Ny',Ny,flg,ierr)
call DMDACreate2d(comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,   &
&     DMDA_STENCIL_STAR,user%mx,user%my,Nx,Ny,one,one,          &
&     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,user%da,ierr)

call DMDAGetInfo(user%da,PETSC_NULL_INTEGER,user%mx,user%my,    &
&               PETSC_NULL_INTEGER,                             &
&               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
&               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
&               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
&               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
&               PETSC_NULL_INTEGER,ierr)

!
!   Visualize the distribution of the array across the processors
!
call DMView(user%da,PETSC_VIEWER_DRAW_WORLD,ierr)


!  Get local grid boundaries (for 2-dimensional DMDA)
call DMDAGetCorners(user%da,user%xs,user%ys,PETSC_NULL_INTEGER, &
&     user%xm,user%ym,PETSC_NULL_INTEGER,ierr)
call DMDAGetGhostCorners(user%da,user%gxs,user%gys,             &
&     PETSC_NULL_INTEGER,user%gxm,user%gym,                     &
&     PETSC_NULL_INTEGER,ierr)

!  Here we shift the starting indices up by one so that we can easily
!  use the Fortran convention of 1-based indices (rather 0-based indices).
user%xs  = user%xs+1
user%ys  = user%ys+1
user%gxs = user%gxs+1
user%gys = user%gys+1

user%ye  = user%ys+user%ym-1
user%xe  = user%xs+user%xm-1
user%gye = user%gys+user%gym-1
user%gxe = user%gxs+user%gxm-1

call DMCreateGlobalVector(user%da,X,ierr)
call VecDuplicate(X,F,ierr)
call DMCreateLocalVector(user%da,localF,ierr)

call VecGetLocalSize(X,m,ierr)
!call DMCreateMatrix(user%da,MATAIJ,J,ierr)

call MatCreateAIJ(comm,m,m,N,N,ifive,PETSC_NULL_INTEGER,ithree, &
&                 PETSC_NULL_INTEGER,J,ierr)
  
!Set runtime options (e.g., -ksp_monitor -ksp_rtol <rtol> -ksp_type <type>)
  
call KSPSetFromOptions(ksp,ierr)

call ComputeRHS(F,user,ierr)
call ComputeMatrix(J,user,ierr)

call KSPSetOperators(ksp,J,J,SAME_NONZERO_PATTERN,ierr)
call KSPSolve(ksp,F,X,ierr)

call KSPGetIterationNumber(ksp,lin_its,ierr)
print*,'linear solve iterations = ',lin_its

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     VTK print out
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
call PetscViewerCreate(comm, viewer, ierr)
call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK, ierr)
call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK, ierr)
call PetscViewerFileSetName(viewer, "laplacien.vtk", ierr)
call DMDASetUniformCoordinates(user%da, xmin, xmax, ymin, ymax, &
                             PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)
call DMView(user%da,viewer,ierr)
call VecView(X,viewer,ierr)
call PetscViewerDestroy(viewer, ierr)

call MatDestroy(J,ierr)
call VecDestroy(X,ierr)
call VecDestroy(localF,ierr)
call VecDestroy(F,ierr)
call KSPDestroy(ksp,ierr)
call DMDestroy(user%da,ierr)
call PetscFinalize(ierr)
end

subroutine ComputeMatrix(lap,user,ierr)
#include <finclude/petscdef.h>
use petsc
use f90module
implicit none

Mat            :: lap
PetscErrorCode :: ierr
PetscInt       :: i,j
PetscInt       :: row,ione
PetscInt       :: col(5),ifive
PetscScalar    :: two,one
PetscScalar    :: v(5),hx,hy,hxdhy
PetscScalar    :: hydhx
type (userctx) :: user

ione   = 1
ifive  = 5
one    = 1.d0
two    = 2.d0
hx     = one/(mx-1)
hy     = one/(my-1)
hxdhy  = hx/hy        
hydhx  = hy/hx

do j=user%ys,user%ye
   row = (j - user%gys)*user%gxm + user%xs - user%gxs - 1
   do i=user%xs,user%xe
      row = row + 1
      if (i == 1 .or. j == 1 .or. i == user%mx .or. j == user%my) then
         col(1) = row
         v(1)   = one
         call MatSetValuesLocal(lap,ione,row,ione,col,v,INSERT_VALUES,ierr)
      else
         v(1) = -hxdhy
         v(2) = -hydhx
         v(3) = two*(hydhx + hxdhy)
         v(4) = -hydhx
         v(5) = -hxdhy
         col(1) = row - user%gxm
         col(2) = row - 1
         col(3) = row
         col(4) = row + 1
         col(5) = row + user%gxm
         call MatSetValuesLocal(lap,ione,row,ifive,col,v,INSERT_VALUES,ierr)
      endif
   end do
end do

call MatAssemblyBegin(lap,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(lap,MAT_FINAL_ASSEMBLY,ierr)
return 
end

subroutine ComputeRHS(F,user,ierr)
#include <finclude/petscdef.h>
use petsc

type(Vec)           :: X,F
PetscErrorCode      :: ierr
type (userctx)      :: user
PetscScalar,pointer :: lf_v(:)

call VecGetArrayF90(F,lf_v,ierr)

call ComputeRHSLocal(lf_v,user,ierr)

call VecRestoreArrayF90(F,lf_v,ierr)

end subroutine ComputeRHS

subroutine ComputeRHSLocal(x,f,user,ierr)
#include <finclude/petscsysdef.h>
use petscsys
use f90module
type (userctx) :: user
PetscScalar    :: f(user%xs:user%xe,user%ys:user%ye)
PetscErrorCode :: ierr
PetscScalar    :: two,one,hx,hy,hxdhy,hydhx
PetscInt       :: i,j

one    = 1.0
two    = 2.0
hx     = one/dble(user%mx-1)
hy     = one/dble(user%my-1)
hxdhy  = hx/hy
hydhx  = hy/hx

do j=user%ys,user%ye
   do i=user%xs,user%xe
      if (i == 1 .or. j == 1 .or. i == user%mx .or. j == user%my) then
         f(i,j) = 1.0
      else
         f(i,j) = 1.0
      endif
   end do
end do
ierr = 0
return
end
