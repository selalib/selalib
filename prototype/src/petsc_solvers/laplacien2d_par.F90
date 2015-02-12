#define PETSC_USE_FORTRAN_MODULES 1
program main
#include <finclude/petscdef.h>
use petsc
use petsc_module
implicit none

MPI_Comm            :: comm
Vec                 :: X,F,U
Mat                 :: J
KSP                 :: ksp
PetscInt            :: Nx,Ny,mx,my,ifive,ithree
PetscBool           :: flg

PetscInt            :: one
PetscInt            :: its
PetscScalar         :: mone
PetscErrorCode      :: ierr
type(userctx)       :: user
PetscViewer         :: viewer
PetscScalar         :: xmin, xmax, ymin, ymax
double precision    :: norm
PetscLogDouble      :: t1, t2
PetscInt            :: psize, prank
character(len=72)   :: chaine

mone           = -1.d0
one            = 1
ifive          = 5
ithree         = 3

call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
call PetscGetCPUTime(t1, ierr)
comm = PETSC_COMM_WORLD
call MPI_Comm_rank(PETSC_COMM_WORLD,prank,ierr)
call MPI_Comm_size(PETSC_COMM_WORLD,psize,ierr)

user%rank = prank

write(chaine,"('Hello World from node',2i3)") prank, psize
call PetscSynchronizedPrintf(PETSC_COMM_WORLD, trim(chaine)//'\n', ierr)
call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr)


mx = 32
my = 32
call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-mx',user%mx,flg,ierr)
call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-my',user%my,flg,ierr)

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
call VecDuplicate(X,U,ierr)
call DMCreateMatrix(user%da,MATAIJ,J,ierr)

!Multigrid
!call DMSetFunction(da,ComputeRHS,ierr)
!call DMSetJacobian(da,ComputeMatrix,ierr)
!call KSPSetDM(ksp,da,ierr)
!call KSPSetFromOptions(ksp,ierr)
!call KSPSolve(ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
!call KSPGetSolution(ksp,x,ierr)
  
call KSPSetFromOptions(ksp,ierr)
call ComputeRHS(F,user,ierr)
call ComputeSolution(U,user,ierr)
call ComputeMatrix(J,user,ierr)
call KSPSetOperators(ksp,J,J,SAME_NONZERO_PATTERN,ierr)
call KSPSolve(ksp,F,X,ierr)
call KSPGetIterationNumber(ksp,its,ierr)

!     VTK print out
xmin = 0.; xmax = 1.
ymin = 0.; ymax = 1.
call PetscViewerCreate(comm, viewer, ierr)
call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK, ierr)
call PetscViewerFileSetName(viewer, "laplacien.vtk", ierr)
call DMDASetUniformCoordinates(user%da, xmin, xmax, ymin, ymax, &
                             PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)
call DMView(user%da,viewer,ierr)
call VecView(X,viewer,ierr)
call VecView(U,viewer,ierr)
call VecView(F,viewer,ierr)
call PetscViewerDestroy(viewer, ierr)

!  Check the error
call VecAXPY(X,mone,U,ierr)
call VecNorm(X,NORM_2,norm,ierr)
if (user%rank .eq. 0) then
   if (norm .gt. 1.e-12) then
       write(6,100) norm,its
   else
      write(6,110) its
   endif
endif
100 format('Norm of error ',e11.4,' iterations ',i5)
110 format('Norm of error < 1.e-12,iterations ',i5)

call PetscGetCPUTime(t2, ierr)
write(chaine,"('KSP code took', f15.3,' CPU seconds \n')") (t2-t1)*psize
call PetscPrintf(PETSC_COMM_WORLD, chaine, ierr)

call MatDestroy(J,ierr)
call VecDestroy(X,ierr)
call VecDestroy(U,ierr)
call VecDestroy(F,ierr)
call KSPDestroy(ksp,ierr)
call DMDestroy(user%da,ierr)
call PetscFinalize(ierr)
end

subroutine ComputeMatrix(lap,user,ierr)
#include <finclude/petscdef.h>
use petsc
use petsc_module
implicit none

Mat            :: lap
PetscErrorCode :: ierr
PetscInt       :: i,j
PetscInt       :: row
PetscInt       :: ione
PetscInt       :: col(5)
PetscInt       :: ifive
PetscScalar    :: two,one
PetscScalar    :: v(5)
PetscScalar    :: hydhx,hxdhy
PetscScalar    :: hx,hy
type (userctx) :: user

ione   = 1
ifive  = 5
one    = 1.d0
two    = 2.d0
hx     = one/(dble(user%mx-1))
hy     = one/(dble(user%my-1))
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
use petsc_module
implicit none

Vec                 :: F
PetscErrorCode      :: ierr
type (userctx)      :: user
PetscScalar,pointer :: lf_v(:)

call VecGetArrayF90(F,lf_v,ierr)

call ComputeRHSLocal(lf_v,user,ierr)

call VecRestoreArrayF90(F,lf_v,ierr)

end subroutine ComputeRHS

subroutine ComputeRHSLocal(f,user,ierr)
#include <finclude/petscsysdef.h>
use petscsys
use petsc_module
implicit none
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
      f(i,j) = 8*PETSC_PI*PETSC_PI*sin(2*PETSC_PI*(i-1)*hx)*sin(2*PETSC_PI*(j-1)*hy)*hx*hy
   end do
end do
ierr = 0
return
end

subroutine ComputeSolution(U,user,ierr)
#include <finclude/petscdef.h>
use petsc
use petsc_module
implicit none
Vec                 :: U
PetscErrorCode      :: ierr
type (userctx)      :: user
PetscScalar,pointer :: pu(:)

call VecGetArrayF90(U,pu,ierr)
call ComputeSolutionLocal(pu,user,ierr)
call VecRestoreArrayF90(U,pu,ierr)

ierr = 0
return
end

subroutine ComputeSolutionLocal(u,user,ierr)
#include <finclude/petscsysdef.h>
use petscsys
use petsc_module
implicit none
type (userctx) :: user
PetscScalar    :: u(user%xs:user%xe,user%ys:user%ye)
PetscErrorCode :: ierr
PetscScalar    :: one,hx,hy
PetscInt       :: i,j

one    = 1.0
hx     = one/dble(user%mx-1)
hy     = one/dble(user%my-1)

do j=user%ys,user%ye
   do i=user%xs,user%xe
      u(i,j) = sin(2*PETSC_PI*(i-1)*hx)*sin(2*PETSC_PI*(j-1)*hy)
   end do
end do
ierr = 0

return
end
