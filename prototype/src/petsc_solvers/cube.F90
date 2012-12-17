!
!   Laplacian in 3D. Modeled by the partial differential equation
!
!   with boundary conditions
!
!   u = 1 for x,y,z = 1
!
program cube

#include <finclude/petscdef.h>
use petsc
implicit none

KSP                :: ksp
DM                 :: da
Vec                :: x, b
PetscInt           :: i1=1,i3=-32
PetscLogDouble     :: t1, t2
PetscErrorCode     :: ierr

external ComputeRHS,ComputeMatrix

call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
call PetscGetCPUTime(t1, ierr)

!Creates an object that will manage the communication of three-dimensional 
!regular array data that is distributed across some processors.
call DMDACreate3d(PETSC_COMM_WORLD, &
                  DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, &
                  DMDA_STENCIL_STAR,i3,i3,i3, &
                  PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, &
                  i1,i1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)  

call DMSetFunction(da,ComputeRHS,ierr)
call DMSetJacobian(da,ComputeMatrix,ierr)   


call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
call KSPSetDM(ksp,da,ierr)
call KSPSetFromOptions(ksp,ierr)
call KSPSolve(ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
call KSPGetSolution(ksp,x,ierr)
call PlotSolution( da, x, ierr)

call KSPDestroy(ksp,ierr)
call DMDestroy(da,ierr)
call PetscGetCPUTime(t2, ierr)
call PetscFinalize(ierr)
print"('Code took', g15.3,' CPU seconds')", t2-t1

end program cube


subroutine ComputeRHS(da, x, b,ierr)

implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>

PetscErrorCode :: ierr
PetscInt       :: mz,mx,my
Vec            :: x, b
DM             :: da
PetscScalar    :: HxHyHz

call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,mz,                        &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
                 PETSC_NULL_INTEGER,ierr)

HxHyHz = 1d0 / (mx-1d0) / (my-1d0) / (mz-1d0)
call VecSet(b, HxHyHz, ierr)

return
end


subroutine ComputeMatrix(da, x, jac, str, ierr)

implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscdmda.h>


Mat            :: jac
PetscErrorCode :: ierr
DM             :: da
PetscInt       :: i,j,k,l
PetscInt       :: xm,ym,zm,xs,ys,zs,i1,i7
PetscScalar    :: v(7)
MatStencil     :: row(4),col(4,7)
PetscInt       :: mx,my,mz
PetscScalar    :: Hx,Hy,Hz,HxHydHz,HzHxdHy,HyHzdHx
Vec            :: x
MatStructure   :: str

i1 = 1
i7 = 7

call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,mz,        &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                 PETSC_NULL_INTEGER,ierr)

Hx = 1d0 / (mx-1)
Hy = 1d0 / (my-1)
Hz = 1d0 / (mz-1)

HxHydHz = Hx*Hy/Hz
HzHxdHy = Hz*Hx/Hy
HyHzdHx = Hy*Hz/Hx

call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr)        
  
do k=zs,zs+zm-1    
   do j=ys,ys+ym-1
      do i=xs,xs+xm-1

         row(MatStencil_i) = i  ; row(MatStencil_j) = j; row(MatStencil_k) = k 

         if ( i==0 .or. i==mx-1 .or. j==0 .or. j==my-1 .or. k==0 .or. k==mz-1 ) then
            !v(1) = Hx*Hy*Hz
            v(1) = 2.d0*(HxHydHz + HzHxdHy + HyHzdHx)
            call MatSetValuesStencil(jac,i1,row,i1,row,v,INSERT_VALUES,ierr)
         else
            col(MatStencil_i,1)=i  ; col(MatStencil_j,1)=j  ; col(MatStencil_k,1)=k-1
            col(MatStencil_i,2)=i  ; col(MatStencil_j,2)=j  ; col(MatStencil_k,2)=k+1
            col(MatStencil_i,3)=i  ; col(MatStencil_j,3)=j-1; col(MatStencil_k,3)=k
            col(MatStencil_i,4)=i  ; col(MatStencil_j,4)=j+1; col(MatStencil_k,4)=k
            col(MatStencil_i,5)=i-1; col(MatStencil_j,5)=j  ; col(MatStencil_k,5)=k
            col(MatStencil_i,6)=i+1; col(MatStencil_j,6)=j  ; col(MatStencil_k,6)=k
            col(MatStencil_i,7)=i  ; col(MatStencil_j,7)=j  ; col(MatStencil_k,7)=k

            v(1) = HxHydHz 
            v(2) = HxHydHz
            v(3) = HzHxdHy
            v(4) = HzHxdHy
            v(5) = HyHzdHx
            v(6) = HyHzdHx 
            v(7) = - 2d0*(HxHydHz+HzHxdHy+HyHzdHx)
            call MatSetValuesStencil(jac,i1,row,i7,col,v,INSERT_VALUES,ierr)
         end if

      end do
   end do
end do
call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
return 
end


subroutine PlotSolution( da, x, ierr)
implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscviewer.h>
#include <finclude/petscdmda.h>

PetscErrorCode   :: ierr
DM               :: da
PetscInt         :: mx, my, mz
PetscViewer      :: viewer
PetscInt         :: i, j, k, l, csize
PetscInt         :: is, js, ks, im, jm, km
Vec              :: x
PetscScalar      :: xmin = 0., xmax = 1.
PetscScalar      :: ymin = 0., ymax = 1.
PetscScalar      :: zmin = 0., zmax = 1.


call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,mz,                       &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,ierr)


call DMDASetUniformCoordinates(da, xmin, xmax, ymin,    &
                                 ymax, zmin, zmax, ierr)

call PetscViewerCreate(PETSC_COMM_WORLD, viewer, ierr)
call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK, ierr)
call PetscViewerFileSetName(viewer, "cube.vtk", ierr)
call DMView(da,viewer,ierr)
call VecView(x,viewer,ierr)
call PetscViewerDestroy(viewer, ierr)

end subroutine PlotSolution
