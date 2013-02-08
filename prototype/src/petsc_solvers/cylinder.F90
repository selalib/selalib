!
!   Laplacian in 3D. Modeled by the partial differential equation
!
!   in a cylindrical conducting pipe
!
!   with boundary conditions
!
!   u = 1 for r = 1, periodic in theta and z
!
!   This uses multigrid to solve the linear system

program cylinder

#include <finclude/petscdef.h>
use petsc
implicit none

KSP                :: ksp
Mat                :: A
PetscErrorCode     :: ierr
DM                 :: da
Vec                :: x
PetscInt,parameter :: i1=1,i3=16, i5=5
PetscLogDouble     :: t1, t2

external    ComputeRHS, ComputeMatrix

call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
call PetscGetCPUTime(t1, ierr)

!Creates an object that will manage the communication of three-dimensional 
!regular array data that is distributed across some processors.
call DMDACreate3d(PETSC_COMM_WORLD,  &
                  DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_PERIODIC,DMDA_BOUNDARY_PERIODIC, &
                  DMDA_STENCIL_STAR,       &
                  i3,3*i3,i3,     & !global dimensions in each direction
                  PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, & !corresponding number of proc
                  i1,i1, &
                  PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)  

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

end program cylinder

subroutine ComputeRHS(da,x,b,ierr)

implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>


PetscErrorCode :: ierr
PetscInt       :: mr,ma,mz
Vec            :: x,b
DM             :: da
PetscScalar    :: r, Hr, Ha, Hz
PetscScalar, pointer    :: rho(:)
PetscInt       :: is, js, ks, im, jm, km
PetscInt       :: i, j, k, l

call DMDAGetInfo(da,PETSC_NULL_INTEGER,mr,ma,mz,                       &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

call DMDAGetCorners(da,is,js,ks,im,jm,km,ierr)

Hr = 1d0 / mr
Ha = 1d0 / (ma-1)
Hz = 1d0 / (mz-1)

call VecGetArrayF90(b,rho,ierr)

l = 0
do k=ks,ks+km-1
   do j=js,js+jm-1
      do i=is,is+im-1
         l = l+1
         r = (i+1)*Hr*(i+1)*Hr
         if (r <= 0.2) then
            rho(l) = - 4d0 * PETSC_PI * r
         else
            rho(l) = 0.
         end if
      end do
   end do
end do

call VecRestoreArrayF90(b,rho,ierr)

return
end


subroutine ComputeMatrix(da, x, jac,str, ierr)

implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscdmda.h>

Mat            :: jac
Vec            :: x
PetscErrorCode :: ierr
DM             :: da
PetscInt       :: i,j,k,mr,ma,mz
PetscInt       :: xm,ym,zm,xs,ys,zs,i1,i2,i5,i6,i7
PetscScalar    :: v(7),Hr,Ha,Hz,HrHr,HaHa,HzHz
PetscScalar    :: rc, rl, rr
MatStencil     :: row(4),col(4,7)
MatStructure   :: str

i1 = 1; i2 = 2; i5 = 5; i6 = 6; i7 = 7

call DMDAGetInfo(da,PETSC_NULL_INTEGER,mr,ma,mz,                       &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
                 PETSC_NULL_INTEGER,ierr)

Hr = 1d0 / mr
Ha = 1d0 / (ma-1)
Hz = 1d0 / (mz-1)
HrHr = Hr*Hr; HaHa = Ha*Ha; HzHz = Hz*Hz

call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr)        
  
do k=zs,zs+zm-1    
   do j=ys,ys+ym-1
      do i=xs,xs+xm-1

         row(MatStencil_i) = i  ; row(MatStencil_j) = j; row(MatStencil_k) = k 

         col(MatStencil_i,1)=i  ; col(MatStencil_j,1)=j  ; col(MatStencil_k,1)=k-1
         col(MatStencil_i,2)=i  ; col(MatStencil_j,2)=j  ; col(MatStencil_k,2)=k+1
         col(MatStencil_i,3)=i  ; col(MatStencil_j,3)=j-1; col(MatStencil_k,3)=k
         col(MatStencil_i,4)=i  ; col(MatStencil_j,4)=j+1; col(MatStencil_k,4)=k
         col(MatStencil_i,5)=i-1; col(MatStencil_j,5)=j  ; col(MatStencil_k,5)=k
         col(MatStencil_i,6)=i+1; col(MatStencil_j,6)=j  ; col(MatStencil_k,6)=k
         col(MatStencil_i,7)=i  ; col(MatStencil_j,7)=j  ; col(MatStencil_k,7)=k

         rl = i * i * Hr * Hr
         rc = (i+1) * (i+1) * Hr * Hr
         rr = (i+2) * (i+2) * Hr * Hr

         v(1) = rc / HzHz
         v(2) = rc / HzHz
         v(3) = 1d0 / (rc*HaHa)
         v(4) = 1d0 / (rc*HaHa)

         if (i == 0) then

            v(5) = - (2.d0/(rc*HaHa)+2d0*rc/HzHz + (rr+rc)/(rr-rc)/rc)
            v(6) =   (rr+rc)/(rr-rc)/rc
            col(MatStencil_i,5)=i
            col(MatStencil_i,6)=i+1
            call MatSetValuesStencil(jac,i1,row,i6,col,v,INSERT_VALUES,ierr)

         else if (i == mr-1) then

            v(5) = 1.d0
            col(MatStencil_i,5)=i
            call MatSetValuesStencil(jac,i1,row,i5,col,v,INSERT_VALUES,ierr)

         else

            v(5) = (rc+rl)/(rc-rl)/(rr-rl)
            v(6) = (rr+rc)/(rr-rc)/(rr-rl)
            v(7) = -2.d0/(rc*HaHa)-2d0*rc/HzHz -v(5)-v(6)
   
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
#include <finclude/petscdmda.h>
#include <finclude/petscvec.h90>
#include <finclude/petscviewer.h>

PetscErrorCode   :: ierr
DM               :: da
PetscScalar      :: r, a, z
PetscInt         :: mr, ma, mz
PetscViewer      :: viewer
PetscInt         :: i, j, k, l
PetscInt         :: is, js, ks, im, jm, km
Vec              :: x, c
PetscScalar      :: rmin = .2, rmax = .4
PetscScalar      :: amin = 0., amax = 2.*PETSC_PI
PetscScalar      :: zmin = 0., zmax = 1.

character(len=72)      :: chaine
PetscScalar, pointer   :: coor(:)

call DMDAGetInfo(da,PETSC_NULL_INTEGER,mr,ma,mz,          &
               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

call DMDASetUniformCoordinates(da, rmin, rmax, amin,    &
                                   amax, zmin, zmax, ierr)

call DMDAGetCorners(da,is,js,ks,im,jm,km,ierr)

call DMDAGetCoordinates(da, c, ierr)
call VecGetArrayF90(c,coor,ierr)
!  Create coordinates vector data
l = 0
do i= is,is+im-1
   do j=js,js+jm-1
      do k=ks,ks+km-1
         r = coor(l+1)
         a = coor(l+2)
         z = coor(l+3)
         l = l+1; coor(l) = r * cos(a)
         l = l+1; coor(l) = r * sin(a)
         l = l+1; coor(l) = z
      end do
   end do
end do

call VecRestoreArrayF90(c,coor,ierr)

call PetscViewerCreate(PETSC_COMM_WORLD, viewer, ierr)
call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK, ierr)
call PetscViewerFileSetName(viewer, "cylinder.vtk", ierr)
call DMDASetCoordinates(da, c, ierr)
call DMView(da,viewer,ierr)
call VecView(x,viewer,ierr)
call PetscViewerDestroy(viewer, ierr)

write(chaine,"('mr,ma,mz =',3i4,1X,i6,'\n')")mr,ma,mz, mr*ma*mz
call PetscPrintf(PETSC_COMM_WORLD, trim(chaine), ierr)

call PetscViewerDestroy(viewer, ierr)

end subroutine PlotSolution

