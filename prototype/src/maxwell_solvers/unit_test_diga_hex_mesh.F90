program test_maxwell_dg_hex_mesh

#include "sll_working_precision.h"
#include "sll_memory.h"
use hex_mesh
use sll_constants
use sll_maxwell_diga_hex_mesh

implicit none

type(maxwell_dg_hex_mesh)   :: maxwell
type(hex_mesh_2d), pointer  :: mesh
sll_int32                   :: num_cells
sll_int32                   :: error
sll_int32                   :: i, k
sll_int32                   :: degree
sll_real64                  :: A(5)
sll_real64                  :: B(5)
sll_real64                  :: C(5)
sll_int32                   :: istep
sll_int32                   :: nstep = 100
sll_real64                  :: dt = 0.01
sll_real64                  :: time
sll_real64, pointer         :: S_Ex(:,:)
sll_real64, pointer         :: X_Ex(:,:)
sll_real64, pointer         :: tmp_Ex(:,:)
sll_real64, pointer         :: tmp_Ey(:,:)
sll_real64, pointer         :: tmp_Bz(:,:)
sll_real64, pointer         :: tmp_Po(:,:)

num_cells = 20

print *, ""
print *, "Creating a mesh with 40 cells, mesh coordinates written in ./hex_mesh_coo.txt"
mesh => new_hex_mesh_2d(num_cells)
call sll_display(mesh)
call write_hex_mesh_2d(mesh,"hex_mesh_coo.txt")
print *, ""

degree = 1
call initialize(maxwell, mesh, degree)

!Low storage Runge Kutta order 4

A(1) = 0
A(2) = - 567301805773D0/1357537059087D0
A(3) = - 2404267990393D0/2016746695238D0
A(4) = - 3550918686646D0/2091501179385D0
A(5) = - 1275806237668D0/842570457699D0

B(1) = 1432997174477D0/9575080441755D0
B(2) = 5161836677717D0/13612068292357D0
B(3) = 1720146321549D0/2090206949498D0
B(4) = 3134564353537D0/4481467310338D0
B(5) = 2277821191437D0/14882151754819D0

C(1) = 0D0
C(2) = 1432997174477D0/9575080441755D0
C(3) = 2526269341429D0/6820363962896D0
C(4) = 2006345519317D0/3224310063776D0
C(5) = 2802321613138D0/2924317926251D0

do i = 1, 5
   print*, C(i)
end do
time = 0.0

SLL_CLEAR_ALLOCATE(tmp_Ex(maxwell%n_ddl,mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(tmp_Ey(maxwell%n_ddl,mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(tmp_Bz(maxwell%n_ddl,mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(tmp_Po(maxwell%n_ddl,mesh%num_triangles), error)

do istep = 1, nstep


   !maxwell%Ex = sin(time)*maxwell%x_ddl*sin(sll_pi*maxwell%y_ddl)
   !maxwell%Ey = sin(time)*maxwell%y_ddl*sin(sll_pi*maxwell%x_ddl)
   !maxwell%Bz = (cos(time)-1)*(sll_pi*maxwell%y*cos(sll_pi*maxwell%x_ddl) &
   !            -sll_pi*maxwell%x_ddl*cos(sll_pi*maxwell%y_ddl))

      
!   call rksetup()
!   maxwell%D_Ex = maxwell%Ex
!   call accumulate(1._f64/6._f64)
!   call rkstage(0.5_f64)
!   maxwell%D_Ex = maxwell%Ex
!   call accumulate(1._f64/3._f64)
!   call rkstage(0.5_f64)
!   maxwell%D_Ex = maxwell%Ex
!   call accumulate(1._f64/3._f64)
!   call rkstage(1.0_f64)
!   maxwell%D_Ex = maxwell%Ex
!   call accumulate(1._f64/6._f64)
!   call rkstep()
   call set_charge_and_currents(time)

   do k = 1, 5

      tmp = maxwell%D_Ex - maxwell%Jx
      call solve(maxwell, mesh)
      call set_charge_and_currents(time+C(k)*dt)

      maxwell%D_Ex = maxwell%Ex
      maxwell%D_Ex = A(k)*tmp + dt * maxwell%D_Ex

      !maxwell%D_Ex = A(k)*maxwell%D_Ex + dt * (maxwell%D_Ex - maxwell%Jx)
      !maxwell%D_Ey = A(k)*maxwell%D_Ey + dt * (maxwell%D_Ey - maxwell%Jy)
      !maxwell%D_Bz = A(k)*maxwell%D_Bz + dt * (maxwell%D_Bz)
      !maxwell%D_Po = A(k)*maxwell%D_Po + dt * (maxwell%D_Po + maxwell%Ro)  ! xi=1

      maxwell%Ex = maxwell%Ex + B(k) * maxwell%D_Ex 

      !maxwell%Ey = maxwell%Ey + B(k) * maxwell%D_Ey 
      !maxwell%Bz = maxwell%Bz + B(k) * maxwell%D_Bz
      !maxwell%Po = maxwell%Po + B(k) * maxwell%D_Po

   end do

   time = time + dt

   !error = maxval(abs(maxwell%Ex - sin(time)*maxwell%x_ddl*sin(sll_pi*maxwell%y_ddl)))

   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g15.3,' s, ')",advance="no") time
   write(*,"(' Ex error = ',g25.15)") maxval(maxwell%Ex)- exp(1d0)

end do

contains

subroutine set_charge_and_currents(t)

   sll_real64, intent(in) :: t

   maxwell%Jx = ((cos(t)-1)*(sll_pi*cos(sll_pi*maxwell%x_ddl) &
               +sll_pi*sll_pi*maxwell%x_ddl*sin(sll_pi*maxwell%y_ddl)) &
               -cos(t)*maxwell%x_ddl*sin(sll_pi*maxwell%y_ddl))
   maxwell%Jy = ((cos(t)-1)*(sll_pi*cos(sll_pi*maxwell%y_ddl) &
               +sll_pi*sll_pi*maxwell%y_ddl*sin(sll_pi*maxwell%x_ddl)) &
               -cos(t)*maxwell%y_ddl*sin(sll_pi*maxwell%x_ddl))
   maxwell%Ro = sin(t)*(sin(sll_pi*maxwell%y_ddl)+sin(sll_pi*maxwell%x_ddl))

end subroutine set_charge_and_currents



subroutine rksetup()

   S_Ex = 0.0_f64

   X_Ex = maxwell%Ex 

end subroutine rksetup

subroutine rkstage(coef)

   sll_real64, intent(in) :: coef

   maxwell%Ex = X_Ex + coef * dt * maxwell%D_Ex

end subroutine rkstage

subroutine accumulate(coef)

   sll_real64, intent(in) :: coef

   S_Ex =  S_Ex + coef * maxwell%D_Ex

end subroutine accumulate

subroutine rkstep()

   maxwell%Ex = X_Ex + dt * S_Ex

end subroutine rkstep

end program test_maxwell_dg_hex_mesh
