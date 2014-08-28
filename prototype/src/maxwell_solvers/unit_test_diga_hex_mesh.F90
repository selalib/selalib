! Solve Maxwell System without sources using discontinuous
! Galerkin numerical method on hexagonal mesh.
! At t=0 we initialize Bz as a gaussian in the center of the domain
! and see what happens.
! We still need to find an analytic solution to check precision order.
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
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
sll_int32                   :: nstep = 1000
sll_real64                  :: cfl, dt 
sll_real64                  :: time
sll_real64, pointer         :: S_Ex(:,:)
sll_real64, pointer         :: X_Ex(:,:)
sll_real64, pointer         :: tmp_Ex(:,:)
sll_real64, pointer         :: tmp_Ey(:,:)
sll_real64, pointer         :: tmp_Bz(:,:)
sll_real64, pointer         :: tmp_Po(:,:)

num_cells = 10

print *, ""
print *, "Creating a mesh with 40 cells, mesh coordinates written in ./hex_mesh_coo.txt"
mesh => new_hex_mesh_2d(num_cells,             &
                        0.0_f64,               &
                        0.0_f64,               &
                        0.5_f64*sqrt(3.0_f64), &
                        0.5_f64,               &
                       -0.5_f64*sqrt(3.0_f64), &
                        0.5_f64,               &
                        0.0_f64,               &
                        1.0_f64,               &
                        10.0_f64 )

call sll_display(mesh)
call write_hex_mesh_2d(mesh,"hex_mesh_coo.txt")
call write_mtv(mesh)
print *, ""

degree = 2

call initialize(maxwell, mesh, degree)

cfl = 0.2
dt = cfl/sqrt(2./(mesh%delta/(degree+1))**2)

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

time = 0.0

SLL_CLEAR_ALLOCATE(tmp_Ex(maxwell%n_ddl,mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(tmp_Ey(maxwell%n_ddl,mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(tmp_Bz(maxwell%n_ddl,mesh%num_triangles), error)
SLL_CLEAR_ALLOCATE(tmp_Po(maxwell%n_ddl,mesh%num_triangles), error)

maxwell%Ex = 0d0
maxwell%Ey = 0d0
maxwell%Bz = exp(-(maxwell%x_ddl**2+maxwell%y_ddl**2))
maxwell%Jx = 0d0
maxwell%Jy = 0d0

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
!   call set_charge_and_currents(time)


   do k = 1, 5

      tmp_Ex = maxwell%D_Ex - maxwell%Jx
      tmp_Ey = maxwell%D_Ey - maxwell%Jy
      tmp_Bz = maxwell%D_Bz 
      !tmp_Po = maxwell%D_Po !+ maxwell%Ro

      ! Compute D_Ex, D_Ey, D_Bz, D_Po
      call solve(maxwell, mesh)
      !call set_charge_and_currents(time+C(k)*dt)

      maxwell%D_Ex = A(k)*tmp_Ex + dt * (maxwell%D_Ex - maxwell%Jx)
      maxwell%D_Ey = A(k)*tmp_Ey + dt * (maxwell%D_Ey - maxwell%Jy)
      maxwell%D_Bz = A(k)*tmp_Bz + dt * maxwell%D_Bz 
      !maxwell%D_Po = A(k)*tmp_Po + dt * maxwell%D_Po !+ maxwell%Ro)  ! xi=0

      maxwell%Ex = maxwell%Ex + B(k) * maxwell%D_Ex 
      maxwell%Ey = maxwell%Ey + B(k) * maxwell%D_Ey 
      maxwell%Bz = maxwell%Bz + B(k) * maxwell%D_Bz
      !maxwell%Po = maxwell%Po + B(k) * maxwell%D_Po

   end do

   time = time + dt
   if (mod(istep, 2) == 0) call plot_simple(maxwell, mesh)
   !call plot_double(maxwell, mesh)

   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g15.3,' s, ')",advance="no") time
   write(*,*)
   !write(*,"(' Ex error = ',g25.15)") maxval(abs(maxwell%Ex-sin(time)*maxwell%x_ddl*sin(sll_pi*maxwell%y_ddl)))

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

subroutine plot_simple( this, mesh )

   type(maxwell_dg_hex_mesh), intent(in) :: this
   type(hex_mesh_2d),         intent(in) :: mesh
   sll_int32, save                       :: iplot = 0
   sll_int32                             :: idl
   sll_int32                             :: iel
   character(len=4)                      :: cplot
   character(len=15)                     :: dat_file
   character(len=15)                     :: gnu_file

   iplot = iplot+1
   call int2string(iplot, cplot)

   dat_file = "fields_"//cplot//".dat"
   write(*,"(10x, 'Fichier de sortie GNUplot ', a)") dat_file

   gnu_file = "Ex_hex_mesh.gnu"
   open(83, file = gnu_file, position="append")
   if ( iplot == 1 ) rewind(83)
   write(83,"(a,g15.3,a)")"set title 'Time = ",time, "'"
   write(83,"(a)")"splot '"//dat_file//"' w lp "
   close(83)
   gnu_file = "Ey_hex_mesh.gnu"
   open(83, file = gnu_file, position="append")
   if ( iplot == 1 ) rewind(83)
   write(83,"(a,g15.3,a)")"set title 'Time = ",time, "'"
   write(83,"(a)")"splot '"//dat_file//"' u 1:2:4 w lp "
   close(83)
   gnu_file = "Bz_hex_mesh.gnu"
   open(83, file = gnu_file, position="append")
   if ( iplot == 1 ) rewind(83)
   write(83,"(a,g15.3,a)")"set title 'Time = ",time, "'"
   write(83,"(a)")"splot '"//dat_file//"' u 1:2:5 w lp "
   close(83)

   open(94, file = dat_file)

   do iel = 1, mesh%num_triangles

      do idl = 1, this%n_ddl
         write(94,*) sngl(this%x_ddl(idl,iel)), &
                     sngl(this%y_ddl(idl,iel)), &
                     sngl(this%Ex(idl,iel)),    &
                     sngl(this%Ey(idl,iel)),    &
                     sngl(this%Bz(idl,iel))
      end do

      idl = 1
      write(94,*) sngl(this%x_ddl(idl,iel)), &
                  sngl(this%y_ddl(idl,iel)), &
                  sngl(this%Ex(idl,iel)),    &
                  sngl(this%Ey(idl,iel)),    &
                  sngl(this%Bz(idl,iel))
      write(94,*)
      write(94,*)
   
   end do

   close(94)

end subroutine plot_simple

subroutine plot_double( this, mesh )

   type(maxwell_dg_hex_mesh), intent(in) :: this
   type(hex_mesh_2d),         intent(in) :: mesh
   sll_int32, save :: iplot = 0
   sll_int32       :: idl, iel
   character(len=4) :: cplot
   character(len=07) :: dat_file
   character(len=15) :: gnu_file

   iplot = iplot+1
   call int2string(iplot, cplot)

   dat_file = "Ex_"//cplot
   gnu_file = "Ex_hex_mesh.gnu"

   write(*,"(10x, 'Fichier de sortie GNUplot ', a)") dat_file//".dat"

   open(83, file = gnu_file, position="append")
   if ( iplot == 1 ) rewind(83)
   write(83,"(a,g15.3,a)")"set title 'Time = ",time, "'"
   write(83,"(a)")"splot '"//dat_file//"' w lp, '"//dat_file//"' u 1:2:4 w l"
   close(83)

   open(94, file = dat_file)

   do iel = 1, mesh%num_triangles

      do idl = 1, this%n_ddl
         write(94,*) sngl(this%x_ddl(idl,iel)), &
                     sngl(this%y_ddl(idl,iel)), &
                     sngl(this%Ex(idl,iel)),    &
                     sngl(-sin(time)*this%x_ddl(idl,iel)*sin(sll_pi*this%y_ddl(idl,iel)))
      end do

      idl = 1
      write(94,*) sngl(this%x_ddl(idl,iel)), &
                  sngl(this%y_ddl(idl,iel)), &
                  sngl(this%Ex(idl,iel)),    &
                  sngl(-sin(time)*this%x_ddl(idl,iel)*sin(sll_pi*this%y_ddl(idl,iel)))
      write(94,*)
      write(94,*)
   
   end do

   close(94)

end subroutine plot_double


subroutine write_mtv(mesh)

type(hex_mesh_2d), pointer :: mesh
real(8) :: coor(2,mesh%num_pts_tot)
integer :: ntri(3,mesh%num_triangles)
real(8) :: x1
real(8) :: y1
integer :: is1, is2, is3

open( 10, file="hex_mesh.mtv")

!--- Trace du maillage ---

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Maillage' "

do i = 1, mesh%num_triangles

   x1 = mesh%center_cartesian_coord(1, i)
   y1 = mesh%center_cartesian_coord(2, i)

   call get_cell_vertices_index( x1, y1, mesh, is1, is2, is3)

   ntri(1,i) = is1
   ntri(2,i) = is2
   ntri(3,i) = is3

   coor(1,is1) = mesh%global_to_x1(is1)
   coor(2,is1) = mesh%global_to_x2(is1)
   coor(1,is2) = mesh%global_to_x1(is2)
   coor(2,is2) = mesh%global_to_x2(is2)
   coor(1,is3) = mesh%global_to_x1(is3)
   coor(2,is3) = mesh%global_to_x2(is3)

   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)

end do

!--- Numeros des noeuds et des triangles

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des noeuds et des triangles' "

do i = 1, mesh%num_triangles
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_triangles
   x1 = (  coor(1,ntri(1,i))  &
         + coor(1,ntri(2,i))  &
     + coor(1,ntri(3,i))    )/3.
   y1 = (  coor(2,ntri(1,i))  &
         + coor(2,ntri(2,i))  &
     + coor(2,ntri(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(f8.5)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(f8.5)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

do i = 1, mesh%num_pts_tot
   x1 = coor(1,i)
   y1 = coor(2,i)
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=5 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

!--- Numeros des noeuds 

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des noeuds' "

do i = 1, mesh%num_triangles
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_pts_tot
   x1 = coor(1,i)
   y1 = coor(2,i)
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=5 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

!--- Numeros des triangles

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des triangles' "

do i = 1, mesh%num_triangles
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_triangles
   x1 = (  coor(1,ntri(1,i))  &
         + coor(1,ntri(2,i))  &
     + coor(1,ntri(3,i))    )/3.
   y1 = (  coor(2,ntri(1,i))  &
         + coor(2,ntri(2,i))  &
     + coor(2,ntri(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

write(10,*)"$END"
close(10)
   
end subroutine write_mtv

end program test_maxwell_dg_hex_mesh
