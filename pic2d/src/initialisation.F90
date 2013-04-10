module initialisation
#include "sll_working_precision.h"
use zone

implicit none

sll_int32, private :: i, j

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialisation des champs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init( tm )

type (tm_mesh_fields) :: tm
sll_real64 :: aux1, aux2

select case (nomcas)

   case ("viry__")
      call statio( tm, 0.d0 )
      tm%ex = 0.d0; tm%ey = 0.d0; tm%bz = 0.d0
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0
!-- onde perodique en x et t 
   case ("statio")
      call statio( tm, 0.d0 )
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

!-- onde entrant dans un guide d'onde
   case ("entran")
      call entran( tm, 0.d0 )
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

   case ("Eincom")
      call Eincom( tm, 0.d0 )
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

   case ("Sincom")
      call Sincom( tm, 0.d0 )
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

   case ("Nincom")
      call Nincom( tm, 0.d0 )
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

   case("faisce")
      tm%ex = 0.d0; tm%ey = 0.d0; tm%bz = 0.d0
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

!-- gaussienne stationnaire
   case("gaussv")
      tm%ex = 0.d0; tm%ey = 0.d0; tm%bz = 0.d0
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0

!-- gaussienne perturbee
   case("plasma")
      tm%ex = 0.d0; tm%ey = 0.d0; tm%bz = 0.d0
      tm%jx = 0.d0; tm%jy = 0.d0; tm%r0 = 0.d0; tm%r1 = 0.d0
      do i=0,nx-1
         aux1 = alpha/kx * sin(kx*x(i))
         aux2 = alpha * cos(kx*x(i))
         do j=0,ny
            tm%ex(i,j) = aux1
            tm%r1(i,j) = aux2
         enddo
      enddo
      
end select

end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sol_exacte( sol, time )

type (tm_mesh_fields) :: sol
sll_real64 :: time

select case (nomcas)

!-- onde perodique en x et t 
   case ("statio")
      call statio( sol, time )

!-- onde entrant dans un guide d'onde
   case ("entran")
      call entran( sol, time )

   case ("Eincom")
      call Eincom( sol, time )

   case ("Sincom")
      call Sincom( sol, time )

   case ("Nincom")
      call Nincom( sol, time )

end select

end subroutine sol_exacte

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine statio( sol, time )

type (tm_mesh_fields) :: sol
sll_real64, parameter :: mm=4, nn=2
sll_real64 :: time
sll_real64 :: omega, alpha, beta

! cas periodique en temps et en espace
!
! Ex = - c**2*n*pi/omega/dimy*cos(m*pi*x/dimx)*sin(n*pi*y/dimy)*sin(omega t)
! Ey =   c**2*m*pi/omega/dimx*sin(m*pi*x/dimx)*cos(n*pi*y/dimy)*sin(omega t)
! Bz =                    cos(m*pi*x/dimx)*cos(n*pi*y/dimy)*cos(omega t)

alpha = mm*pi/dimx
beta  = nn*pi/dimy
omega = sqrt( alpha*alpha + beta*beta )

do i=0,nx-1
do j=0,ny
    sol%ex(i,j) = - beta/omega * cos(alpha * 0.5*(x(i)+x(i+1))) &
         & * sin(beta * y(j)) * sin(omega * time)
end do
end do

do i=0,nx
do j=0,ny-1
    sol%ey(i,j) =  alpha/omega * sin(alpha * x(i)) &
         & * cos(beta * 0.5*(y(j)+y(j+1))) * sin(omega * time)
end do  
end do  

do i=0,nx-1
do j=0,ny-1
    sol%bz(i,j) = cos(alpha * 0.5*(x(i)+x(i+1))) &
         & * cos(beta * 0.5*(y(j)+y(j+1))) * cos(omega * time)
end do  
end do  

end subroutine statio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine entran( sol, time )

type (tm_mesh_fields) ::  sol
sll_real64 :: time
sll_int32, parameter :: mm=1 !si modif, faire aussi dans CL (maxwell.f90)
sll_real64 :: alpha, omega

!onde entrant dans le guide d'onde    
!Ex = - c**2*n*pi/(w*dimy)*sin(n*pi*y/dimy)*sin(wt-kx)  
!Ey =   c**2*k/w*cos(n*pi*y/dimy)*cos(wt-kx)
!Bz =   cos(n*pi*y/dimy)*cos(wt-kx)

alpha = 2*mm*pi/dimx
omega = alpha
 
do i=0,nx-1
   do j=0,ny-1
      sol%ex(i,j) = 0.d0
      
      sol%ey(i,j) = cos(omega*time - alpha*x(i))
   
      sol%bz(i,j) = cos(omega*time - alpha*0.5*(x(i)+x(i+1)))
   end do
   sol%ex(i,ny) = 0.d0
end do

do j=0,ny-1
   sol%ey(nx,j) = cos(omega*time - alpha*dimx)
end do


end subroutine entran


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Eincom( sol, time )

type (tm_mesh_fields) ::  sol
sll_real64 :: time
sll_int32, parameter :: mm=1 !si modif, faire aussi dans CL (maxwell.f90)
sll_real64 :: alpha

!onde entrant dans le guide d'onde par l'est  

alpha = 2*mm*pi/dimx
 
do i=0,nx-1
   do j=0,ny-1
      sol%ex(i,j) = 0.d0
      sol%ey(i,j) = cos(alpha*(time + x(i)))
      sol%bz(i,j) = -cos(alpha*(time + 0.5*(x(i)+x(i+1))))
   end do
   sol%ex(i,ny) = 0.d0
end do

do j=0,ny-1
   sol%ey(nx,j) = cos(alpha*(time + dimx))
end do

end subroutine Eincom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Sincom( sol, time )

type (tm_mesh_fields) ::  sol
sll_real64 :: time
sll_int32, parameter :: mm=1 !si modif, faire aussi dans CL (maxwell.f90)
sll_real64 :: alpha

!onde entrant dans le guide d'onde par le sud 

alpha = 2*mm*pi/dimy
 
do i=0,nx-1
   do j=0,ny-1
      sol%ex(i,j) = cos(alpha*(time - y(j)))
      sol%ey(i,j) = 0.d0
      sol%bz(i,j) = -cos(alpha*(time - 0.5*(y(j)+y(j+1))))
   end do
   sol%ex(i,ny) = cos(alpha*(time - dimy))
end do

do j=0,ny-1
   sol%ey(nx,j) = 0.d0
end do

end subroutine Sincom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Nincom( sol, time )

type (tm_mesh_fields) ::  sol
sll_real64 :: time
sll_int32, parameter :: mm=1 !si modif, faire aussi dans CL (maxwell.f90)
sll_real64 :: alpha

!onde entrant dans le guide d'onde par le nord

alpha = 2*mm*pi/dimy
 
do i=0,nx-1
   do j=0,ny-1
      sol%ex(i,j) = cos(alpha*(time + y(j)))
      sol%ey(i,j) = 0.d0
      sol%bz(i,j) = cos(alpha*(time + 0.5*(y(j)+y(j+1))))
   end do
   sol%ex(i,ny) = cos(alpha*(time + dimy))
end do

do j=0,ny-1
   sol%ey(nx,j) = 0.d0
end do

end subroutine Nincom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module initialisation
