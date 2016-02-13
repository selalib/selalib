!  Copyright INRIA
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
!**************************************************************
!> @copyright INRIA
!> @ingroup maxwell_solvers
!> @brief
!> Implements the Maxwell solver in 2D with FDTD method
!>
!> @details
!> Equation
!>
!>\f$ \displaystyle \frac{\partial E_x}{\partial t} = {c^2} \frac{\partial B_z}{\partial y} \f$,
!>
!>\f$\displaystyle \frac{\partial E_y}{\partial t} = -{c^2} \frac{\partial B_z}{\partial x} \f$
!>
!>\f$\displaystyle \frac{\partial B_z}{\partial t} =  \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x}  \f$.
!>
!>FDTD scheme
!>
!>\f$\displaystyle B_{z}^{n+1/2} = B_z^{n-1/2} - \Delta t \big( \frac{\partial E_y^n}{\partial x}
!>- \frac{\partial E_x^n}{\partial y} \big)\f$
!>
!>\f$\displaystyle E_x^{n+1} = E_x^{n} + c^2\Delta t  \frac{\partial B_z^{n+1/2}}{\partial y} \f$
!>
!>\f$\displaystyle E_y^{n+1} = E_y^{n} - c^2\Delta t  \frac{\partial B_z^{n+1/2}}{\partial x} \f$
!>
module sll_m_maxwell_2d_fdtd

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  implicit none

  public :: &
    sll_o_create, &
    sll_t_maxwell_2d_fdtd, &
    sll_o_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Initialize maxwell solver 2d with FDTD scheme
interface sll_o_create
 module procedure initialize_maxwell_2d_fdtd
 module procedure initialize_maxwell_2d_fdtd_alt
end interface sll_o_create
!> Solve maxwell solver 2d with FDTD scheme
interface sll_o_solve
 module procedure solve_maxwell_2d_fdtd
end interface sll_o_solve
!> Solve Ampere-Maxwell equation
interface sll_solve_ampere
 module procedure ampere_2d_fdtd
end interface sll_solve_ampere
!> Solve Faraday equation
interface sll_solve_faraday
 module procedure ampere_2d_fdtd
end interface sll_solve_faraday


!> @brief Object with data to solve Maxwell equation 
!> Maxwell in TE mode: (Ex,Ey,Bz)
type :: sll_t_maxwell_2d_fdtd
  private
  sll_int32  :: nc_eta1      !< x cells number
  sll_int32  :: nc_eta2      !< y cells number
  sll_int32  :: polarization !< TE or TM
  sll_real64 :: e_0          !< electric conductivity
  sll_real64 :: mu_0         !< magnetic permeability
  sll_real64 :: c            !< speed of light
  sll_real64 :: eta1_min     !< left side 
  sll_real64 :: eta1_max     !< right side
  sll_real64 :: delta_eta1   !< step size
  sll_real64 :: eta2_min     !< bottom side
  sll_real64 :: eta2_max     !< top side
  sll_real64 :: delta_eta2   !< step size
  sll_int32  :: i1           !< first indice of the block dimension 1
  sll_int32  :: j1           !< last indice of the block dimension 1
  sll_int32  :: i2           !< first indice of the block dimesnion 2
  sll_int32  :: j2           !< last indice of the block dimension 2
  sll_real64 :: dx           !< step size along dimension 1
  sll_real64 :: dy           !< step size along dimension 2
end type sll_t_maxwell_2d_fdtd


contains

!>Initilialize the maxwell solver
subroutine initialize_maxwell_2d_fdtd_alt(self, x1, x2, nc_x, &
                                      y1, y2, nc_y, polarization )

   type(sll_t_maxwell_2d_fdtd) :: self         !< maxwell solver object
   sll_real64            :: x1           !< first incidice along x
   sll_real64            :: y1           !< last indice along x
   sll_real64            :: x2           !< first indice along y
   sll_real64            :: y2           !< last indice along y
   sll_int32             :: nc_x         !< size step along y
   sll_int32             :: nc_y         !< size step along y
   sll_int32             :: polarization !< TE or TM

   self%c   = 1.0_f64
   self%e_0 = 1.0_f64
   
   self%dx  = (x2-x1)/nc_x
   self%dy  = (y2-y1)/nc_y
   self%i1  = 1
   self%j1  = nc_x+1
   self%i2  = 1
   self%j2  = nc_y+1

   self%polarization = polarization

end subroutine initialize_maxwell_2d_fdtd_alt

!>Initilialize the maxwell solver
subroutine initialize_maxwell_2d_fdtd(self, i1, j1, i2, j2, dx, dy, polarization )

   type(sll_t_maxwell_2d_fdtd) :: self        !< maxwell solver object
   sll_int32          :: i1             !< first incidice along x
   sll_int32          :: j1             !< last indice along x
   sll_int32          :: i2             !< first indice along y
   sll_int32          :: j2             !< last indice along y
   sll_real64         :: dx             !< size step along x
   sll_real64         :: dy             !< size step along y
   sll_int32          :: polarization   !< TE or TM
   !sll_int32        :: error           !< error code

   self%c   = 1.0_f64
   self%e_0 = 1.0_f64
   
   self%dx = dx
   self%dy = dy
   self%i1 = i1
   self%j1 = j1
   self%i2 = i2
   self%j2 = j2
   self%polarization = polarization

end subroutine initialize_maxwell_2d_fdtd

!> self routine exists only for testing purpose. Use ampere and faraday
!> in your appication.
subroutine solve_maxwell_2d_fdtd(self, fx, fy, fz, dt)

   type(sll_t_maxwell_2d_fdtd)         :: self !< maxwell object
   sll_real64, dimension(:,:) :: fx   !< Ex or Bx
   sll_real64, dimension(:,:) :: fy   !< Ey or By
   sll_real64, dimension(:,:) :: fz   !< Bz or Ez
   sll_real64, intent(in)     :: dt   !< time step

   call faraday_2d_fdtd(self, fx, fy, fz, 0.5*dt)   
   call bc_periodic_2d_fdtd(self, fx, fy, fz, dt)
   call ampere_2d_fdtd(self, fx, fy, fz, dt) 
   call bc_periodic_2d_fdtd(self, fx, fy, fz, dt)
   call faraday_2d_fdtd(self, fx, fy, fz, 0.5*dt)   
   call bc_periodic_2d_fdtd(self, fx, fy, fz, dt)

end subroutine solve_maxwell_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Solve Faraday equation
subroutine faraday_2d_fdtd( self, fx, fy, fz, dt )

type(sll_t_maxwell_2d_fdtd)                 :: self !< Maxwell object
sll_real64, dimension(:,:), target :: fx   !< Ex or Bx
sll_real64, dimension(:,:), target :: fy   !< Ey or By
sll_real64, dimension(:,:), target :: fz   !< Bz or Ez
sll_real64, intent(in) :: dt               !< time step
sll_int32 :: i, j
sll_real64 :: dex_dy, dey_dx, dez_dx, dez_dy
sll_real64 :: dx, dy
sll_int32  :: i1, j1, i2, j2
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: ez
sll_real64, dimension(:,:), pointer :: bx
sll_real64, dimension(:,:), pointer :: by
sll_real64, dimension(:,:), pointer :: bz

dx  = self%dx
dy  = self%dy

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 

i1 = self%i1
j1 = self%j1
i2 = self%i2
j2 = self%j2

if (self%polarization == TE_POLARIZATION) then

   ex => fx; ey => fy; bz => fz
   do i=i1,j1-1
   do j=i2,j2-1
      dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
      dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
      bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
   end do
   end do

end if

if (self%polarization == TM_POLARIZATION) then

   bx => fx; by => fy; ez => fz
   do i=i1,j1
   do j=i2+1,j2
      dez_dy  = (ez(i,j)-ez(i,j-1)) / dy
      bx(i,j) = bx(i,j) - dt * dez_dy 
   end do
   end do
   
   do i=i1+1,j1
   do j=i2,j2
      dez_dx  = (ez(i,j)-ez(i-1,j)) / dx
      by(i,j) = by(i,j) + dt * dez_dx 
   end do
   end do

end if

end subroutine faraday_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Solve ampere-maxwell equation with FDTD scheme
subroutine ampere_2d_fdtd( self, fx, fy, fz, dt, jx, jy )

type(sll_t_maxwell_2d_fdtd) :: self !< Maxwell object
sll_int32 :: i1, j1, i2, j2
sll_real64, dimension(:,:), intent(inout), target :: fx !< Ex or Bx
sll_real64, dimension(:,:), intent(inout), target :: fy !< Ey or By
sll_real64, dimension(:,:), intent(inout), target :: fz !< Bz or Ez
sll_real64, dimension(:,:), intent(in), optional :: jx  !< Jx current
sll_real64, dimension(:,:), intent(in), optional :: jy  !< Jy current
sll_int32 :: i, j
sll_real64 :: dbz_dx, dbz_dy, dbx_dy, dby_dx
sll_real64, intent(in) :: dt !< time step
sll_real64 :: dx, dy
sll_real64 :: csq
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: ez
sll_real64, dimension(:,:), pointer :: bx
sll_real64, dimension(:,:), pointer :: by
sll_real64, dimension(:,:), pointer :: bz

i1 = self%i1
j1 = self%j1
i2 = self%i2
j2 = self%j2

csq = self%c * self%c
dx  = self%dx
dy  = self%dy

!*** Calcul du champ electrique E au temps n+1
!*** sur les points internes du maillage
!*** Ex aux points (i+1/2,j)
!*** Ey aux points (i,j+1/2)

if (self%polarization == TE_POLARIZATION) then

   ex => fx; ey => fy; bz => fz

   do i=i1,j1
   do j=i2+1,j2
      dbz_dy  = (bz(i,j)-bz(i,j-1)) / dy
      ex(i,j) = ex(i,j) + dt*csq*dbz_dy 
   end do
   end do
   
   do i=i1+1,j1
   do j=i2,j2
      dbz_dx  = (bz(i,j)-bz(i-1,j)) / dx
      ey(i,j) = ey(i,j) - dt*csq*dbz_dx 
   end do
   end do

   if (present(jx) .and. present(jy)) then

      ex(i1:j1,i2+1:j2) = ex(i1:j1,i2+1:j2) - dt * jx(i1:j1,i2+1:j2) / self%e_0
      ey(i1+1:j1,i2:j2) = ey(i1+1:j1,i2:j2) - dt * jy(i1+1:j1,i2:j2) / self%e_0

   endif

end if

if (self%polarization == TM_POLARIZATION) then

   bx => fx; by => fy; ez => fz

   do i=i1,j1-1
   do j=i2,j2-1
      dbx_dy  = (bx(i,j+1)-bx(i,j)) / dy
      dby_dx  = (by(i+1,j)-by(i,j)) / dx
      ez(i,j) = ez(i,j) - dt * (dbx_dy - dby_dx)
   end do
   end do

   if (present(jx) .and. .not. present(jy)) then

      ez(i1:j1-1,i2:j2-1) = ez(i1:j1-1,i2:j2-1) - dt * jx(i1:j1-1,i2:j2-1) / self%e_0

   endif

end if

end subroutine ampere_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Set boundary conditions 
subroutine bc_periodic_2d_fdtd(self, fx, fy, fz, dt)

type(sll_t_maxwell_2d_fdtd) :: self !< maxwell solver object
sll_int32 :: i1, j1, i2, j2
sll_real64, dimension(:,:), intent(inout) :: fx !< Ex or Bx
sll_real64, dimension(:,:), intent(inout) :: fy !< Ey or By
sll_real64, dimension(:,:), intent(inout) :: fz !< Bz or Ez
sll_real64 :: dbz_dx, dbz_dy, dez_dx, dez_dy
sll_real64, intent(in) :: dt !< time step
sll_real64 :: dx, dy
sll_int32 :: i, j
sll_real64 :: csq

i1 = self%i1
j1 = self%j1
i2 = self%i2
j2 = self%j2

csq = self%c * self%c
dx  = self%dx
dy  = self%dy

fz(i1:j1-1,j2) = fz(i1:j1-1,i2)
fz(j1,i2:j2-1) = fz(i1,i2:j2-1)

fz(j1,j2) = fz(i1,i2)

if ( self%polarization == TE_POLARIZATION) then
   do i = i1, j1
      dbz_dy = (fz(i,i2)-fz(i,j2-1)) / dy
      fx(i,i2) = fx(i,j2) + dt*csq*dbz_dy 
   end do
     
   do j = i2, j2
      dbz_dx = (fz(i1,j)-fz(j1-1,j)) / dx
      fy(i1,j) = fy(j1,j) - dt*csq*dbz_dx 
   end do
end if

if ( self%polarization == TM_POLARIZATION) then
   do i = i1, j1
      dez_dy = (fz(i,i2)-fz(i,j2-1)) / dy
      fx(i,i2) = fx(i,j2) - dt*csq*dez_dy 
   end do
     
   do j = i2, j2
      dez_dx = (fz(i1,j)-fz(j1-1,j)) / dx
      fy(i1,j) = fy(j1,j) + dt*csq*dez_dx 
   end do
end if

end subroutine bc_periodic_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!PN DEFINED BUT NOT USED
!PN !> Set periodic bounday conditions
!PN subroutine bc_metallic_2d_fdtd(self, fx, fy, fz, side)
!PN 
!PN type(sll_t_maxwell_2d_fdtd) :: self !< maxwell object
!PN sll_int32, intent(in) :: side !< which domain edge
!PN sll_real64, dimension(:,:), target  :: fx !< Ex or Bx
!PN sll_real64, dimension(:,:), target  :: fy !< Ey or By
!PN sll_real64, dimension(:,:), target  :: fz !< Bz or Ez
!PN sll_real64, dimension(:,:), pointer :: bx, by, ez
!PN sll_real64, dimension(:,:), pointer :: ex, ey, bz
!PN sll_int32 :: i1, j1, i2, j2
!PN 
!PN i1 = self%i1
!PN j1 = self%j1
!PN i2 = self%i2
!PN j2 = self%j2
!PN 
!PN if (self%polarization == TE_POLARIZATION) then
!PN    ex => fx; ey => fy; bz => fz
!PN    select case(side)
!PN    case(SOUTH)
!PN       ex(i1:j1,i2) = 0._f64
!PN    case(NORTH)
!PN       bz(i1:j1,j2) = bz(i1:j1,j2-1)
!PN    case(WEST)
!PN       ey(i1,i2:j2) = 0._f64
!PN    case(EAST)
!PN       bz(j1,i2:j2) = bz(j1-1,i2:j2)
!PN    end select
!PN end if
!PN 
!PN if (self%polarization == TM_POLARIZATION) then
!PN    bx => fx; by => fy; ez => fz
!PN    select case(side)
!PN    case(SOUTH)
!PN       bx(i1:j1,i2) = 0.0_f64
!PN    case(NORTH)
!PN       ez(i1:j1,j2) = ez(i1:j1,j2-1)
!PN    case(WEST)
!PN       by(i1,i2:j2) = 0._f64
!PN    case(EAST)
!PN       ez(j1,i2:j2) = ez(j1-1,i2:j2)
!PN    end select
!PN end if
!PN 
!PN end subroutine bc_metallic_2d_fdtd
!PN 
!PN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PN 
!PN !> Bundary conditions
!PN subroutine bc_silver_muller_2d_fdtd( self, ex, ey, bz, ccall, dt )
!PN 
!PN type(sll_t_maxwell_2d_fdtd) :: self !< maxwell object
!PN sll_int32, intent(in) :: ccall !< domain edge (N,S,E,W)
!PN sll_int32 :: i1, j1, i2, j2
!PN sll_real64 :: a11,a12,a21,a22,b1,b2,dis
!PN sll_int32 :: i, j
!PN sll_real64, intent(in) :: dt !< time step
!PN sll_real64 :: dx, dy, c, csq
!PN sll_real64, dimension(:,:), pointer :: ex !< x electric field
!PN sll_real64, dimension(:,:), pointer :: ey !< y electric field
!PN sll_real64, dimension(:,:), pointer :: bz !< z magnetic field
!PN 
!PN c   = self%c
!PN csq = c * c
!PN dx  = self%dx
!PN dy  = self%dy
!PN 
!PN i1 = self%i1
!PN j1 = self%j1
!PN i2 = self%i2
!PN j2 = self%j2
!PN 
!PN !Conditions de Silver-Muller
!PN !------------------------------------
!PN !Ey = -c Bz sur la frontiere ouest
!PN !Ey =  c Bz sur la frontiere est
!PN !Ex = -c Bz sur la frontiere nord
!PN !Ex =  c Bz sur la frontiere sud
!PN    
!PN !On effectue le calcul de B sur les points fictifs du maillage
!PN !simultanement avec la prise en compte des conditions limites sur
!PN !E. Pour cela, on resout sur les points frontieres, l'equation de la
!PN !condition limite en moyennant en temps pour E et en espace pour B puis
!PN !l'equation d'Ampere
!PN 
!PN select case (ccall)
!PN 
!PN case (NORTH)
!PN    !Frontiere Nord : Ex = -c Bz 
!PN    do i = i1, j1
!PN          
!PN       a11 = 1.0_f64; a12 = + c
!PN       a21 = 1.0_f64/dt; a22 = - csq / dy
!PN       b1  = - ex(i,j2) - c * bz(i,j2-1)
!PN       b2  =   ex(i,j2)/dt - csq/dy*bz(i,j2-1)
!PN          
!PN       dis = a11*a22-a21*a12 
!PN          
!PN       !ex(i,j2) = (b1*a22-b2*a12)/dis
!PN       bz(i,j2) = (a11*b2-a21*b1)/dis
!PN          
!PN    end do
!PN       
!PN case (SOUTH)
!PN 
!PN    !Frontiere Sud : Ex =  c Bz
!PN    do i = i1, j1
!PN          
!PN       a11 = 1.0_f64; a12 = - c
!PN       a21 = 1.0_f64/dt; a22 = csq / dy
!PN       b1  = - ex(i,i2) + c * bz(i,i2+1)
!PN       b2  = ex(i,i2)/dt + csq / dy * bz(i,i2+1) 
!PN          
!PN       dis = a11*a22-a21*a12 
!PN          
!PN       ex(i,i2) = (b1*a22-b2*a12)/dis
!PN       !bz(i,i2) = (a11*b2-a21*b1)/dis
!PN          
!PN    end do
!PN       
!PN case (EAST)
!PN 
!PN    !Frontiere Est : Ey =  c Bz
!PN    do j = i2, j2
!PN          
!PN       a11 = 1.0_f64; a12 = - c
!PN       a21 = 1.0_f64/dt; a22 = + csq / dx
!PN       b1  = - ey(j1,j) + c * bz(j1-1,j)
!PN       b2  = ey(j1,j)/dt + csq/dx*bz(j1-1,j) 
!PN          
!PN       dis = a11*a22-a21*a12 
!PN          
!PN       !ey(j1,j) = (b1*a22-b2*a12)/dis
!PN       bz(j1,j) = (a11*b2-a21*b1)/dis
!PN       
!PN    end do
!PN       
!PN case (WEST)
!PN 
!PN    !Frontiere Ouest : Ey = -c Bz
!PN    do j = i2, j2
!PN       
!PN       a11 = 1.0_f64; a12 = + c
!PN       a21 = 1.0_f64/dt; a22 = - csq / dx
!PN       b1  = - ey(i1,j) - c * bz(i1+1,j)
!PN       b2  =   ey(i1,j)/dt - csq/dx*bz(i1+1,j) 
!PN       
!PN       dis = a11*a22-a21*a12 
!PN    
!PN       ey(i1,j) = (b1*a22-b2*a12)/dis
!PN       !bz(i1,j) = (a11*b2-a21*b1)/dis
!PN       
!PN    end do
!PN 
!PN end select
!PN 
!PN end subroutine bc_silver_muller_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_m_maxwell_2d_fdtd
