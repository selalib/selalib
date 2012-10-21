module sll_vlasov2d

#include "selalib.h"

 use used_precision
 use geometry_module

 implicit none
 private
 public :: new, dealloc, densite_charge,transposexv,transposevx,densite_courant
 public :: advection_x1, advection_x2, advection_x3, advection_x4

 sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4
 sll_real64 :: x1_min, x2_min, x3_min, x4_min
 sll_real64 :: x1_max, x2_max, x3_max, x4_max
 sll_real64 :: delta_x1, delta_x2, delta_x3, delta_x4

 type, public :: vlasov2d
   sll_real64, dimension(:,:,:,:), pointer :: ft
   type(cubic_spline_1d_interpolator) :: interp_x1
   type(cubic_spline_1d_interpolator) :: interp_x2
   type(cubic_spline_1d_interpolator) :: interp_x3
   type(cubic_spline_1d_interpolator) :: interp_x4
   type (geometry) :: geomx, geomv
   logical :: transposed       ! permet de definir si f ou ft est derniere 
                               ! fonction de distribution mise a jour
   sll_int32 :: jstartx, jendx
   sll_int32 :: jstartv, jendv  ! definition de la bande de calcul
 end type vlasov2d

 sll_int32, private :: i, j, k, l

 interface new
   module procedure new_vlasov2d
 end interface

 interface dealloc
   module procedure dealloc_vlasov2d
 end interface

contains

 subroutine new_vlasov2d(this,geomx,geomv,error, jstartx, jendx, &
                        jstartv, jendv,vz)

  type(vlasov2d),intent(out)      :: this
  type(geometry),intent(in)       :: geomx
  type(geometry),intent(in)       :: geomv
  sll_real64, optional            :: vz
  sll_int32, intent(out)          :: error
  sll_int32, intent(in), optional ::  jstartx
  sll_int32, intent(in), optional ::  jendx
  sll_int32, intent(in), optional ::  jstartv
  sll_int32, intent(in), optional ::  jendv

  error = 0

  ! on commence par utiliser la fonction f(x,y,vx,vy)
  this%transposed=.false.
  if (.not.(present(jstartx))) then
     this%jstartx = 1
  else
     this%jstartx = jstartx
  end if
  if (.not.(present(jendx))) then
     this%jendx = geomx%ny
  else
     this%jendx = jendx
  end if
  if (.not.(present(jstartv))) then 
      this%jstartv = 1
  else
     this%jstartv = jstartv
  end if
  if (.not.(present(jendv))) then
     this%jendv = geomv%nx
  else
     this%jendv = jendv
  end if

  nc_x1  = geomx%nx
  nc_x2  = geomx%ny
  nc_x3  = geomv%nx
  nc_x4  = geomv%ny

  x1_min = geomx%x0
  x1_max = geomx%x1
  x2_min = geomx%y0
  x2_max = geomx%y1
  x3_min = geomv%x0
  x3_max = geomv%x1
  x4_min = geomv%y0
  x4_max = geomv%y1

  delta_x1 = geomx%dx
  delta_x2 = geomx%dy
  delta_x3 = geomv%dx
  delta_x4 = geomv%dy
  
  this%geomx=geomx
  this%geomv=geomv

  SLL_ALLOCATE(this%ft(geomv%nx,geomv%ny,geomx%nx,this%jstartx:this%jendx),error)

  call this%interp_x1%initialize( nc_x1, x1_min, x1_max, PERIODIC_SPLINE)
  call this%interp_x2%initialize( nc_x2, x2_min, x2_max, PERIODIC_SPLINE)
  call this%interp_x3%initialize( nc_x3, x3_min, x3_max, HERMITE_SPLINE)
  call this%interp_x4%initialize( nc_x4, x4_min, x4_max, HERMITE_SPLINE)

 end subroutine new_vlasov2d

 subroutine dealloc_vlasov2d(this)
  type(vlasov2d),intent(out) :: this
  sll_int32 :: error
  SLL_DEALLOCATE(this%ft,error)
 end subroutine dealloc_vlasov2d

 subroutine advection_x1(this,f,dt)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt

  sll_real64 :: depx   ! deplacement par rapport au maillage
  sll_real64 :: vx     ! vitesse du point courant
  sll_int32 :: iv, jv ! indices de boucle

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  do l = 1, nc_x4
  do k = 1, nc_x3 
  do j = 1, nc_x2
  do i = 1, nc_x1
     vx = x3_min+(iv-1)*this%geomv%dx
     depx = vx*dt
  end do
  end do
  end do
  end do

 end subroutine advection_x1

 subroutine advection_x2(this,f,dt)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt

  sll_real64 :: depx, depy   ! deplacement par rapport au maillage
  sll_real64 :: vx, vy       ! vitesse du point courant
  sll_int32 :: iv, jv ! indices de boucle

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed)
  do l = 1, nc_x4
  do k = 1, nc_x3 
  do j = 1, nc_x2
  do i = 1, nc_x1
     vy = this%geomv%y0+(jv-1)*this%geomv%dy
     depy = vy*dt
     !call interpole(this%splinex,f(:,:,iv,jv),depx,depy,jv==0)
  end do
  end do
  end do
  end do

 end subroutine advection_x2

 subroutine advection_x3(this,fx,dt)

  type(vlasov2d),intent(inout) :: this

  sll_real64, dimension(:,:), intent(in) :: fx
  sll_real64, intent(in) :: dt

  sll_real64 :: depvx, depvy   ! deplacement par rapport au maillage
  sll_int32 :: i, j ,iv, jv, im1! indices de boucle

  SLL_ASSERT(this%transposed) 

  do l = 1, nc_x4
  do k = 1, nc_x3 
  do j = 1, nc_x2
     depvx =  fx(i,j)*dt
     !call this%interp_x3%compute_interpolants( f(:,j,k,l) )
     !f(:,j,k,l) = interp_x1%interpolate_array_disp(nc_x1,f(:,j,k,l),depvx ) 
  end do
  end do
  end do

 end subroutine advection_x3

 subroutine advection_x4(this,fy,dt)

  type(vlasov2d),intent(inout) :: this

  sll_real64, dimension(:,:), intent(in) :: fy
  sll_real64, intent(in) :: dt

  sll_real64 :: depvx, depvy   ! deplacement par rapport au maillage
  sll_int32 :: i, j ,iv, jv, im1! indices de boucle

  SLL_ASSERT(this%transposed) 

  do l = 1, nc_x4
  do k = 1, nc_x3 
  do i = 1, nc_x1
     depvy = fy(i,j)*dt
     !call this%interp_x3%compute_interpolants( f(:,j,k,l) )
     !f(i,:,k,l) = interp_x1%interpolate_array_disp(nc_x1,f(:,j,k,l),depvx ) 
  end do
  end do
  end do

 end subroutine advection_x4

!>------------------------------------------------
!> calcule la densite de charge rho a partir de ft
!> en fait le moment d'ordre 0 de f. Les constantes
!> ne sont pas introduites
!>------------------------------------------------
!> Poisson n'est pas parallele on transmet donc rho
!> a tous les processeurs
subroutine densite_charge(this, rho)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: rho
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locrho
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   SLL_ASSERT(this%transposed)
   
   rho(:,:) = 0.
   locrho(:,:) = 0.
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do jv=1,this%geomv%ny-1
            do iv=1,this%geomv%nx-1 
               locrho(i,j) = locrho(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) 
            end do
         end do
      end do
   end do
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locrho,rho,c,MPI_REAL8,MPI_SUM,comm,error)

end subroutine densite_charge

!>------------------------------------------------
!> calcule la densite de courant jx et jy a partir de ft
!> en fait le moment d'ordre 0 de f. Les constantes
!> ne sont pas introduites
!>------------------------------------------------
subroutine densite_courant(this, jx, jy)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: jx, jy
   sll_real64 :: vx, vy       ! vitesse du point courant
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjx
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjy
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   ! verifier que la transposition est a jour
   if (.not.(this%transposed)) &
      stop 'densite_courant: on travaille sur ft et pas f'
   !    rho(:,this%jstartx:this%jendx)=0.
   jx(:,:) = 0.; jy(:,:) = 0.
   locjx(:,:) = 0.; locjy(:,:) = 0.
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do jv=1,this%geomv%ny-1
            vy = this%geomv%y0+(jv-1)*this%geomv%dy
            do iv=1,this%geomv%nx-1 
               vx = this%geomv%x0+(iv-1)*this%geomv%dx
               locjx(i,j) = locjx(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) * vx
               locjy(i,j) = locjy(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) * vy
            end do
         end do
      end do
   end do

   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjx,jx,c, MPI_REAL8,MPI_SUM,comm,error)
   call mpi_allreduce(locjy,jy,c, MPI_REAL8,MPI_SUM,comm,error)

end subroutine densite_courant

!>---------------------------------------------
!> transpose la fonction de distribution
!> ATTENTION: cet fonction fait intervenir des 
!> communications entre les processeurs.
!>---------------------------------------------
subroutine transposexv(this,f)
   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,:),intent(in) :: f
   sll_int32 :: sizexy, sizevxvy
   sll_int32 :: my_num, num_threads

   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   sizexy = this%geomx%nx * this%geomx%ny
   sizevxvy = this%geomv%nx * this%geomv%ny
   call transpose(f, this%ft, sizexy, sizevxvy, num_threads)

this%transposed=.true.

end subroutine transposexv

subroutine transposevx(this,f)

   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,:),intent(out) :: f

   sll_int32 :: sizexy, sizevxvy
   sll_int32 :: my_num, num_threads

   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   sizexy = this%geomx%nx * this%geomx%ny
   sizevxvy = this%geomv%nx * this%geomv%ny
   call transpose(this%ft,f, sizevxvy, sizexy, num_threads)

   this%transposed=.false.

end subroutine transposevx


!      do j=1,loc_sz_x2
!         do i=1,loc_sz_x1
!            do l=1,sim%mesh4d%num_cells4
!               alpha = -sim%efield_split(i,j)%ex*0.5_f64*sim%dt
!               ! interpolate_array_disp() has an interface that must be changed
!               sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
!                    sim%nc_x3, &
!                    sim%f_x3x4(i,j,:,l), &
!                    alpha )
!            end do
!         end do
!      end do

!      ! Continue with vy...(x4)
!      do j=1,loc_sz_x2
!         do i=1,loc_sz_x1
!            do k=1,sim%mesh4d%num_cells3
!               alpha = -sim%efield_split(i,j)%ey*0.5_f64*sim%dt
!               ! interpolate_array_disp() has an interface that must be changed
!               sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
!                    sim%nc_x4, &
!                    sim%f_x3x4(i,j,k,:), &
!                    alpha )
!            end do
!         end do
!      end do

!      ! full time step advection in 'x' (x1)
!      do l=1,loc_sz_x4
!         do k=1,loc_sz_x3
!            do j=1,sim%mesh4d%num_cells2
!               vmin = sim%mesh4d%x3_min
!               delta = sim%mesh4d%delta_x3
!               alpha = (vmin + (j-1)*delta)*sim%dt
!               call sim%interp_x1%compute_interpolants( sim%f_x1x2(:,j,k,l) )
!               ! interpolate_array_disp() has an interface that must be changed
!               sim%f_x1x2(:,j,k,l) = sim%interp_x1%interpolate_array_disp( &
!                    sim%nc_x1, &
!                    sim%f_x1x2(:,j,k,l), &
!                    alpha )
!            end do
!         end do
!      end do
!      
!      ! full time step advection in 'y' (x2)
!      do l=1,loc_sz_x4
!         do k=1,loc_sz_x3
!            do i=1,sim%mesh4d%num_cells1
!               vmin = sim%mesh4d%x4_min
!               delta = sim%mesh4d%delta_x4
!               alpha = (vmin + (j-1)*delta)*sim%dt
!               call sim%interp_x1%compute_interpolants( sim%f_x1x2(i,:,k,l) )
!               ! interpolate_array_disp() has an interface that must be changed
!               sim%f_x1x2(i,:,k,l) = sim%interp_x2%interpolate_array_disp( &
!                    sim%nc_x2, &
!                    sim%f_x1x2(i,:,k,l), &
!                    alpha )
!            end do
!         end do
!      end do
!      


end module sll_vlasov2d
