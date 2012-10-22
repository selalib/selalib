module sll_vlasov2d

#include "selalib.h"

 use used_precision
 use geometry_module
 use diagnostiques_module

 implicit none
 private
 public :: new, dealloc, densite_charge,transposexv,transposevx,densite_courant
 public :: advection_x1, advection_x2, advection_x3, advection_x4
 public :: thdiag

 type, public :: vlasov2d
   sll_real64, dimension(:,:,:,:), pointer :: ft
   type(cubic_spline_1d_interpolator) :: interp_x1
   type(cubic_spline_1d_interpolator) :: interp_x2
   type(cubic_spline_1d_interpolator) :: interp_x3
   type(cubic_spline_1d_interpolator) :: interp_x4
   type (geometry) :: geomx, geomv
   logical :: transposed      
   sll_int32 :: jstartx, jendx
   sll_int32 :: jstartv, jendv
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
                        jstartv, jendv)

  type(vlasov2d),intent(out)      :: this
  type(geometry),intent(in)       :: geomx
  type(geometry),intent(in)       :: geomv
  sll_int32, intent(out)          :: error
  sll_int32, intent(in)           :: jstartx
  sll_int32, intent(in)           :: jendx
  sll_int32, intent(in)           :: jstartv
  sll_int32, intent(in)           :: jendv

  sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_real64 :: x1_min, x2_min, x3_min, x4_min
  sll_real64 :: x1_max, x2_max, x3_max, x4_max

  error = 0

  this%transposed=.false.

  this%jstartx = jstartx
  this%jendx   = jendx
  this%jstartv = jstartv
  this%jendv   = jendv

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

  this%geomx=geomx
  this%geomv=geomv

  SLL_ALLOCATE(this%ft(geomv%nx,geomv%ny,geomx%nx,this%jstartx:this%jendx),error)

  call this%interp_x1%initialize( nc_x1, x1_min, x1_max, PERIODIC_SPLINE)
  call this%interp_x2%initialize( nc_x2, x2_min, x2_max, PERIODIC_SPLINE)
  call this%interp_x3%initialize( nc_x3, x3_min, x3_max, PERIODIC_SPLINE)
  call this%interp_x4%initialize( nc_x4, x4_min, x4_max, PERIODIC_SPLINE)

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
  sll_int32  :: nc_x1, nc_x2, nc_x3
  sll_real64 :: alpha

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x2    = this%geomx%ny
  nc_x3    = this%geomv%nx

  do l=this%jstartv,this%jendv
     do k=1,nc_x3
        alpha = this%geomv%xgrid(k)*dt
        do j=1,nc_x2
           f(:,j,k,l) = this%interp_x1%interpolate_array_disp( &
                                        nc_x1, f(:,j,k,l), alpha )
        end do
     end do
  end do

 end subroutine advection_x1

 subroutine advection_x2(this,f,dt)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x2, nc_x3
  sll_real64 :: alpha

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed)

  nc_x1    = this%geomx%nx
  nc_x2    = this%geomx%ny
  nc_x3    = this%geomv%nx

  do l=this%jstartv,this%jendv
    alpha = this%geomv%ygrid(l)*dt
    do k=1,nc_x3
        do i=1,nc_x1
           f(i,:,k,l) = this%interp_x2%interpolate_array_disp( &
                    nc_x2, f(i,:,k,l), alpha )
       end do
    end do
  end do

 end subroutine advection_x2

 subroutine advection_x3(this,ex,dt)

  type(vlasov2d),intent(inout) :: this

  sll_real64, dimension(:,:), intent(in) :: ex
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x3, nc_x4
  sll_real64 :: alpha

  SLL_ASSERT(this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny

  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        do l=1,nc_x4
           alpha = ex(i,j)*dt
           this%ft(:,l,i,j) = this%interp_x3%interpolate_array_disp( &
                              nc_x3, this%ft(:,l,i,j), alpha )
       end do
    end do
 end do

 end subroutine advection_x3

 subroutine advection_x4(this,ey,dt)

  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:), intent(in) :: ey
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x3, nc_x4
  sll_real64 :: alpha

  SLL_ASSERT(this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny

  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        do k=1,nc_x3
           alpha = ey(i,j)*dt
           this%ft(k,:,i,j) = this%interp_x4%interpolate_array_disp( &
                              nc_x4, this%ft(k,:,i,j), alpha )
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

   SLL_ASSERT(this%transposed)

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

subroutine thdiag(this,f,nrj,t)

   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,this%jstartv:),intent(in) :: f
   !sll_int32 :: error
   sll_real64, intent(in) :: t,nrj   ! current time
   ! variables locales
   !sll_int32 :: i,iv, j,jv
   !sll_real64 :: x, vx, y, vy
   !sll_real64,dimension(7) :: diagloc
   !sll_real64,dimension(11) :: auxloc
   sll_int32 :: my_num, num_threads
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: comm

   comm   = sll_world_collective%comm
   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   if (my_num == MPI_MASTER) then
      diag=0.
      aux=0.
   end if

!   diagloc = 0._wp
!   auxloc  = 0._wp
!   do i = 1,this%geomx%nx
!      x = this%geomx%x0+(i-1)*this%geomx%dx
!      do j = 1,this%geomx%ny
!         y= this%geomx%y0+(j-1)*this%geomx%dy
!         do iv=1,this%geomv%nx
!            vx = this%geomv%x0+(iv-1)*this%geomv%dx
!            do jv=this%jstartv,this%jendv
!               vy = this%geomv%y0+(jv-1)*this%geomv%dy
!               diagloc(2) = diagloc(2) + f(i,j,iv,jv)*f(i,j,iv,jv)
!
!               auxloc(1) = auxloc(1) + f(i,j,iv,jv)         ! avg(f)
!               auxloc(2) = auxloc(2) + x*f(i,j,iv,jv)       ! avg(x)
!               auxloc(3) = auxloc(3) + vx*f(i,j,iv,jv)      ! avg(vx)
!               auxloc(4) = auxloc(4) + x*x*f(i,j,iv,jv)     ! avg(x^2)
!               auxloc(5) = auxloc(5) + vx*vx*f(i,j,iv,jv)   ! avg(vx^2)
!               auxloc(6) = auxloc(6) + x*vx*f(i,j,iv,jv)    ! avg(x*vx)
!               auxloc(7) = auxloc(7) + y*f(i,j,iv,jv)       ! avg(y)
!               auxloc(8) = auxloc(8) + vy*f(i,j,iv,jv)      ! avg(vy)
!               auxloc(9) = auxloc(9) + y*y*f(i,j,iv,jv)     ! avg(y^2)
!               auxloc(10) = auxloc(10) + vy*vy*f(i,j,iv,jv) ! avg(vy^2)
!               auxloc(11) = auxloc(11) + y*vy*f(i,j,iv,jv)  ! avg(y*vy)
!
!            end do
!         end do
!      end do
!   end do
!   auxloc=auxloc!*this%geomx%dx*this%geomx%dy*this%geomv%dx*this%geomv%dy
!
!   call mpi_reduce(auxloc,aux,11,MPI_REAL8,MPI_SUM,0,  &
!                   comm, error)
!   call mpi_reduce(diagloc(2),diag(2),1,MPI_REAL8,MPI_SUM,0, &
!                   comm, error)
!
if (my_num==MPI_MASTER) then
   aux(13)=t
   aux(12)=nrj
   write(*,"('time ', g8.3,' test nrj',f10.5)") t, nrj
   call time_history("thf","(13(1x,e15.6))",aux(1:13))
end if

end subroutine thdiag

end module sll_vlasov2d
