module sll_vlasov4d

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

 use used_precision
 use geometry_module
 use diagnostiques_module
 use sll_cubic_splines
 use sll_cubic_spline_interpolator_1d
 use sll_quintic_spline_interpolator_1d
 use sll_cubic_spline_interpolator_2d

 implicit none
 private
 public :: new, dealloc, densite_charge,transposexv,transposevx,densite_courant
 public :: advection_x1, advection_x2, advection_x3, advection_x4
 public :: thdiag

 type, public :: vlasov2d
   sll_real64, dimension(:,:,:,:), pointer :: ft
   type (geometry) :: geomx, geomv
   logical :: transposed      
   sll_int32 :: jstartx, jendx
   sll_int32 :: jstartv, jendv
   class(sll_interpolator_1d_base), pointer :: interp_x1
   class(sll_interpolator_1d_base), pointer :: interp_x2
   class(sll_interpolator_1d_base), pointer :: interp_x3
   class(sll_interpolator_1d_base), pointer :: interp_x4
#ifdef _TWO_D
   class(sll_interpolator_2d_base), pointer :: interp_x
   class(sll_interpolator_2d_base), pointer :: interp_v
#endif
 end type vlasov2d

#ifdef _QUINTIC
 type(cubic_spline_1d_interpolator),   target :: spl_x1
 type(cubic_spline_1d_interpolator),   target :: spl_x2
 type(quintic_spline_1d_interpolator), target :: spl_x3
 type(quintic_spline_1d_interpolator), target :: spl_x4
#else
 type(cubic_spline_1d_interpolator), target :: spl_x1
 type(cubic_spline_1d_interpolator), target :: spl_x2
 type(cubic_spline_1d_interpolator), target :: spl_x3
 type(cubic_spline_1d_interpolator), target :: spl_x4
#endif

#ifdef _TWO_D
 type(cubic_spline_2d_interpolator),   target :: spl_x
 type(cubic_spline_2d_interpolator),   target :: spl_v
#endif

 sll_int32, private :: i, j, k, l

 interface new
   module procedure initialize_vlasov2d
 end interface

 interface dealloc
   module procedure dealloc_vlasov2d
 end interface

contains

 subroutine initialize_vlasov2d(this,geomx,geomv,error, jstartx, jendx, &
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

#ifdef TWO_D
  call this%interp_x%initialize( nc_x1, nc_x2, &
                                    x1_min, x1_max, x2_min, x2_max, &
                                    PERIODIC_SPLINE, PERIODIC_SPLINE)
  call this%interp_v%initialize( nc_x3, nc_x4, &
                                    x3_min, x3_max, x4_min, x4_max, &
                                    PERIODIC_SPLINE, PERIODIC_SPLINE)
  this%interp_x => spl_x
  this%interp_v => spl_v
#else

  call spl_x1%initialize( nc_x1, x1_min, x1_max, SLL_PERIODIC)
  call spl_x2%initialize( nc_x2, x2_min, x2_max, SLL_PERIODIC)
  call spl_x3%initialize( nc_x3, x3_min, x3_max, SLL_PERIODIC)
  call spl_x4%initialize( nc_x4, x4_min, x4_max, SLL_PERIODIC)

  this%interp_x1 => spl_x1
  this%interp_x2 => spl_x2
  this%interp_x3 => spl_x3
  this%interp_x4 => spl_x4

#endif

 end subroutine initialize_vlasov2d

 subroutine dealloc_vlasov2d(this)
  type(vlasov2d),intent(out) :: this
  sll_int32 :: error
  SLL_DEALLOCATE(this%ft,error)
 end subroutine dealloc_vlasov2d

 subroutine advection_x1(this,f,dt)
  type(vlasov2d) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x2, nc_x3
  sll_real64 :: alpha, x3_min, delta_x3

  ! verifier que la transposition est a jours
  SLL_ASSERT( .not. this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x2    = this%geomx%ny
  nc_x3    = this%geomv%nx

  x3_min   = this%geomv%x0
  delta_x3 = this%geomv%dx

  do l=this%jstartv,this%jendv
     do k=1,nc_x3
        alpha = (x3_min +(k-1)*delta_x3)*dt
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
  sll_real64 :: x4_min, delta_x4
  sll_real64 :: alpha

  SLL_ASSERT( .not. this%transposed)

  nc_x1    = this%geomx%nx
  nc_x2    = this%geomx%ny
  nc_x3    = this%geomv%nx

  x4_min   = this%geomv%y0
  delta_x4 = this%geomv%dy

  do l=this%jstartv,this%jendv
    alpha = (x4_min +(l-1)*delta_x4)*dt
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
#ifdef _QUINTIC
  sll_real64 :: depvx(this%geomv%nx)
  sll_real64 :: delta_x3
  sll_real64 :: x3_min
  sll_real64 :: x3_max
#endif

  SLL_ASSERT(this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny

#ifdef _QUINTIC

  delta_x3 = this%geomv%dx

  x3_min = this%geomv%x0
  x3_max = this%geomv%x1

  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        alpha = ex(i,j)*dt
        do k=1,nc_x3
           depvx(k) = x3_min + (k-1)*delta_x3 - alpha
           if(depvx(k) < x3_min) then
              depvx(k) = depvx(k) + x3_max-x3_min
           else if (depvx(k) > x3_max) then
              depvx(k) = depvx(k) - x3_max+x3_min
           end if
        end do 
        do l=1,nc_x4
           this%ft(:,l,i,j) = this%interp_x3%interpolate_array( &
                               nc_x3, this%ft(:,l,i,j), depvx)
        end do
      end do
   end do
#else
  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        alpha = ex(i,j)*dt
        do l=1,nc_x4
           this%ft(:,l,i,j) = this%interp_x3%interpolate_array_disp( &
                              nc_x3, this%ft(:,l,i,j), alpha )
       end do
    end do
 end do
#endif

 end subroutine advection_x3

 subroutine advection_x4(this,ey,dt)

  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:), intent(in) :: ey
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x3, nc_x4
  sll_real64 :: alpha
#ifdef _QUINTIC
  sll_real64 :: depvy(this%geomv%ny)
  sll_real64 :: delta_x4
  sll_real64 :: x4_min
  sll_real64 :: x4_max
#endif

  SLL_ASSERT(this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny

#ifdef _QUINTIC

  delta_x4 = this%geomv%dy

  x4_min = this%geomv%y0
  x4_max = this%geomv%y1

  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        alpha = ey(i,j)*dt
        do l=1,nc_x4
           depvy(l) = x4_min + (l-1)*delta_x4 - alpha
           if(depvy(l) < x4_min) then
              depvy(l) = depvy(l) + x4_max-x4_min
           else if (depvy(l) > x4_max) then
              depvy(l) = depvy(l) - x4_max+x4_min
           end if
        end do 
        do k=1,nc_x3
           this%ft(k,:,i,j) = this%interp_x4%interpolate_array( &
                               nc_x4, this%ft(k,:,i,j), depvy)
        end do
      end do
   end do
#else

  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        do k=1,nc_x3
           alpha = ey(i,j)*dt
           this%ft(k,:,i,j) = this%interp_x4%interpolate_array_disp( &
                              nc_x4, this%ft(k,:,i,j), alpha )
       end do
    end do
 end do

#endif

 end subroutine advection_x4

#ifdef TWO_D

 subroutine advection_x(this,f,dt)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x2
  sll_real64 :: depx, depy 

  SLL_ASSERT(.not. this%transposed)

  nc_x1    = this%geomx%nx
  nc_x2    = this%geomx%ny

  do l=1,this%geomv%nx
     depy = (this%geomv%y0+(l-1)*this%geomv%dy)*dt
     do k=this%jstartv,this%jendv
        depx = (this%geomv%x0+(k-1)*this%geomv%dx)*dt
        f(:,:,k,l) = this%interp_x%interpolate_array_disp( &
                     nc_x1, nc_x2, f(:,:,k,l), depx, depy )
     end do
  end do

 end subroutine advection_x

 subroutine advection_v(this,fx,fy,dt)

  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:), intent(in) :: fx, fy
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x3, nc_x4
  sll_real64 :: depvx, depvy  

  SLL_ASSERT(this%transposed)

  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny

  do j=this%jstartx,this%jendx
    do i=1,this%geomx%nx
      depvx = fx(i,j)*dt
      depvy = fy(i,j)*dt
      this%ft(:,:,i,j) = this%interp_v%interpolate_array_disp( &
                         nc_x3, nc_x4, this%ft(:,:,i,j), depx, depy )
    end do
   end do

 end subroutine advection_v

#endif

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

subroutine densite_courant(this, jx, jy)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: jx, jy
   sll_real64 :: vx, vy 
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

subroutine thdiag(this,f,nrj,t,jstartv)

   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,jstartv:),intent(in) :: f
   sll_int32 :: jstartv
   !sll_int32 :: error
   sll_real64, intent(in) :: t,nrj   ! current time
   ! variables locales
   !sll_int32 :: i,iv, j,jv
   !sll_real64 :: x, vx, y, vy
   !sll_real64,dimension(7) :: diagloc
   sll_real64,dimension(11) :: auxloc
   sll_int32 :: my_num, num_threads
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: comm, error
   sll_real64 :: cell_volume

   comm   = sll_world_collective%comm
   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   auxloc = 0.0
   cell_volume = this%geomx%dx * this%geomx%dy * this%geomv%dx * this%geomv%dy
   auxloc(1) = cell_volume * sum(f) ! avg(f)
   auxloc(2) = cell_volume * sum(abs(f)) ! L1 norm
   auxloc(3) = cell_volume * sum(f*f) ! L2 norm
   
   call mpi_reduce(auxloc,aux(1:11),11,MPI_REAL8,MPI_SUM,MPI_MASTER,comm, error)

   if (my_num == MPI_MASTER) then
      diag=0.
      aux=0.
      aux(13)=t
      aux(12)=nrj
      write(*,"('time ', g8.3,' test nrj',g15.5)") t, nrj
      call time_history("thf","(13(1x,e15.6))",aux(1:13))
   end if

end subroutine thdiag

end module sll_vlasov4d
