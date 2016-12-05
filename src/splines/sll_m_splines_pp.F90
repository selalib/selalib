
!> @ingroup particle_groups
!> @author Benedikt Perse, Katharina Kormann, IPP
!> @brief Splines in pp form
!> @details Only implements splines on a uniform grid
module sll_m_splines_pp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none 

  public :: &
       sll_t_spline_pp_1d, &
       sll_t_spline_pp_2d, &
       sll_t_spline_pp_3d, &
       sll_s_spline_pp_b_to_pp_1d, &
       sll_s_spline_pp_b_to_pp_2d, &
       sll_s_spline_pp_b_to_pp_3d, &
       sll_s_spline_pp_horner_m_1d, &
       sll_s_spline_pp_horner_m_2d, &
       sll_s_spline_pp_horner_m_3d, &
       sll_s_spline_pp_horner_primitive_1d, &
       sll_f_spline_pp_horner_1d, &
       sll_f_spline_pp_horner_2d, &
       sll_f_spline_pp_horner_3d, &
       sll_s_spline_pp_init_1d, &
       sll_s_spline_pp_init_2d, &
       sll_s_spline_pp_init_3d, &
       sll_s_spline_pp_free_1d, &
       sll_s_spline_pp_free_2d, &
       sll_s_spline_pp_free_3d
  
  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  sll_real64, parameter :: inv_2  = 1._f64/2._f64
  sll_real64, parameter :: inv_3  = 1._f64/3._f64
  sll_real64, parameter :: inv_4  = 1._f64/4._f64
  sll_real64, parameter :: inv_6  = 1._f64/6._f64
  sll_real64, parameter :: inv_8  = 1._f64/8._f64
  sll_real64, parameter :: inv_12 = 1._f64/12._f64
  sll_real64, parameter :: inv_18 = 1._f64/18._f64
  sll_real64, parameter :: inv_20 = 1._f64/20._f64
  sll_real64, parameter :: inv_24 = 1._f64/24._f64
  sll_real64, parameter :: inv_30 = 1._f64/30._f64
  sll_real64, parameter :: inv_36 = 1._f64/36._f64
  sll_real64, parameter :: inv_48 = 1._f64/48._f64
  sll_real64, parameter :: inv_72 = 1._f64/72._f64
  sll_real64, parameter :: inv_120 = 1._f64/120._f64
  sll_real64, parameter :: inv_144 = 1._f64/144._f64
  sll_real64, parameter :: inv_720 = 1._f64/720._f64
  !>arbitrary degree 1d spline
  type :: sll_t_spline_pp_1d

     sll_int32 :: degree !< degree of 1d spline
     sll_real64, allocatable :: poly_coeffs(:,:) !< poly_coeffs[i,j] coefficient of x^{deg+1-j} for ith B-spline function  size= (degree+1, degree+1)
     sll_real64, allocatable :: poly_coeffs_fp(:,:) !< poly_coeffs[i,j] coefficient of x^{deg+1-j} for ith B-spline function  size= (degree+1, degree+1)
     sll_int32 :: n_cells !< number of gridcells
     sll_real64, allocatable :: scratch_b(:) !< scratch data for b_to_pp-converting
     sll_real64, allocatable :: scratch_p(:) !< scratch data for b_to_pp-converting
     
  end type sll_t_spline_pp_1d
  !>arbitrary degree 2d spline
  type :: sll_t_spline_pp_2d
     type(sll_t_spline_pp_1d) spline1 !< first 1d spline
     type(sll_t_spline_pp_1d) spline2 !< second 1d spline

  end type sll_t_spline_pp_2d
  
  !>arbitrary degree 3d spline
  type :: sll_t_spline_pp_3d
     type(sll_t_spline_pp_1d) spline1 !< first 1d spline
     type(sll_t_spline_pp_1d) spline2 !< second 1d spline
     type(sll_t_spline_pp_1d) spline3 !< third 1d spline
     
     sll_real64, allocatable :: scratch_b(:)  !< scratch data for b_to_pp-converting
     sll_real64, allocatable :: scratch_p(:)  !< scratch data for b_to_pp-converting
     sll_real64, allocatable :: scratch_pp(:) !< scratch data for b_to_pp-converting
     
  end type sll_t_spline_pp_3d
  
contains
  
  !> Convert 1d spline in B form to spline in pp form with periodic boundary conditions
  subroutine sll_s_spline_pp_b_to_pp_1d( self, n_cells, b_coeffs, pp_coeffs)
    type( sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline 
    sll_int32, intent(in) ::  n_cells !< number of gridcells
    sll_real64,intent(in) :: b_coeffs(n_cells)   !< coefficients of spline in B-form
    sll_real64,intent(out):: pp_coeffs(self%degree+1,n_cells)  !< coefficients of spline in pp-form
    sll_int32 :: i
       
    do i=1, self%degree   
       call sll_s_spline_pp_b_to_pp_1d_cell(self,(/b_coeffs(n_cells-self%degree+i:n_cells),b_coeffs(1:i)/),pp_coeffs(:,i))
    end do
    
    do i=self%degree+1,n_cells
       call sll_s_spline_pp_b_to_pp_1d_cell(self,b_coeffs(i-self%degree:i),pp_coeffs(:,i))
    end do

  end subroutine sll_s_spline_pp_b_to_pp_1d

  !> Convert 1d spline in B form in a cell to spline in pp form with periodic boundary conditions
  subroutine sll_s_spline_pp_b_to_pp_1d_cell( self, b_coeffs, pp_coeffs )
    type( sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline 
    sll_real64, intent( in )  :: b_coeffs(self%degree+1)   !< coefficients of spline in B-form
    sll_real64, intent( out ) :: pp_coeffs(self%degree+1)  !< coefficients of spline in pp-form
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: degp1
    degp1 = self%degree+1
    pp_coeffs=0.0_f64 
    do i=1, degp1
       do j=1, degp1
          pp_coeffs(j)=pp_coeffs(j)+b_coeffs(i)* self%poly_coeffs(j,i)
       end do
    end do
    
  end subroutine sll_s_spline_pp_b_to_pp_1d_cell
 

  !> Convert 2d spline in B form to spline in pp form with periodic boundary conditions   
  subroutine sll_s_spline_pp_b_to_pp_2d( self, n_cells, b_coeffs, pp_coeffs)
    type( sll_t_spline_pp_2d), intent(inout)::  self !< arbitrary degree 2d spline 
    sll_int32,  intent(in) ::  n_cells(2) !< number of gridcells
    sll_real64, intent(in) :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
    sll_real64, intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1),n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
    sll_int32 :: i,j
    sll_int32 :: degree1,degree2
    degree1= self%spline1%degree
    degree2= self%spline2%degree
     
    do j=1, n_cells(2)
       do i=1, n_cells(1)
          call sll_s_spline_pp_b_to_pp_2d_cell(self%spline1, self%spline2, n_cells, b_coeffs, pp_coeffs, i,j)
       end do
    end do
    
  end subroutine sll_s_spline_pp_b_to_pp_2d
  
  !> Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions 
  subroutine sll_s_spline_pp_b_to_pp_2d_cell(spline1,spline2,n_cells, b_coeffs, pp_coeffs,i,j)
    type( sll_t_spline_pp_1d), intent(inout)::  spline1 !< arbitrary degree 1d spline
    type( sll_t_spline_pp_1d), intent(inout)::  spline2 !< arbitrary degree 1d spline
    sll_int32, intent(in)    :: n_cells(2) !< number of gridcells
    sll_real64,intent(in)    :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((spline1%degree+1)*(spline2%degree+1),n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
    sll_int32, intent(in)    :: i,j !< indices 
    sll_int32 :: k,l
    sll_int32 :: degree1,degree2,degp1,degp2
    degree1= spline1%degree
    degree2= spline2%degree
    degp1=degree1+1
    degp2=degree2+1
    !> convert b-coefficients in pp-coefficients in first dimension
    if (i>degree1) then
       if(j>degree2) then
          do l=0,degree2
             spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+l)*n_cells(1):i+(j-degp2+l)*n_cells(1))
             call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
          end do
       else 
          !> use of modulo for boundary cells in second dimension 
          do l=0,degree2
             spline1%scratch_b=b_coeffs(i-degree1+modulo(j-degp2+l,n_cells(2))*n_cells(1):i+modulo(j-degp2+l,n_cells(2))*n_cells(1))
             call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
          end do
       end if
    else 
       !> use of modulo for boundary cells in both dimensions 
       do l=0,degree2
          do k=0, degree1
             spline1%scratch_b(k+1)=b_coeffs(modulo(i-degp1+k,n_cells(1))+1+modulo(j-degp2+l,n_cells(2))*n_cells(1))
          end do

          call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
       end do
    end if
    
    !> convert b-coefficients in pp_coefficients in second dimension
    do k=1, degp1
       do l=1,degp2
        spline2%scratch_b(l)=pp_coeffs(k+degp1*(l-1),i+n_cells(1)*(j-1))
       end do
       call sll_s_spline_pp_b_to_pp_1d_cell(spline2, spline2%scratch_b, spline2%scratch_p)
       do l=1,degp2
          pp_coeffs(k+degp1*(l-1),i+n_cells(1)*(j-1))=spline2%scratch_p(l)
       end do
    end do
  end subroutine sll_s_spline_pp_b_to_pp_2d_cell
  

  !> Convert 3d spline in B form to spline in pp form with periodic boundary conditions   
  subroutine sll_s_spline_pp_b_to_pp_3d( self, n_cells, b_coeffs, pp_coeffs)
    type( sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 2d spline 
    sll_int32,  intent(in) ::  n_cells(3) !< number of gridcells
    sll_real64, intent(in) :: b_coeffs(n_cells(1)*n_cells(2),n_cells(2))   !< coefficients of spline in B-form
    sll_real64, intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
    sll_int32 :: i,j,k
    sll_int32 :: degree1,degree2,degree3
    degree1= self%spline1%degree
    degree2= self%spline2%degree
    degree3= self%spline3%degree
     
    do k=1, n_cells(3)
       do j=1, n_cells(2)
          do i=1, n_cells(1)
             call sll_s_spline_pp_b_to_pp_3d_cell(self, n_cells, b_coeffs, pp_coeffs, i,j,k)
          end do
       end do
    end do
    
  end subroutine sll_s_spline_pp_b_to_pp_3d

  !> Convert 3d spline in B form in a cell to spline in pp form with periodic boundary conditions 
  subroutine sll_s_spline_pp_b_to_pp_3d_cell(self,n_cells, b_coeffs, pp_coeffs,i,j,k)
    type( sll_t_spline_pp_3d), intent(inout)::  self !< arbitrary degree 3d spline 
    sll_int32, intent(in)    :: n_cells(3) !< number of gridcells
    sll_real64,intent(in)    :: b_coeffs(n_cells(1)*n_cells(2)*n_cells(3))   !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1)*(self%spline3%degree+1),n_cells(1)*n_cells(2)*n_cells(3))  !< coefficients of spline in pp-form
    sll_int32, intent(in)    :: i,j,k !< indices 
    sll_int32 :: l,m,n
    sll_int32 :: degree1,degree2,degree3,degp1,degp2,degp3
    degree1= self%spline1%degree
    degree2= self%spline2%degree
    degree3= self%spline3%degree
    degp1=degree1+1
    degp2=degree2+1
    degp3=degree3+1
    !> convert b-coefficients in pp-coefficients in first dimension
    if (i>degree1) then
       if(j>degree2) then
          if(k>degree3) then
             do m=0, degree2
                do n=0, degree3
                   self%spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+m)*n_cells(1)+(k-degp3+n)*n_cells(1)*n_cells(2):i+(j-degp2+m)*n_cells(1)+(k-degp3+n)*n_cells(1)*n_cells(2))
                   call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2))) 
                end do
             end do
          else
             !> use of modulo for boundary cells in third dimension 
             do m=0, degree2
                do n=0, degree3
                   self%spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+m)*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2):i+(j-degp2+m)*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2))
                   call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
                end do
             end do
          end if
       else
          !> use of modulo for boundary cells in second and third dimension 
          do m=0, degree2
             do n=0, degree3
                self%spline1%scratch_b=b_coeffs(i-degree1+modulo(j-degp2+m,n_cells(2))*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2):i+modulo(j-degp2+m,n_cells(2))*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2))
                call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
             end do
          end do
       end if
    else 
       !> use of modulo for boundary cells in all three dimensions 
       do m=0, degree2
          do n=0, degree3
             do l=0,degree1
                self%spline1%scratch_b(l+1)=b_coeffs(modulo(i-degp1+l,n_cells(1))+1+modulo(j-degp2+m,n_cells(2))*n_cells(1)+modulo(k-degp3+n,n_cells(3))*n_cells(1)*n_cells(2))
             end do
             call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, self%spline1%scratch_b, self%scratch_pp(1+m*degp1+n*degp1*degp2:degp1*(1+m+n*degp2)))
          end do
       end do
    end if
   
    !> convert b-coefficients in pp-coefficients in second dimension  
    do l=1, degp1
       do n=1, degp3
          do m=1, degp2
             self%spline2%scratch_b(m)=self%scratch_pp(l+degp1*(m-1)+degp1*degp2*(n-1))
          end do
          call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, self%spline2%scratch_b, self%spline2%scratch_p)
          do m=1, degp2
             self%scratch_pp(l+degp1*(m-1)+degp1*degp2*(n-1))=self%spline2%scratch_p(m)
          end do
       end do
    end do

    !> convert b-coefficients in pp-coefficients in first dimension
    do l=1, degp1
       do m=1, degp2
          do n=1, degp3
             self%spline3%scratch_b(n)=self%scratch_pp(l+degp1*(m-1)+degp1*degp2*(n-1))
          end do
          call sll_s_spline_pp_b_to_pp_1d_cell(self%spline3, self%spline3%scratch_b, self%spline3%scratch_p)
          do n=1, degp3
             pp_coeffs(l+degp1*(m-1)+degp1*degp2*(n-1),i+n_cells(1)*(j-1)+n_cells(1)*n_cells(2)*(k-1))=self%spline3%scratch_p(n)
          end do
       end do
    end do
 
  end subroutine sll_s_spline_pp_b_to_pp_3d_cell
  
  
  !> Perform a 1d hornerschema on the poly_coeffs
  subroutine sll_s_spline_pp_horner_m_1d(self, val, degree, x)
    type( sll_t_spline_pp_1d), intent(in)::  self !< arbitrary degree 1d spline 
    sll_real64, intent(out):: val(:) !< array of values
    sll_int32, intent(in)  :: degree !< degree of the spline
    sll_real64, intent(in) :: x !< point at which we evaluate our spline
    sll_int32 :: i
    
    do i=1, size(val)
       val(i)=sll_f_spline_pp_horner_1d(degree, self%poly_coeffs, x, i)
    end do
  end subroutine sll_s_spline_pp_horner_m_1d

  !> Perform two times a 1d hornerschema on the poly_coeffs
  subroutine sll_s_spline_pp_horner_m_2d(self, val, degree, x)
    type( sll_t_spline_pp_2d), intent(in)::  self !< arbitrary degree 2d spline 
    sll_real64, intent(out):: val(:,:) !< array of values
    sll_int32, intent(in)  :: degree(2) !< degree of the spline
    sll_real64, intent(in) :: x(2) !< point at which we evaluate our spline
    sll_int32 :: i
      
    do i=1, size(val,1)
       val(i,1)=sll_f_spline_pp_horner_1d(degree(1), self%spline1%poly_coeffs, x(1), i)
       val(i,2)=sll_f_spline_pp_horner_1d(degree(2), self%spline2%poly_coeffs, x(2), i)
    end do
  end subroutine sll_s_spline_pp_horner_m_2d

    !> Perform three times a 1d hornerschema on the poly_coeffs
  subroutine sll_s_spline_pp_horner_m_3d(self, val, degree, x)
    type( sll_t_spline_pp_3d), intent(in)::  self !< arbitrary degree 3d spline 
    sll_real64, intent(out):: val(:,:) !< array of values
    sll_int32, intent(in)  :: degree(3) !< degree of the spline
    sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
    sll_int32 :: i
      
    do i=1, size(val,1)
       val(i,1)=sll_f_spline_pp_horner_1d(degree(1), self%spline1%poly_coeffs, x(1), i)
       val(i,2)=sll_f_spline_pp_horner_1d(degree(2), self%spline2%poly_coeffs, x(2), i)
       val(i,3)=sll_f_spline_pp_horner_1d(degree(3), self%spline3%poly_coeffs, x(3), i)
    end do
  end subroutine sll_s_spline_pp_horner_m_3d

  
  !> Perform a 1d hornerschema on the pp_coeffs evaluate at x
  subroutine sll_s_spline_pp_horner_primitive_1d(val, degree, pp_coeffs, x)
    sll_real64, intent(out):: val(:) !< array of values
    sll_int32, intent(in)  :: degree !< degree of the spline
    sll_real64, intent(in) :: pp_coeffs(:,:)  !< coefficients of spline in pp-form
    sll_real64, intent(in) :: x !< point at which we evaluate our spline
    sll_int32 :: i
    
    do i=1, size(val)
       val(i)=sll_f_spline_pp_horner_1d(degree, pp_coeffs, x, i)*x
    end do
  end subroutine sll_s_spline_pp_horner_primitive_1d
  
  
  !> Perform a 1d hornerschema on the pp_coeffs at index
  function sll_f_spline_pp_horner_1d(degree, pp_coeffs, x, index) result(res)
    sll_int32,  intent(in) :: degree !< degree of the spline
    sll_real64, intent(in) :: pp_coeffs(:,:)  !< coefficients of spline in pp-form
    sll_real64, intent(in) :: x !< point at which we evaluate our spline
    sll_int32,  intent(in) :: index !< index of cell in which is x
    sll_real64 :: res !< value of the splinefunction at point x
    sll_int32 :: i
    
    res=pp_coeffs(1,index)
    do i=1,degree
       res=res*x+pp_coeffs(i+1,index)
    end do
  end function sll_f_spline_pp_horner_1d
  
  !> Perform a 2d hornerschema on the pp_coeffs at the indices
  function sll_f_spline_pp_horner_2d(degree, pp_coeffs, x, indices, n_cells) result(res)
    sll_int32,  intent(in) :: degree(2) !< degree of the spline
    sll_real64, intent(in) :: pp_coeffs(:,:)  !< coefficients of spline in pp-form
    sll_real64, intent(in) :: x(2) !< point at which we evaluate our spline
    sll_int32,  intent(in) :: indices(2) !< indices of cell in which is x
    sll_int32,  intent(in) :: n_cells(2) !< number of gridcells
    sll_real64 :: res !< value of the splinefunction at point x
    sll_real64 :: pp_coeffs_1d(degree(2)+1,1)
    sll_int32  :: i
    !> Perform a 1d hornerschema in the first dimension
    do i=0,degree(2)
       pp_coeffs_1d(i+1,1)=sll_f_spline_pp_horner_1d(degree(1), pp_coeffs(1+i*(degree(1)+1):(degree(1)+1)*(i+1),1+(indices(2)-1)*n_cells(1):n_cells(1)*indices(2)), x(1),indices(1))
    end do
    !> Perform a 1d hornerschema in the second dimension
    res=sll_f_spline_pp_horner_1d(degree(2), pp_coeffs_1d, x(2),1)
  end function sll_f_spline_pp_horner_2d
  
  !> Perform a 3d hornerschema on the pp_coeffs at the indices
  function sll_f_spline_pp_horner_3d (degree, pp_coeffs, x, indices, n_cells) result(res)
    sll_int32,  intent(in) :: degree(3) !< degree of the spline
    sll_real64, intent(in) :: pp_coeffs(:,:)  !< coefficients of spline in pp-form
    sll_real64, intent(in) :: x(3) !< point at which we evaluate our spline
    sll_int32,  intent(in) :: indices(3) !< indices of cell in which is x
    sll_int32,  intent(in) :: n_cells(3) !< number of gridcells
    sll_real64 :: res !< value of the splinefunction at point x
    sll_real64 :: pp_coeffs_2d((degree(2)+1)*(degree(3)+1),1)
    sll_int32  :: i,j
    sll_int32  :: degp1,degp2
    degp1=degree(1)+1
    degp2=degree(2)+1
    !> Perform a 1d hornerschema in the first dimension
    do j=0, degree(3)
       do i=0, degree(2)
          pp_coeffs_2d(1+i+j*degp2,1)= sll_f_spline_pp_horner_1d(degree(1),pp_coeffs(1+i*degp1+j*degp1*degp2:degp1+i*degp1+j*degp1*degp2,1+n_cells(1)*(indices(2)-1+(indices(3)-1)*n_cells(2)):n_cells(1)*(indices(2)+(indices(3)-1)*n_cells(2))),x(1),indices(1))
       end do
    end do
    !> Perform a 2d hornerschema in the second and third dimension
    res=sll_f_spline_pp_horner_2d(degree(2:3),pp_coeffs_2d,x(2:3), [1,1],[1,1])

  end function sll_f_spline_pp_horner_3d
  

  !> Initialize \a sll_t_spline_pp_1d object (Set \a poly_coeffs depending on \a  degree)
  subroutine sll_s_spline_pp_init_1d( self, degree, n_cells)
    type( sll_t_spline_pp_1d), intent(out) ::  self  !< arbitrary degree 1d spline
    sll_int32, intent(in) :: degree !<degree of spline
    sll_int32, intent(in) :: n_cells !< number of gridcells 
    sll_int32 :: ierr  
   
    SLL_ASSERT(n_cells >= degree )
    self%degree = degree
    self%n_cells=n_cells
    SLL_ALLOCATE(self%poly_coeffs(degree+1,degree+1), ierr)
    SLL_ASSERT( ierr == 0)
    SLL_ALLOCATE(self%poly_coeffs_fp(degree+1,degree+1), ierr)
    SLL_ASSERT( ierr == 0)
    SLL_ALLOCATE(self%scratch_b(degree+1), ierr)
    SLL_ASSERT( ierr == 0)
    SLL_ALLOCATE(self%scratch_p(degree+1), ierr)
    SLL_ASSERT( ierr == 0)
    
    select case( self%degree )
    case (1) 
       self%poly_coeffs = reshape( (/-1._f64, 1._f64, 1._f64, 0._f64/), (/2, 2/) ) 
       self%poly_coeffs_fp = reshape( (/-inv_2, 1._f64, inv_2, 0._f64/), (/2,2 /) )
       
    case (2) 
       self%poly_coeffs = reshape((/inv_2, -1._f64, inv_2&
            , -1._f64, 1._f64, inv_2&
            , inv_2, 0._f64, 0._f64/), (/3,3/))
       
        self%poly_coeffs_fp = reshape((/inv_6, -inv_2, inv_2&
             , -inv_3, inv_2, inv_2&
             , inv_6, 0._f64, 0._f64/), (/3,3/))

    case(3) 
       self%poly_coeffs = reshape((/-inv_6, inv_2, -inv_2 , inv_6&
            , inv_2, -1._f64, 0._f64, 4._f64*inv_6&
            , -inv_2, inv_2, inv_2, inv_6&
            , inv_6, 0._f64, 0._f64, 0._f64/), (/4, 4/))
       
       self%poly_coeffs_fp = reshape((/-inv_24, inv_6, -inv_4 , inv_6&
            , inv_8, -inv_3, 0._f64, 4._f64*inv_6&
            , -inv_8, inv_6, inv_4, inv_6&
            , inv_24, 0._f64, 0._f64, 0._f64/), (/4, 4/) )

    case(4) 
       self%poly_coeffs= reshape((/inv_24,-inv_6,inv_4,-inv_6,inv_24&
            ,-inv_6,inv_2,-inv_4,-inv_2,11._f64*inv_24&
            ,inv_4,-inv_2,-inv_4,inv_2,11._f64*inv_24&
            ,-inv_6,inv_6,inv_4,inv_6,inv_24&
            ,inv_24,0._f64,0._f64,0._f64,0._f64/), (/5,5/)) 

       self%poly_coeffs_fp = reshape((/inv_120,-inv_24,inv_12,-inv_12,inv_24&
            ,-inv_30,inv_8,-inv_12,-inv_2,11._f64*inv_24&
            ,inv_20,-inv_8,-inv_12,inv_4,11._f64*inv_24&
            , -inv_30,inv_24,inv_12,inv_12,inv_24&
            ,inv_120,0._f64,0._f64,0._f64,0._f64/), (/5,5/))

    case(5) 
       self%poly_coeffs =reshape((/-inv_120,inv_24,-inv_12,inv_12,-inv_24,inv_120&
            ,inv_24,-inv_6,inv_6,inv_6,-5._f64*inv_12, 26._f64*inv_120& 
            ,-inv_12,inv_4,0._f64,-inv_2,0._f64,11._f64*inv_20&
            ,inv_12,-inv_6,-inv_6,inv_6,5._f64*inv_12,26._f64*inv_120& 
            ,-inv_24,inv_24,inv_12,inv_12,inv_24,inv_120&
            ,inv_120,0._f64,0._f64,0._f64,0._f64,0._f64/),(/6,6/))

     self%poly_coeffs_fp = reshape((/-inv_720,inv_120,-inv_48,inv_36,-inv_48,inv_120&
          ,inv_144,-inv_30,inv_24,inv_18,-5._f64*inv_24, 26._f64*inv_120& 
          ,-inv_72,inv_20,0._f64,-inv_6,0._f64,11._f64*inv_20&
          ,inv_72,-inv_30,-inv_24,inv_18,5._f64*inv_24,26._f64*inv_120& 
          ,-inv_144,inv_120,inv_48,inv_36,inv_48,inv_120&
          ,inv_720,0._f64,0._f64,0._f64,0._f64,0._f64/),(/6,6/)) 


    case default
       SLL_ERROR('sll_s_spline_pp_init','Not implemented')
    end select
  end subroutine sll_s_spline_pp_init_1d

  !> Initialize \a sll_t_spline_pp_2d object
  subroutine sll_s_spline_pp_init_2d(self,degree,n_cells)
    type(sll_t_spline_pp_2d) self !< arbitrary degree 2d spline
    sll_int32, intent(in) :: degree(2) !< degrees of the 1d splines
    sll_int32, intent(in) :: n_cells(2) !< number of gridcells in every dimension
           
    call sll_s_spline_pp_init_1d(self%spline1, degree(1),n_cells(1))
    call sll_s_spline_pp_init_1d(self%spline2, degree(2),n_cells(2))

  end subroutine sll_s_spline_pp_init_2d

  !> Initialize \a sll_t_spline_pp_3d object
  subroutine sll_s_spline_pp_init_3d(self,degree,n_cells)
    type(sll_t_spline_pp_3d) self !< arbitrary degree 3d spline
    sll_int32, intent(in) :: degree(3) !< degrees of the 1d splines
    sll_int32, intent(in) :: n_cells(3) !< number of gridcells in every dimension
    sll_int32 :: ierr  
        
    call sll_s_spline_pp_init_1d(self%spline1, degree(1),n_cells(1))
    call sll_s_spline_pp_init_1d(self%spline2, degree(2),n_cells(2))
    call sll_s_spline_pp_init_1d(self%spline3, degree(3),n_cells(3))
    SLL_ALLOCATE(self%scratch_b((degree(2)+1)*(degree(3)+1)), ierr)
    SLL_ASSERT( ierr == 0)
    SLL_ALLOCATE(self%scratch_p((degree(2)+1)*(degree(3)+1)), ierr)
    SLL_ASSERT( ierr == 0)
    SLL_ALLOCATE(self%scratch_pp((degree(1)+1)*(degree(2)+1)*(degree(3)+1)), ierr)
    SLL_ASSERT( ierr == 0)

  end subroutine sll_s_spline_pp_init_3d
  

  !> Destructor 1d  
  subroutine sll_s_spline_pp_free_1d(self)
    type( sll_t_spline_pp_1d),intent(inout) ::  self !< arbitrary degree 1d spline
   
    deallocate(self%poly_coeffs)
    deallocate(self%poly_coeffs_fp)
    deallocate(self%scratch_b)
    deallocate(self%scratch_p)
    
  end subroutine sll_s_spline_pp_free_1d

  !> Destructor 2d  
  subroutine sll_s_spline_pp_free_2d(self)
    type( sll_t_spline_pp_2d),intent(inout) ::  self !< arbitrary degree 2d spline

    call sll_s_spline_pp_free_1d(self%spline1)
    call sll_s_spline_pp_free_1d(self%spline2)
    
  end subroutine sll_s_spline_pp_free_2d

  !> Destructor 3d  
  subroutine sll_s_spline_pp_free_3d(self)
    type( sll_t_spline_pp_3d),intent(inout) ::  self !< arbitrary degree 3d spline

    call sll_s_spline_pp_free_1d(self%spline1)
    call sll_s_spline_pp_free_1d(self%spline2)
    call sll_s_spline_pp_free_1d(self%spline3)
    deallocate(self%scratch_b)
    deallocate(self%scratch_p)
    deallocate(self%scratch_pp)

  end subroutine sll_s_spline_pp_free_3d
  
end module sll_m_splines_pp
