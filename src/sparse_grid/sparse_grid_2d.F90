module sparse_grid_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_module_periodic_interpolator_1d
use sll_arbitrary_degree_splines
use sll_module_lagrange_interpolator_1d
use sll_sparse_grid_interpolator
use, intrinsic :: iso_c_binding


implicit none

  type fft_fg_2d
     type(C_PTR)  :: bw
     complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: in
     real(C_DOUBLE), dimension(:,:), pointer :: out
     type(C_PTR) :: p_in, p_out
     integer(C_SIZE_T) :: sz
  end type fft_fg_2d

#ifdef STDF95
type :: sparse_grid_interpolator_2d
#else
type, extends(sparse_grid_interpolator) :: sparse_grid_interpolator_2d
#endif
sll_int32, dimension(:,:), pointer  :: index
type(fft_fg_2d) :: fft_object_fg

#ifdef STDF95
#else
contains
  procedure, pass(interpolator) :: initialize => initialize_sg2d! Initialization routine
  procedure, pass(interpolator) :: interpolate_const_disp
 procedure, pass(interpolator) :: interpolate_value=>interpolate_value_sg ! Compute the value of the sparse grid interpolant at position eta
 procedure, pass(interpolator) :: interpolate_disp_nconst_in_1d ! Interpolate along one (x)-direction with displacement non-constant in one (v)-direction
 procedure, pass(interpolator) :: interpolate_disp_linnconst_in_1d ! Interpolate along one (x)-direction with displacement non-constant in one (v)-direction
 procedure, pass(interpolator) :: fg_to_sg
 procedure, pass(interpolator) :: SPFFT ! Compute the Sparse grid FFT coefficients
 procedure, pass(interpolator) :: interpolate_array_disp_sgfft ! Compute value at displaced grid points using trigonometric interpolation (based on SG FFT)
 procedure, pass(interpolator) :: sg_to_fg_sgfft
 procedure, pass(interpolator) :: filter_highest
 procedure, pass(interpolator) :: filter
 procedure, pass(interpolator) :: linear_filter
#endif
end type sparse_grid_interpolator_2d

contains

!------------------------------------------------------------------------------!
!!!! SGFFT routines !!!!!

subroutine interpolate_array_disp_sgfft(interpolator,dim, displacment_in,data_in,data_out)
  class(sparse_grid_interpolator_2d),  intent(inout)       :: interpolator
  sll_int32, intent(in) :: dim
  sll_real64, intent(in) :: displacment_in
  sll_comp64, dimension(:), intent(inout)   :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_real64:: displacement
  
  displacement = displacment_in*2.0_f64*sll_pi/interpolator%length(dim)
  call Displace(interpolator,dim,displacement,data_in);
  call ISPFFT(interpolator,data_in,data_out);
  
end subroutine interpolate_array_disp_sgfft


subroutine sg_to_fg_sgfft(interpolator,data_in,fourier_sg,data_out)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: fourier_sg
  sll_real64, dimension(:,:), intent(out) :: data_out
  sll_int32 :: i1,i2

  call SPFFT(interpolator,data_in,fourier_sg);
  interpolator%fft_object_fg%in = 0.0_f64
  call sg_to_fg_complex(interpolator, fourier_sg,&
       interpolator%fft_object_fg%in);
  do i1=1,32
     do i2 = 1,32
        write(35,*) real(interpolator%fft_object_fg%in(i1,i2)), aimag(interpolator%fft_object_fg%in(i1,i2))
     end do
  end do
 ! print*, sum(fourier_sg)
!#ifdef FFTW_F2003
!  call fftw_execute_dft_c2r(interpolator%fft_object_fg%bw, &
!       interpolator%fft_object_fg%in, data_out);
!  call rfftwnd_one_complex_to_real(interpolator%plan_fg,fourier_fg,data_out);
!#endif  

end subroutine sg_to_fg_sgfft

!!!! End SGFFT routines !!!!!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!!!!!!!!!! Initialization routines !!!!!!!!

  subroutine initialize_sg2d( &
    interpolator, &
    levels, &
    order, &
    interpolation, &
    interpolation_type, &
    eta_min, &
    eta_max, &
    boundary, &
    modified)

#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in)          :: eta_min
    sll_real64, dimension(:),  intent(in)         :: eta_max
    sll_int32, dimension(:),intent(in)            :: levels
    sll_int32, intent(in)                         :: order, interpolation
    sll_int32, intent(in)                         :: interpolation_type
    sll_int32, intent(in)                         :: modified, boundary
    sll_int32                                     :: i,j,k1,k2,l1,l2,l,counter
    sll_int32                                     :: ierr
    sll_int32, dimension(:) ,allocatable          :: novec,lvec,kvec
    sll_int32, dimension(2) :: sizefg
    sll_int32 :: no1, no2

    interpolator%dim = 2;
    interpolator%modified = modified;
    interpolator%boundary = boundary;
    interpolator%max_level = levels(1);
    do l=2,interpolator%dim
       interpolator%max_level = max(levels(l),interpolator%max_level);
    end do
    if (interpolator%modified == 1) then
       interpolator%max_level = interpolator%max_level + 1;
    end if

    interpolator%size_basis = 0;
    if (interpolator%boundary == 0) then
       do l = 0, interpolator%max_level
          do l1 = max(0, l-levels(2)) , min(l, levels(1))
             l2 = l-l1
             interpolator%size_basis = &
                  interpolator%size_basis + &
                  max(2**(l1-1),1)*max(2**(l2-1),1);
          end do
       end do
    else
       do l = 0, interpolator%max_level
          do l1 = max(0, l-levels(2)) , min(l, levels(1))
             l2 = l-l1
             if( l1 == 0) then
                no1 = 2;
             else
                no1 = 2**(l1-1);
             end if
             if( l2 == 0) then
                no2 = 2;
             else
                no2 = 2**(l2-1);
             end if
             interpolator%size_basis = &
                  interpolator%size_basis + no1*no2;
          end do
       end do
    end if

    call  initialize_sg( interpolator, levels, order, interpolation,& 
    interpolation_type, eta_min, eta_max);

    SLL_ALLOCATE(lvec(interpolator%dim),ierr);
    SLL_ALLOCATE(kvec(interpolator%dim),ierr);
    SLL_ALLOCATE(novec(interpolator%dim),ierr);

    SLL_ALLOCATE(interpolator%index(0:interpolator%levels(1),0:interpolator%levels(2)),ierr);
    ! Set the hierarchy of the grid
    counter = 1
    do l = 0, interpolator%max_level
       interpolator%level_mapping(l) = counter;
       do l1 = max(0, l-levels(2)) , min(l, levels(1))
          novec(1) = max(2**(l1-1),1)
          if (interpolator%boundary == 1 .AND. l1 == 0) then
             novec(1) = 2;
          end if
          lvec(1) = l1;
          l2 = l-l1
          novec(2) = max(2**(l2-1),1)
          if (interpolator%boundary == 1 .AND. l2 == 0) then
             novec(2) = 2;
          end if
          lvec(2) = l2;
          interpolator%index(l1,l2) = counter
          do k1=0,novec(1)-1
             kvec(1) = k1
             do k2=0,novec(2)-1
                kvec(2) = k2
                interpolator%hierarchy(counter)%index_on_level(1) = k1;
                interpolator%hierarchy(counter)%index_on_level(2) = k2;
                do j=1,interpolator%dim
                   if (interpolator%boundary == 0) then
                      call set_hierarchy_info&
                           (interpolator,counter,j,lvec,kvec,novec);
                   elseif (interpolator%boundary == 1) then
                      call set_hierarchy_info_boundary&
                           (interpolator,counter,j,lvec,kvec,novec);
                   end if
                end do
                counter = counter +1;
             end do
          end do
       end do
    end do
    interpolator%level_mapping(interpolator%max_level+1) = counter;

    ! Now rescale all the coordinates to the actual mesh size
    do i=1,interpolator%size_basis
       do j=1,interpolator%dim
          interpolator%hierarchy(i)%coordinate(j) = &
               interpolator%hierarchy(i)%coordinate(j)* &
               interpolator%length(j) + &
               interpolator%eta_min(j)
       end do
    end do

!!$    if (fft_fg_init == 'true') then
!!$       ! Initialize fft_object_fg
!!$       sizefg(1) = 2**interpolator%levels(1);
!!$       sizefg(2) = 2**interpolator%levels(2);
!!$       interpolator%fft_object_fg%sz = int(sizefg(1)*sizefg(2),C_SIZE_T)
!!$   
!!$       interpolator%fft_object_fg%p_in = &
!!$            fftw_alloc_complex(interpolator%fft_object_fg%sz)
!!$       call c_f_pointer(interpolator%fft_object_fg%p_in,&
!!$            interpolator%fft_object_fg%in,[sizefg(1), sizefg(2)])
!!$       interpolator%fft_object_fg%p_out = &
!!$            fftw_alloc_real(interpolator%fft_object_fg%sz) 
!!$       call c_f_pointer(interpolator%fft_object_fg%p_out,&
!!$            interpolator%fft_object_fg%out,[sizefg(1), sizefg(2)])
!!$
!!$       interpolator%fft_object_fg%bw = fftw_plan_dft_c2r_2d&
!!$            (sizefg(1), &
!!$            sizefg(2), interpolator%fft_object_fg%in, &
!!$            interpolator%fft_object_fg%out, FFTW_PATIENT);
!!$       
!!$       !fftw2d_create_plan(interpolator%dim, &
!!$       !      interpolator%fft_object_fg%out,&
!!$       !      FFTW_FORWARD, FFTW_MEASURE);
!!$       print*, 'Size fg: ',interpolator%fft_object_fg%sz
!!$    end if

  end subroutine initialize_sg2d

subroutine set_hierarchy_info(interpolator,counter,cdim,lvecin,kvecin,novecin)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
  sll_int32 :: ld ! current level
  sll_int32 :: kd ! current index within level
  sll_int32,intent(in) :: cdim ! current dimension
  sll_int32,intent(in) :: counter ! counter for node
  sll_int32, dimension(:), intent(in) :: lvecin,kvecin,novecin
  sll_int32, dimension(:), allocatable :: lvec,kvec,novec
  sll_int32 :: jj,stride

  SLL_ALLOCATE(lvec(interpolator%dim),jj);
  SLL_ALLOCATE(kvec(interpolator%dim),jj);
  SLL_ALLOCATE(novec(interpolator%dim),jj);

  do jj=1,interpolator%dim
     lvec(jj) = lvecin(jj);
     kvec(jj) = kvecin(jj);
     novec(jj) = novecin(jj);
  end do
  ld = lvec(cdim);
  kd = kvec(cdim);

  interpolator%hierarchy(counter)%level(cdim) = ld;

  stride = cdim*2-1
  if (ld==0) then
     interpolator%hierarchy(counter)%coordinate(cdim) = 0
     interpolator%hierarchy(counter)%parent(stride) = &
          counter
     interpolator%hierarchy(counter)%parent(stride+1) = &
          counter
     interpolator%hierarchy(counter)%function_type(cdim) = 0
  else
     interpolator%hierarchy(counter)%coordinate(cdim) = &
          1.0_f64/(2.0_f64**ld)+kd*1.0_f64/(2.0_f64**(ld-1))


     lvec(cdim) = lvec(cdim)-1;
     novec(cdim) = max(novec(cdim)/2,1);

     ! This one is actually only a neighboring point not the direct parent.
     kvec(cdim) = modulo((kd+1)/2,max(2**(ld-2),1));
     interpolator%hierarchy(counter)%parent(&
          modulo(kd,2)+stride) = &
          interpolator%hierarchy(&
          interpolator%index(lvec(1),lvec(2))+&
          kvec(1)*novec(2)+kvec(2))%parent(stride)
     ! This is the actual parent.
     kvec(cdim) = kd/2;
     interpolator%hierarchy(counter)%parent(&
          modulo(kd+1,2)+stride) = &
          interpolator%index(lvec(1),lvec(2))+&
          kvec(1)*novec(2)+kvec(2)
     ! Now tell my parent that I am his child.
     interpolator%hierarchy(&
          interpolator%hierarchy(counter)%parent(&
          modulo(kd+1,2)+stride))%children(modulo(kd,2)+stride) = counter
     if(ld==1) then
        interpolator%hierarchy(&
             interpolator%hierarchy(counter)%parent(&
             modulo(kd+1,2)+stride))%children(modulo(kd+1,2)+stride) = counter
     end if
     if (interpolator%order == 1) then
        interpolator%hierarchy(counter)%&
             function_type(cdim) = -1
     elseif (ld==1 .OR. interpolator%order==2) then
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 1+modulo(kd,2)
     elseif (ld==2 .OR. interpolator%order==3)then !(order==3) then!
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 3 + modulo(kd,4)
     else
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 7 + modulo(kd,8)
     end if
  end if

end subroutine set_hierarchy_info


subroutine set_hierarchy_info_boundary(interpolator,counter,cdim,lvecin,kvecin,novecin)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
  sll_int32 :: ld ! current level
  sll_int32 :: kd ! current index within level
  sll_int32,intent(in) :: cdim ! current dimension
  sll_int32,intent(in) :: counter ! counter for node
  sll_int32, dimension(:), intent(in) :: lvecin,kvecin,novecin
  sll_int32, dimension(:), allocatable :: lvec,kvec,novec
  sll_int32 :: jj,stride

  SLL_ALLOCATE(lvec(interpolator%dim),jj);
  SLL_ALLOCATE(kvec(interpolator%dim),jj);
  SLL_ALLOCATE(novec(interpolator%dim),jj);

  do jj=1,interpolator%dim
     lvec(jj) = lvecin(jj);
     kvec(jj) = kvecin(jj);
     novec(jj) = novecin(jj);
  end do
  ld = lvec(cdim);
  kd = kvec(cdim);

  interpolator%hierarchy(counter)%level(cdim) = ld;

  stride = cdim*2-1
  if (ld==0) then
     interpolator%hierarchy(counter)%coordinate(cdim) = real(kd,f64);
     interpolator%hierarchy(counter)%parent(stride) = &
          counter
     interpolator%hierarchy(counter)%parent(stride+1) = &
          counter
     interpolator%hierarchy(counter)%function_type(cdim) = -1
  else
     interpolator%hierarchy(counter)%coordinate(cdim) = &
          1.0_f64/(2.0_f64**ld)+kd*1.0_f64/(2.0_f64**(ld-1))


     lvec(cdim) = lvec(cdim)-1;
     novec(cdim) = novec(cdim)/2;

     if(ld==1) then
        kvec(cdim) = kd/2;
        novec(cdim) = 2;
        interpolator%hierarchy(counter)%parent(stride) = &
             interpolator%index(lvec(1),lvec(2))+&
             kvec(1)*novec(2)+kvec(2)
        if (cdim==1) then
           interpolator%hierarchy(counter)%parent(&
                stride+1) = interpolator%hierarchy(counter)%parent(&
                stride) + novec(2)
        else
           interpolator%hierarchy(counter)%parent(&
             stride+1) = interpolator%hierarchy(counter)%parent(&
             stride) + 1
        end if
        interpolator%hierarchy(interpolator%hierarchy(counter)%parent(stride))%children&
             (stride) = counter
        !interpolator%hierarchy(interpolator%hierarchy(counter)%parent(stride+1))%&
        !     children(stride) = counter
     else
        ! This one is actually only a neighboring point not the direct parent.
        kvec(cdim) =modulo((kd+1)/2,max(2**(ld-2),1));
        kvec(cdim) = (kd+1)/2;
        if (kvec(cdim)<novec(cdim)) then
           interpolator%hierarchy(counter)%parent(&
                modulo(kd,2)+stride) = &
                interpolator%hierarchy(&
                interpolator%index(lvec(1),lvec(2))+&
                kvec(1)*novec(2)+kvec(2))%parent(stride)
        else
           kvec(cdim) = kvec(cdim)-1;
           interpolator%hierarchy(counter)%parent(&
                modulo(kd,2)+stride) = &
                interpolator%hierarchy(&
                interpolator%index(lvec(1),lvec(2))+&
                kvec(1)*novec(2)+kvec(2))%parent(stride+1)
        end if
        ! This is the actual parent.
        kvec(cdim) = kd/2;
        interpolator%hierarchy(counter)%parent(&
             modulo(kd+1,2)+stride) = &
             interpolator%index(lvec(1),lvec(2))+&
             kvec(1)*novec(2)+kvec(2)
        ! Now tell my parent that I am his child.
        interpolator%hierarchy(&
             interpolator%hierarchy(counter)%parent(&
             modulo(kd+1,2)+stride))%children(modulo(kd,2)+stride) = counter
     end if

     if (interpolator%order == 1) then
        interpolator%hierarchy(counter)%&
             function_type(cdim) = -1
     elseif (ld==1 .OR. interpolator%order==2) then
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 1+modulo(kd,2)
     elseif (ld==2 .OR. interpolator%order==3)then !(order==3) then!
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 3 + modulo(kd,4)
     else
        interpolator%hierarchy(counter)%&
             function_type(cdim) = 7 + modulo(kd,8)
     end if
  end if

end subroutine set_hierarchy_info_boundary


!!!! End initialization routines !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!!! Interpolation functions !!!!!!



function interpolate_value_sg( interpolator,data,  eta ) result(val)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in) :: data
    sll_real64 :: val
    sll_real64, dimension(:), intent(in) :: eta
    if (interpolator%boundary == 0) then
       val =  interpolate_from_hierarchical_surplus(interpolator,data,eta)
    else
       val = interpolate_from_hierarchical_surplus_boundary(interpolator,data,eta)
    end if

  end function interpolate_value_sg


! helper function for interpolate_value

 function interpolate_from_hierarchical_surplus( interpolator,data, eta ) result(val)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
    sll_int32 :: j,l1,l2,level
    sll_real64 :: val
    sll_real64, dimension(:), intent(in) :: data
    sll_real64,dimension(:), intent(in) :: eta
    sll_real64,dimension(2) :: eta_norm
    sll_real64,dimension(2) :: phi
    sll_int32, dimension(2) :: no,l
    sll_int32,dimension(:,:), allocatable :: ind
    sll_real64 :: scale
    sll_int32 :: index

    SLL_ALLOCATE(ind(0:interpolator%max_level,1:interpolator%dim),j)

    val = 0
    ind(0:1,1:interpolator%dim) = 0

    do j=1,interpolator%dim
       eta_norm(j) = (eta(j)-interpolator%eta_min(j))/interpolator%length(j)
       eta_norm(j) = modulo(eta_norm(j),1.0_f64)
 
       scale = 0.5_f64
       do level = 2, interpolator%levels(j)
          ind(level,j) = ind(level-1,j)*2
          if (eta_norm(j)> scale*(ind(level,j)+1)) then
             ind(level,j) = ind(level,j)+1
          end if
          scale = scale*0.5_f64
       end do
    end do


    do level = 0, interpolator%max_level
       do l1 = max(0, level-interpolator%levels(2)) , min(level, interpolator%levels(1))
          l(1) = l1
          no(1) = max(2**(l1-1),1)
          l2 = level-l1
          l(2) = l2
          no(2) = max(2**(l2-1),1)
          index = interpolator%index(l1,l2)+ind(l1,1)*no(2)&
                     +ind(l2,2)
          do j=1,interpolator%dim
             call basis_function(real(2**(max(l(j),1)),f64)*eta_norm(j)&
                  -real(2*ind(l(j),j),f64)-1.0_f64, phi(j),&
                  interpolator%hierarchy(index)%function_type(j))
          end do
          val = val + data(index)&
               *phi(1)*phi(2)
       end do
    end do

  end function interpolate_from_hierarchical_surplus



! interpolation from hierarchical surplus non-periodic

 function interpolate_from_hierarchical_surplus_boundary( interpolator,data, eta ) result(val)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
    sll_int32 :: j,l1,l2,level
    sll_real64 :: val
    sll_real64, dimension(:), intent(in) :: data
    sll_real64,dimension(:), intent(in) :: eta
    sll_real64,dimension(2) :: eta_norm
    sll_real64,dimension(2) :: phi
    sll_int32, dimension(2) :: no,l
    sll_int32,dimension(:,:), allocatable :: ind
    sll_real64 :: scale, phi1a, phi2a
    sll_int32 :: index

    SLL_ALLOCATE(ind(0:interpolator%max_level,1:interpolator%dim),j)

    val = 0
    ind(0:1,1:interpolator%dim) = 0

    do j=1,interpolator%dim
       eta_norm(j) = (eta(j)-interpolator%eta_min(j))/interpolator%length(j)
 
       scale = 0.5_f64
       do level = 2, interpolator%levels(j)
          ind(level,j) = ind(level-1,j)*2
          if (eta_norm(j)> scale*(ind(level,j)+1)) then
             ind(level,j) = ind(level,j)+1
          end if
          scale = scale*0.5_f64
       end do
    end do



    do level = 0, interpolator%max_level
       do l1 = max(0, level-interpolator%levels(2)) , min(level, interpolator%levels(1))
          l(1) = l1
          if (l1 == 0) then
             no(1) = 2;
          else
             no(1) = 2**(l1-1);
          end if
          l2 = level-l1
          l(2) = l2
          if (l2==0) then
             no(2) = 2;
          else
             no(2) = max(2**(l2-1),1);
          end if
          if (l1>0) then
             if (l2>0) then
                index = interpolator%index(l1,l2)+ind(l1,1)*no(2)&
                     +ind(l2,2)

                do j=1,interpolator%dim
                   call basis_function(real(2**l(j),f64)*eta_norm(j)&
                        -real(2*ind(l(j),j),f64)-1.0_f64, phi(j),&
                        interpolator%hierarchy(index)%function_type(j))
                end do
                val = val + data(index)&
                     *phi(1)*phi(2)
             else ! l2=0
                index = interpolator%index(l1,l2)+ind(l1,1)*no(2)
                call basis_function(real(2**l(1),f64)*eta_norm(1)&
                     -real(2*ind(l(1),1),f64)-1.0_f64, phi(1),&
                     interpolator%hierarchy(index)%function_type(1))
                call basis_function(eta_norm(2), phi(2), -1)
                val = val + data(index)*phi(1)*phi(2);
                call basis_function(eta_norm(2)-1.0_f64, phi(2), -1)
                val = val + data(index+1)*phi(1)*phi(2);
             end if
          else !l1 = 0
             if (l2>0) then
                index = interpolator%index(l1,l2)+ind(l2,2)
                call basis_function(real(2**l(2),f64)*eta_norm(2)&
                     -real(2*ind(l(2),2),f64)-1.0_f64, phi(2),&
                     interpolator%hierarchy(index)%function_type(2))
                call basis_function(eta_norm(1), phi(1), -1)
                val = val + data(index)*phi(1)*phi(2);
                call basis_function(eta_norm(1)-1.0_f64, phi(1), -1)
                val = val + data(index+no(2))*phi(1)*phi(2);
             else !l1=l2=0
                index =  interpolator%index(l1,l2);
                 call basis_function(eta_norm(1), phi(1), -1)
                 call basis_function(eta_norm(1)-1.0_f64, phi1a, -1)
                 call basis_function(eta_norm(2), phi(2), -1)
                 call basis_function(eta_norm(2)-1.0_f64, phi2a, -1)
                 val = val + data(index)*phi(1)*phi(2) + &
                      data(index+1)*phi(1)*phi2a + data(index+2)*phi1a*phi(2) + &
                      data(index+3)*phi1a*phi2a
             end if
          end if
       end do
    end do

  end function interpolate_from_hierarchical_surplus_boundary



! Interpolation function for interpolation at (constantly) displaced grid points; displacement only in dimension dim. It is another implementation of the base-class function "interpolate_disp". The advantage is that we can not revisit nodes as we do in the recursive dimension-independently-programmed version.
subroutine interpolate_const_disp(interpolator,dorder,displacement,data_in, data_out,hiera)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
  sll_real64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_real64, intent(in) ::displacement
  sll_int32 :: j,counter,i2,k2
  sll_int32, dimension(:), intent(in) :: dorder
  sll_int32, dimension(2) :: ind_order,no
  sll_int32, dimension(:,:), allocatable :: ind
  logical, intent(in) :: hiera

  SLL_ALLOCATE(ind(interpolator%max_level+1,2), i2);


  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i2 = 0,interpolator%levels(dorder(2))
     ind_order(dorder(2)) = i2
     no(dorder(2)) = max(2**(i2-1),1);
     ind(1,dorder(1)) = 0;
     do k2 = 0,no(dorder(2))-1
        ind(ind_order(dorder(2))+1,dorder(2)) = k2;

        counter = interpolator%index(&
             ind_order(1),ind_order(2))+&
             ind(ind_order(1)+1,1)*no(2)&
             +ind(ind_order(2)+1,2)

        ! Evaluate along dorder(1)-stripe
        call interpolate_disp_1d_periodic&
             (interpolator,displacement,dorder(1),&
             min(interpolator%levels(dorder(1)),interpolator%max_level&
             -ind_order(dorder(2))),counter,data_in,data_out,hiera);
        
     end do
  end do

  if (hiera .EQV. .FALSE.) then
     ! Dehierarchization along dimension dorder(1) only
     do j=interpolator%order,2,-1
        call dehierarchical_part_order&
             (interpolator,data_out,&
             interpolator%dim,2,dorder,j)
     end do

     call dehierarchical_part(interpolator,data_out,&
          interpolator%dim,2,dorder)
  end if

end subroutine Interpolate_const_disp



! Functionality: Interpolates the function values for a displacement on in dimension (periodic b.c. i.e. dimension 1 or 2) where the displacement is allowed to be non-constant in one other dimension (Dirichlet b.c. i.e. dimension 3 or 3).
!INPUT:
! interpolator: interpolator handle
! displacement: Vector containing the values of the displacement (in hierarchical order, one dimensional)
! dorder: Ordering of the dimensions. dorder(1) (=1 or 2) gives the dimension where we want to displace, dorder(2) (=3 or 4) gives the dimension of which the displacement is dependent. dorder(3) = 1 or 2 not dorder(1) and dorder(4) = 3 or 4 not dorder(2).
! data_in: hierarchical surplus of the present function
! data_out: value of the displaced functions

subroutine interpolate_disp_nconst_in_1d(interpolator,displacement,dorder,data_in, data_out)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
  sll_real64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, dimension(:), intent(in) :: dorder
  sll_real64,dimension(:), intent(in) ::displacement
  sll_int32 :: i1,i2,k2,counter,j, index_parent,index_parent_old
  sll_int32, dimension(4) :: no,ind_order
  sll_int32, dimension(:,:), allocatable :: ind
  sll_real64 :: coordinate_self, coordinate_ancestor, factor,disp

  SLL_ALLOCATE(ind(interpolator%max_level+1,2), i1);

  ! Dehierarchization along dimension dorder(1) only
  do j=interpolator%order,2,-1
     call dehierarchical_part_order&
          (interpolator%sparse_grid_interpolator,data_in,1,1,dorder,j)
  end do

  call dehierarchical_part(interpolator%sparse_grid_interpolator,data_in,1,1,dorder)

  ! Interpolation in dorder(1)/dorder(2)-plane
  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i2 = 0,interpolator%levels(dorder(2))
     ind_order(dorder(2)) = i2
     no(dorder(2)) = max(2**(i2-1),1);
     ind(1,dorder(1)) = 0;
     do k2 = 0,no(dorder(2))-1
        ind(ind_order(dorder(2))+1,dorder(2)) = k2;
        disp = displacement(2**(i2-1)+k2+1)

        counter = interpolator%index(&
             ind_order(1),ind_order(2))+&
             ind(ind_order(1)+1,1)*no(2)&
             +ind(ind_order(2)+1,2)

        ! Evaluate along dorder(1)-stripe
        call interpolate_disp_1d_periodic_self&
             (interpolator%sparse_grid_interpolator,disp,dorder(1),&
             min(interpolator%levels(dorder(1)),interpolator%max_level&
             -ind_order(dorder(2))),counter,data_in,data_out)
        ! Evaluate hierarchically along dorder(2) dimension (dorder(1)-stripe-wise)
        index_parent = max(&
             interpolator%hierarchy(counter)%parent(2*dorder(2)-1),&
             interpolator%hierarchy(counter)%parent(2*dorder(2)))
        coordinate_self = interpolator%hierarchy(counter)%coordinate(dorder(2))
        index_parent_old = counter;
        do while(index_parent<index_parent_old)
           coordinate_ancestor = interpolator%hierarchy(index_parent)%&
                coordinate(dorder(2))
           call basis_function((coordinate_self-coordinate_ancestor)/&
                interpolator%length(dorder(2))*&
                2**(interpolator%hierarchy(index_parent)%level(dorder(2))), &
                factor,&
                interpolator%hierarchy(index_parent)%function_type(dorder(2)))
           call interpolate_disp_1d_periodic_for_neighbor&
                (interpolator%sparse_grid_interpolator,disp,factor,&
                dorder(1),min(interpolator%levels(dorder(1)),&
                interpolator%max_level-&
                interpolator%hierarchy(index_parent)%level(dorder(2))),&
                min(interpolator%levels(dorder(1)),interpolator%max_level&
                -ind_order(dorder(2))),index_parent,counter,&
                data_in,data_out)
           index_parent_old = index_parent;
           index_parent =  &
                max(interpolator%hierarchy(index_parent)%parent(2*dorder(2)-1),&
                interpolator%hierarchy(index_parent)%parent(2*dorder(2)))
        end do

     end do
  end do
!  call dehira(interpolator%sparse_grid_interpolator,data_out,2,2,dorder)

end subroutine interpolate_disp_nconst_in_1d

! As previous function but displacement dependent on displacement*coordinate(dorder(2))
subroutine interpolate_disp_linnconst_in_1d(interpolator,displacement,dorder,data_in, data_out)
#ifdef STDF95
    type(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#else	
    class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
#endif
  sll_real64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, dimension(:), intent(in) :: dorder
  sll_real64, intent(in) ::displacement
  sll_int32 :: i1,i2,k2,counter,j, index_parent,index_parent_old
  sll_int32, dimension(4) :: no,ind_order
  sll_int32, dimension(:,:), allocatable :: ind
  sll_real64 :: coordinate_self, coordinate_ancestor, factor,disp

  SLL_ALLOCATE(ind(interpolator%max_level+1,2), i1);

  ! Dehierarchization along dimension dorder(1) only
  do j=interpolator%order,2,-1
     call dehierarchical_part_order&
          (interpolator%sparse_grid_interpolator,data_in,1,1,dorder,j)
  end do

  call dehierarchical_part(interpolator%sparse_grid_interpolator,data_in,1,1,dorder)

  ! Interpolation in dorder(1)/dorder(2)-plane
  ind_order(dorder(1)) = 0
  no(dorder(1)) = 1
  do i2 = 0,interpolator%levels(dorder(2))
     ind_order(dorder(2)) = i2
     no(dorder(2)) = max(2**(i2-1),1);
     ind(1,dorder(1)) = 0;
     do k2 = 0,no(dorder(2))-1
        ind(ind_order(dorder(2))+1,dorder(2)) = k2;

        counter = interpolator%index(&
             ind_order(1),ind_order(2))+&
             ind(ind_order(1)+1,1)*no(2)&
             +ind(ind_order(2)+1,2)
       disp = displacement*interpolator%hierarchy(counter)%coordinate(dorder(2));
        ! Evaluate along dorder(1)-stripe
        call interpolate_disp_1d_periodic_self&
             (interpolator%sparse_grid_interpolator,disp,dorder(1),&
             min(interpolator%levels(dorder(1)),interpolator%max_level-&
             ind_order(dorder(2))),counter,data_in,data_out)
        ! Evaluate hierarchically along dorder(2) dimension (dorder(1)-stripe-wise)
        index_parent = max(&
             interpolator%hierarchy(counter)%parent(2*dorder(2)-1),&
             interpolator%hierarchy(counter)%parent(2*dorder(2)))
        coordinate_self = interpolator%hierarchy(counter)%coordinate(dorder(2))
        index_parent_old = counter;
        do while(index_parent<index_parent_old)
           coordinate_ancestor = interpolator%hierarchy(index_parent)%&
                coordinate(dorder(2))
           call basis_function((coordinate_self-coordinate_ancestor)/&
                interpolator%length(dorder(2))*&
                2**(interpolator%hierarchy(index_parent)%level(dorder(2))), &
                factor,&
                interpolator%hierarchy(index_parent)%function_type(dorder(2)))
           call interpolate_disp_1d_periodic_for_neighbor&
                (interpolator%sparse_grid_interpolator,disp,factor,&
                dorder(1),min(interpolator%levels(dorder(1)),&
                interpolator%max_level-&
                interpolator%hierarchy(index_parent)%level(dorder(2))),&
                min(interpolator%levels(dorder(2)),interpolator%max_level&
                -ind_order(dorder(2))),index_parent,counter,&
                data_in,data_out)
           index_parent_old = index_parent;
           index_parent =  &
                max(interpolator%hierarchy(index_parent)%parent(2*dorder(2)-1),&
                interpolator%hierarchy(index_parent)%parent(2*dorder(2)))
        end do

     end do
  end do
!  call dehira(interpolator%sparse_grid_interpolator,data_out,2,2,dorder)

end subroutine interpolate_disp_linnconst_in_1d



!!$subroutine interpolate_disp_nconst_in_1d(interpolator,displacement,dorder,data_in, data_out)
!!$  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
!!$  sll_real64, dimension(:), intent(inout) :: data_in
!!$  sll_real64, dimension(:), intent(out) :: data_out
!!$  sll_int32, dimension(:),intent(in) :: dorder
!!$  sll_int32 :: dim
!!$  sll_real64,dimension(:), intent(in) ::displacement
!!$  sll_int32 :: j
!!$  
!!$  dim = dorder(1);
!!$
!!$  call dehira(interpolator%sparse_grid_interpolator,data_in,1,1,dorder)
!!$
!!$  call interpolate_disp_recursive(interpolator%sparse_grid_interpolator,interpolator%dim,dim,1,displacement(1), data_in, data_out)
!!$
!!$
!!$  call dehira(interpolator%sparse_grid_interpolator,data_out,2,2,dorder)
!!$

!!$ end subroutine Interpolate_disp_nconst_in_1d



!!!!!!!!! End interpolation functions !!!!!!!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! Derivative interpolation (not implemented yet)


!!$  function interpolate_deriv1_sg2d( interpolator, eta1, eta2 ) result(val)
!!$#ifdef STDF95
!!$    type(sparse_grid_interpolator_2d), intent(in) :: interpolator
!!$#else	
!!$    class(sparse_grid_interpolator_2d), intent(in) :: interpolator
!!$#endif
!!$    sll_real64 :: val
!!$    sll_real64, intent(in) :: eta1
!!$    sll_real64, intent(in) :: eta2
!!$    !val =  interpolate_deriv1_from_hierarchical_surplus(&
!!$    !     interpolator%sparsegrid,eta1,eta2)
!!$    print*, 'Not implemented (interpolate_deriv1)'
!!$
!!$  end function
!!$
!!$  function interpolate_deriv2_sg2d( interpolator, eta1, eta2 ) result(val)
!!$#ifdef STDF95
!!$    type(sparse_grid_interpolator_2d), intent(in) :: interpolator
!!$#else	
!!$    class(sparse_grid_interpolator_2d), intent(in) :: interpolator
!!$#endif
!!$    sll_real64 :: val
!!$    sll_real64, intent(in) :: eta1
!!$    sll_real64, intent(in) :: eta2
!!$    !val =  interpolate_deriv2_from_hierarchical_surplus(&
!!$    !     interpolator%sparsegrid,eta1,eta2)
!!$    print*, 'Not implemented (interpolate_deriv2)'
!!$  end function


!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! Functions to evaluate fg on sg and sg on fg

subroutine fg_to_sg(interpolator,fg_values,sg_values)
sll_real64, dimension(:,:), intent(in) :: fg_values
sll_real64, dimension(:), intent(out) :: sg_values
class(sparse_grid_interpolator_2d), intent(in) :: interpolator
sll_int32 :: j
sll_int32, dimension(2) :: fg_ind

do j=1,interpolator%size_basis
   fg_ind = fg_index(interpolator,j);
   sg_values(j) = fg_values(fg_ind(1),fg_ind(2));
end do

end subroutine fg_to_sg

! Complex version
subroutine sg_to_fg_complex(interpolator,sg_values,fg_values)
sll_comp64, dimension(:,:), intent(out) :: fg_values
sll_comp64, dimension(:), intent(in) :: sg_values
class(sparse_grid_interpolator_2d), intent(in) :: interpolator
sll_int32 :: j
sll_int32, dimension(2) :: fg_ind

do j=1,interpolator%size_basis
   fg_ind = fg_index(interpolator,j);
   fg_values(fg_ind(1),fg_ind(2)) = sg_values(j);
end do

end subroutine sg_to_fg_complex


! Compute the index of a sparse grid node on level "level" with index "index_on_level" on full grid with of max_level
function fg_index(interpolator,sg_index)  
sll_int32, intent(in) :: sg_index
sll_int32, dimension(2) :: fg_index
class(sparse_grid_interpolator_2d), intent(in) :: interpolator
sll_int32 :: j

do j=1,interpolator%dim
   if (interpolator%hierarchy(sg_index)%level(j) == 0) then
      fg_index(j) = 1;
   else
      fg_index(j) = 2**(interpolator%levels(j)-&
           interpolator%hierarchy(sg_index)%level(j))*&
           (1+2*interpolator%hierarchy(sg_index)%index_on_level(j)) + 1;
   end if
end do

end function fg_index



! End functions fg_to_sg and sg_to_fg
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!!!! SGFFT helper functions (Some clean-up needed) !!!!


subroutine ToHierarchical(interpolator,data_in, data_out)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out
  sll_int32 :: i,j
  
! Maybe possible with index in index to make dimension independent
  do i = 0,interpolator%levels(2)
     do j = interpolator%index(0,i),&
          interpolator%index(0,i) + max(2**(i-1),1)-1
        call ToHierarchical1D(interpolator,1,&
             min(interpolator%levels(1),interpolator%max_level-i),&
             j,data_in,data_out)
     end do
  end do
 ! print*, 'Data nach FFT in dim 2'
 ! print*, data_out
 ! print*,''
  do i = 0,interpolator%levels(1)
     do j = interpolator%index(i,0),&
          interpolator%index(i,0) + max(2**(i-1),1)-1
        call ToHierarchical1D_comp(interpolator,2,&
             min(interpolator%levels(2),interpolator%max_level-i),j,data_out)
     end do
  end do
!print*,'Data nach FFT in both dim'
!print*, data_out
!print*,''

end subroutine ToHierarchical


 subroutine ToDehi(interpolator,data_array)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_array
  sll_int32 :: i,j

! Maybe possible with index in index to make dimension independent
  do i = 0,interpolator%levels(1)
     do j = interpolator%index(i,0),&
          interpolator%index(i,0) + max(2**(i-1),1)-1
        call ToDehi1D(interpolator,2,&
             min(interpolator%levels(2),interpolator%max_level-i),j,&
             data_array)
     end do
  end do
  do i = 0,interpolator%levels(2)
     do j = interpolator%index(0,i),&
          interpolator%index(0,i) + max(2**(i-1),1)-1
        call ToDehi1D(interpolator,1,&
             min(interpolator%levels(1),interpolator%max_level-i),j,&
             data_array)
     end do
  end do


end subroutine ToDehi


subroutine ToHira(interpolator,data_array)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_array
  sll_int32 :: i,j

  do i = 0,interpolator%levels(2)
     do j = interpolator%index(0,i),&
          interpolator%index(0,i) + max(2**(i-1),1)-1
        call ToHira1D(interpolator,1,&
             min(interpolator%levels(1),interpolator%max_level-i),j,&
             data_array)
     end do
  end do
  do i = 0,interpolator%levels(1)
     do j = interpolator%index(i,0),&
          interpolator%index(i,0) + max(2**(i-1),1)-1
        call ToHira1D(interpolator,2,&
             min(interpolator%levels(2),interpolator%max_level-i),j,&
             data_array)
     end do
  end do

  
end subroutine ToHira

subroutine ToNodal(interpolator,data_in,data_out)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator  
  sll_comp64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32 :: i,j

  do i = 0,interpolator%levels(1)
     do j = interpolator%index(i,0),&
          interpolator%index(i,0) + max(2**(i-1),1)-1
        call ToNodal1D_comp(interpolator,2,&
             min(interpolator%levels(2),interpolator%max_level-i),&
             j,data_in)
        !call ToNodal1D(interpolator,2,interpolator%sparsegrid%levels-i,&
        !     j,data_in,data_out)
     end do
  end do

  do i = 0,interpolator%levels(2)
     do j = interpolator%index(0,i),&
          interpolator%index(0,i) + max(2**(i-1),1)-1
        call ToNodal1D(interpolator,1,&
             min(interpolator%levels(1),interpolator%max_level-i),&
             j,data_in,data_out)
     end do
  end do
 
end subroutine ToNodal


subroutine Displace(interpolator,dim,displacement,data)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_real64,intent(in) :: displacement
  sll_int32, intent(in) :: dim
  sll_int32 :: i,j
  
  if (dim == 1) then
     ! Maybe possible with index in index to make dimension independent
     do i = 0,interpolator%levels(2)
        do j = interpolator%index(0,i),&
             interpolator%index(0,i) + max(2**(i-1),1)-1
           call Displace1D(interpolator,1,interpolator%levels(1)-i,&
                j,displacement,data)
        end do
     end do
  else
     do i = 0,interpolator%levels(1)
        do j = interpolator%index(i,0),&
             interpolator%index(i,0) + max(2**(i-1),1)-1
           call Displace1D(interpolator,2,interpolator%levels(2)-i,&
                j,displacement,data)
        end do
     end do
  end if


end subroutine DISPLACE

subroutine DisplaceVar(interpolator,alpha1,alpha2,data)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_real64,dimension(:) ,intent(in) :: alpha1,alpha2
  sll_int32 :: i,j,counter
  
! Maybe possible with index in index to make dimension independent
  counter = 1
  do i = 0,interpolator%levels(2)
     do j = interpolator%index(0,i),&
          interpolator%index(0,i) + max(2**(i-1),1)-1
        call Displace1D(interpolator,1,interpolator%levels(1)-i,&
             j,alpha1(counter),data)
        counter = counter +1
     end do
  end do
  counter = 1
  do i = 0,interpolator%levels(1)
     do j = interpolator%index(i,0),&
          interpolator%index(i,0) + max(2**(i-1),1)-1
        call Displace1D(interpolator,2,interpolator%levels(2)-i,&
             j,alpha2(counter),data)
        counter = counter +1 
     end do
  end do


end subroutine DISPLACEVAR

subroutine SPFFT(interpolator,data_in,data_out)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out

  call ToHierarchical(interpolator,data_in,data_out)
  call ToDehi(interpolator,data_out)

end subroutine SPFFT

subroutine ISPFFT(interpolator,data_in, data_out)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out

  call ToHira(interpolator,data_in)
  call ToNodal(interpolator,data_in,data_out)

end subroutine ISPFFT



!!!! End SGFFT helper functions !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! Filter functions !!!!

subroutine filter_highest(interpolator, data)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data
  sll_int32 :: lev, l, j, ind

  lev = interpolator%levels(1);
  do l = 0, interpolator%max_level-lev
     ind = interpolator%index(0,l)
     call extract_periodic(interpolator, 1, lev, ind, data, interpolator%stripe)
     call hierarchical_stripe(interpolator, interpolator%stripe, lev);
     do j=2,2**lev,2
        interpolator%stripe(j) = interpolator%stripe(j)/2;
     end do
     call dehierarchical_stripe(interpolator, interpolator%stripe, lev);
     call insert_periodic(interpolator, 1, lev, ind, interpolator%stripe, data);
  end do

  lev = interpolator%levels(2);
  do l = 0, interpolator%max_level-lev
     ind = interpolator%index(l,0)
     call extract_periodic(interpolator, 2, lev, ind, data, interpolator%stripe)
     call hierarchical_stripe(interpolator, interpolator%stripe, lev);
     do j=2,2**lev,2
        interpolator%stripe(j) = interpolator%stripe(j)*0.5_f64;
     end do
     call dehierarchical_stripe(interpolator, interpolator%stripe, lev);
     call insert_periodic(interpolator, 2, lev, ind, interpolator%stripe, data);
  end do


end subroutine filter_highest


!!$subroutine filter_oscillations_highest(interpolator, data)
!!$  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
!!$  sll_real64, dimension(:), intent(inout) :: data
!!$  sll_int32 :: lev, l, j, ind
!!$
!!$  lev = interpolator%levels(1);
!!$  do l = 0, interpolator%max_level-lev
!!$     ind = interpolator%index(0,l)
!!$     call extract_periodic(interpolator, 1, lev, ind, data, interpolator%stripe)
!!$     !call hierarchical_stripe(interpolator, interpolator%stripe, lev);
!!$     do j=1,2**lev-1
!!$        if ( interpolator%stripe(j)>interpolator%stripe(j+1) ) then
!!$           interpolator%stripe_out(j) = -1;
!!$        else
!!$           interpolator%stripe_out(j) = 1;
!!$        end if
!!$     end do
!!$     p = 1;
!!$     do j=1,2**lev
!!$        if (interpolator%stripe_out(j) == interpolator%stripe_out(j+1)) then
!!$           damp
!!$           p = j+1;
!!$        end if
!!$     end do
!!$     do j=2,2**lev,2
!!$        interpolator%stripe_out(j) = sign(interpolator%stripe(j),1);
!!$     end do
!!$     do j
!!$        interpolator%stripe(j) = interpolator%stripe(j)/2;
!!$     end do
!!$     call dehierarchical_stripe(interpolator, interpolator%stripe, lev);
!!$     call insert_periodic(interpolator, 1, lev, ind, interpolator%stripe, data);
!!$  end do
!!$
!!$  lev = interpolator%levels(2);
!!$  do l = 0, interpolator%max_level-lev
!!$     ind = interpolator%index(l,0)
!!$     call extract_periodic(interpolator, 2, lev, ind, data, interpolator%stripe)
!!$     call hierarchical_stripe(interpolator, interpolator%stripe, lev);
!!$     do j=2,2**lev,2
!!$        interpolator%stripe(j) = interpolator%stripe(j)*0.5_f64;
!!$     end do
!!$     call dehierarchical_stripe(interpolator, interpolator%stripe, lev);
!!$     call insert_periodic(interpolator, 2, lev, ind, interpolator%stripe, data);
!!$  end do
!!$
!!$
!!$end subroutine filter_oscillations_highest


subroutine filter(interpolator, data)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data
  sll_int32 :: j

  call compute_hierarchical_surplus(interpolator,data);
  do j=interpolator%level_mapping(interpolator%max_level)+1, interpolator%size_basis,2
     data(j) = data(j)*0.5_f64;
  end do
  call compute_dehierarchical(interpolator,data);

end subroutine filter

subroutine filter_oscillations(interpolator, data)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data
  sll_int32 :: j

  call compute_hierarchical_surplus(interpolator,data);
  do j=interpolator%level_mapping(interpolator%max_level)+1, interpolator%size_basis,2
     data(j) = data(j)*0.5_f64;
  end do
  call compute_dehierarchical(interpolator,data);

end subroutine filter_oscillations


subroutine linear_filter(interpolator, data, hs, width)
  class(sparse_grid_interpolator_2d), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data, hs
  sll_real64, dimension(2), intent(in) :: width
  sll_real64, dimension(2) :: dx, dxn
  sll_int32 :: j

  do j=1,interpolator%size_basis
     dx = interpolator%hierarchy(j)%coordinate
     dxn(1) = dx(1); dxn(2) = dx(2) + width(2)
     data(j) = data(j)*0.5_f64 + 0.125*interpolator%interpolate_value(hs,dxn);
     dxn(2) = dx(2) - width(2)
     data(j) = data(j) +0.125*interpolator%interpolate_value(hs,dxn);
     dxn(1) = dx(1)+width(1); dxn(2) = dx(2)
     data(j) = data(j) + 0.125*interpolator%interpolate_value(hs,dxn);
     dxn(1) = dx(1) - width(1)
     data(j) = data(j) + 0.125*interpolator%interpolate_value(hs,dxn);
  end do
  


end subroutine linear_filter



!!!! End filter functions !!!!
!------------------------------------------------------------------------------!


end module sparse_grid_2d
