program test_random_deviate_generators
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_gaussian, only: &
    sll_f_gaussian_deviate

  use sll_m_hammersley, only: &
    sll_f_suite_hamm

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define Xmin  0._f64
#define Xmax  5._f64
#define num_particles 8000_i32
#define the_mean  -1._f64
#define the_spatial_variance  4._f64


  sll_int32 :: j, incr_index
  sll_real64 :: pseudo_ran_Hamm_zero_one, pseudo_ran_Hamm_zero_five
  sll_real64 :: pseudo_ran_Hamm_modif_zero_one, z
  sll_real64 :: shift
  logical :: is_passed = .TRUE.


  open(152,file='test_Hamm_Gauss.dat')

  shift = 0.0734_f64
  incr_index = 0
  do j=1, num_particles
     pseudo_ran_Hamm_zero_one = sll_f_suite_hamm(j, 2)
     pseudo_ran_Hamm_zero_five = (Xmax - Xmin) * sll_f_suite_hamm(j, 2) + Xmin! It replaces the 'call random(y)'
     ! This gives a pseudo random sequence between Xmin and Xmax

     z = sll_f_gaussian_deviate() * the_spatial_variance + (the_mean)
     ! This gives a Gaussian centered in 'the_mean'
     pseudo_ran_Hamm_modif_zero_one = sll_f_suite_hamm(j+incr_index, 2) - shift
     do while (pseudo_ran_Hamm_modif_zero_one < 0._f64)
        incr_index = incr_index + 1
        pseudo_ran_Hamm_modif_zero_one = sll_f_suite_hamm(j+incr_index, 2) - shift
     enddo
      pseudo_ran_Hamm_modif_zero_one =  pseudo_ran_Hamm_modif_zero_one / (1._f64 - shift)


     write(152,*) j, pseudo_ran_Hamm_zero_one, pseudo_ran_Hamm_zero_five, &
                  z, pseudo_ran_Hamm_modif_zero_one
  enddo
  close(152)


  write(51,*) "====="
  write(*,*) "-------"
  write(*,*) "End of Program 'unit_test' for 'deviators' "
  if (is_passed) then
     print *, 'PASSED'
  else
     print *, 'FAILED'
  end if
end program test_random_deviate_generators
