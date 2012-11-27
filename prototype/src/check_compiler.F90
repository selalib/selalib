! This program is used to check the compiler and the version 
! We return different exit status to determine the compiler
! To read the following 
! exit status : compiler version
!
! 20 : gfortran <4.6.0
! 21 : gfortran >=4.6.0
!
! 30 : intel <11.0
! 31 : intel >11.0
!
!
!

program check_gcc_version

implicit none

#define gfortran_minimum_version 40600
#define intel_minimum_version 1100
#define ibm_minimum_version 0000

#ifdef __GFORTRAN__
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

PRINT*,'GCC VERSION is ',GCC_VERSION

! Test for GCC < gfortran_minimum_version
#if GCC_VERSION < gfortran_minimum_version
call exit(20)
#else
call exit(21)
#endif
#endif

#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER < intel_minimum_version
call exit(30)
#else
call exit(31)
#endif
#endif

end program check_gcc_version
