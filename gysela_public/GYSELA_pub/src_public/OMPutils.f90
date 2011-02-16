!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!-------------------------------------------------------
! file : OMPutils.f90
! date : 20/03/2003
!  subroutines used for the OpenMP parallelization
!---------------------------------------------
module OMPutils_module
  use prec_const
  implicit none
      
#ifdef _OPENMP
  integer, public :: omp_get_thread_num, omp_get_num_threads
  external omp_get_thread_num, omp_get_num_threads
#endif
      
  type, public :: omp1dvector
     real(RKIND), dimension(:), pointer :: val
  end type omp1dvector
  type, public :: omp1cvector
     complex(CKIND), dimension(:), pointer :: val
  end type omp1cvector
  type, public :: omp2dvector
     real(RKIND), dimension(:,:), pointer :: val
  end type omp2dvector
  type, public :: omp2cvector
     complex(CKIND), dimension(:,:), pointer :: val
  end type omp2cvector
  type, public :: omp2ivector
     integer, dimension(:,:), pointer :: val
  end type omp2ivector
      
  public
  !-> real 1D arrays
  !--> radial direction
  type(omp1dvector), dimension (:), pointer, public :: Romp1_0Nr
  type(omp1dvector), dimension (:), pointer, public :: Romp2_0Nr
  type(omp1dvector), dimension (:), pointer, public :: Romp3_0Nr
  type(omp1dvector), dimension (:), pointer, public :: Romp4_0Nr
  type(omp1dvector), dimension (:), pointer, public :: Romp5_0Nr
  type(omp1dvector), dimension (:), pointer, public :: Romp6_0Nr
  type(omp1dvector), dimension (:), pointer, public :: Romp1_1Nrp1
  type(omp1dvector), dimension (:), pointer, public :: Romp2_1Nrp1
  type(omp1dvector), dimension (:), pointer, public :: Romp3_1Nrp1
  type(omp1dvector), dimension (:), pointer, public :: Romp4_1Nrp1
  !--> theta direction 
  type(omp1dvector), dimension (:), pointer, public :: Romp1_0Ntheta
  type(omp1dvector), dimension (:), pointer, public :: Romp2_0Ntheta
  type(omp1dvector), dimension (:), pointer, public :: Romp3_0Ntheta
  type(omp1dvector), dimension (:), pointer, public :: Romp4_0Ntheta
  type(omp1dvector), dimension (:), pointer, public :: Romp5_0Ntheta
  type(omp1dvector), dimension (:), pointer, public :: Romp6_0Ntheta
  type(omp1dvector), dimension (:), pointer, public :: &
    Romp1_1Nthetap1
  type(omp1dvector), dimension (:), pointer, public :: &
    Romp2_1Nthetap1
  !--> phi direction 
  type(omp1dvector), dimension (:), pointer, public :: Romp1_0Nphi
  type(omp1dvector), dimension (:), pointer, public :: Romp2_0Nphi
  type(omp1dvector), dimension (:), pointer, public :: Romp3_0Nphi
  type(omp1dvector), dimension (:), pointer, public :: Romp4_0Nphi
  type(omp1dvector), dimension (:), pointer, public :: Romp5_0Nphi
  type(omp1dvector), dimension (:), pointer, public :: Romp6_0Nphi
  !--> vparallel direction
  type(omp1dvector), dimension (:), pointer, public :: Romp1_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp2_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp3_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp4_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp5_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp6_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp7_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp8_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp9_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp10_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp11_0Nvpar
  type(omp1dvector), dimension (:), pointer, public :: Romp12_0Nvpar
  !--> used for 2D advections
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_vect_mpstencil
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_fillde_112
  !--> used for cubic splines
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_deriv_01
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_scoefr
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_scoeftheta
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_scoefphi
  type(omp1dvector), dimension (:), pointer, &
    public :: Romp_scoefvpar
  !-> complex 1D arrays
  type(omp1cvector), dimension (:), pointer, &
    public :: Comp_0Nr
      
  !-> real 2D arrays 
  !----> (Nr,Ntheta)
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp1_0Nr_0Ntheta
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp2_0Nr_0Ntheta
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp1_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp2_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp3_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp4_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp5_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp6_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp7_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp8_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp9_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp10_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp11_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp12_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp13_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp14_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp15_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp16_1Nrp1_1Nthetap1
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp17_1Nrp1_1Nthetap1
  !----> (1:Ntheta+1,1:Nr+1)
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp1_1Ntheta_1Nrp1
  !----> (1:Nr+buf,1:Ntheta+buf)
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp1_1Nrpb_1Nthetapb
  !----> (1:nbmoments,Nr)
  type(omp2dvector), dimension (:), pointer, &
    public :: Romp1_1nbmom_0Nr
      
  !-> complex 2D arrays 
  !----> (1:Nr+1,1:Ntheta+1)
  type(omp2cvector), dimension (:), pointer, &
    public :: Comp1_1Nrp1_1Ntheta
  type(omp2cvector), dimension (:), pointer, &
    public :: Comp2_1Nrp1_1Ntheta
  !----> (1:Ntheta,1:Nr+1)
  type(omp2cvector), dimension (:), pointer, &
    public :: Comp1_1Ntheta_1Nrp1
  !----> (1:Ntheta+buf,1:Nr+buf)
  type(omp2cvector), dimension (:), pointer, &
    public :: Comp1_1Nthetapb_1Nrpb
  !----> (1:Nr+buf,1:Ntheta+buf)
  type(omp2cvector), dimension (:), pointer, &
    public :: Comp1_1Nrpb_1Nthetapb
      
  integer, private :: max_dim
      
  !******************************
  contains
  !******************************   
      
  !-----------------------------------------------------
  ! Constructor
  !-----------------------------------------------------
  subroutine init_OMP
    use globals, only : Nr, Ntheta, Nphi, Nvpar, Nmu, &
      Nbthread, stencil, nbmoments_GC, nbmoments_part
    use mem_alloc_module
      
    integer :: tid
    integer :: nbmoments
      
    nbmoments = max(nbmoments_GC,nbmoments_part)
      
    !*** structure allocations ***
    !-> allocation of structures for real 1D arrays
    !--> radial direction
    allocate(Romp1_0Nr(1:Nbthread))
    allocate(Romp2_0Nr(1:Nbthread))
    allocate(Romp3_0Nr(1:Nbthread))
    allocate(Romp4_0Nr(1:Nbthread))
    allocate(Romp5_0Nr(1:Nbthread))
    allocate(Romp6_0Nr(1:Nbthread))
    allocate(Romp1_1Nrp1(1:Nbthread))
    allocate(Romp2_1Nrp1(1:Nbthread))
    allocate(Romp3_1Nrp1(1:Nbthread))
    allocate(Romp4_1Nrp1(1:Nbthread))
    !-> theta direction
    allocate(Romp1_0Ntheta(1:Nbthread))
    allocate(Romp2_0Ntheta(1:Nbthread))
    allocate(Romp3_0Ntheta(1:Nbthread))
    allocate(Romp4_0Ntheta(1:Nbthread))
    allocate(Romp5_0Ntheta(1:Nbthread))
    allocate(Romp6_0Ntheta(1:Nbthread))
    allocate(Romp1_1Nthetap1(1:Nbthread))
    allocate(Romp2_1Nthetap1(1:Nbthread))
    !-> phi direction
    allocate(Romp1_0Nphi(1:Nbthread))
    allocate(Romp2_0Nphi(1:Nbthread))
    allocate(Romp3_0Nphi(1:Nbthread))
    allocate(Romp4_0Nphi(1:Nbthread))
    allocate(Romp5_0Nphi(1:Nbthread))
    allocate(Romp6_0Nphi(1:Nbthread))
    !-> vparallel direction
    allocate(Romp1_0Nvpar(1:Nbthread))
    allocate(Romp2_0Nvpar(1:Nbthread))
    allocate(Romp3_0Nvpar(1:Nbthread))
    allocate(Romp4_0Nvpar(1:Nbthread))
    allocate(Romp5_0Nvpar(1:Nbthread))
    allocate(Romp6_0Nvpar(1:Nbthread))
    allocate(Romp7_0Nvpar(1:Nbthread))
    allocate(Romp8_0Nvpar(1:Nbthread))
    allocate(Romp9_0Nvpar(1:Nbthread))
    allocate(Romp10_0Nvpar(1:Nbthread))
    allocate(Romp11_0Nvpar(1:Nbthread))
    allocate(Romp12_0Nvpar(1:Nbthread))
    !-> used for 2D advections
    allocate(Romp_vect_mpstencil(1:Nbthread))
    allocate(Romp_fillde_112(1:Nbthread))
    !--> used for cubic splines
    allocate(Romp_deriv_01(1:Nbthread))
    allocate(Romp_scoefr(1:Nbthread))
    allocate(Romp_scoeftheta(1:Nbthread))
    allocate(Romp_scoefphi(1:Nbthread))
    allocate(Romp_scoefvpar(1:Nbthread))
    !-> allocation of structures for complex 1D arrays
    allocate(Comp_0Nr(1:Nbthread))
    !-> allocation of structures for real 2D arrays (Nr,Ntheta)
    allocate(Romp1_0Nr_0Ntheta(1:Nbthread))
    allocate(Romp2_0Nr_0Ntheta(1:Nbthread))
    allocate(Romp1_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp2_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp3_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp4_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp5_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp6_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp7_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp8_1Nrp1_1Nthetap1(1:Nbthread))
    allocate(Romp9_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp10_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp11_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp12_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp13_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp14_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp15_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp16_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp17_1Nrp1_1Nthetap1(0:MAX(Nmu,Nbthread)))
    allocate(Romp1_1Ntheta_1Nrp1(1:Nbthread))
    allocate(Romp1_1Nrpb_1Nthetapb(1:Nbthread))
    allocate(Romp1_1nbmom_0Nr(1:Nbthread))
    !-> allocation of structures for complex 2D arrays (Nr,Ntheta)
    allocate(Comp1_1Nrp1_1Ntheta(1:Nbthread))
    allocate(Comp2_1Nrp1_1Ntheta(1:Nbthread))
    allocate(Comp1_1Nrpb_1Nthetapb(1:Nbthread))
    allocate(Comp1_1Nthetapb_1Nrpb(1:Nbthread))
    allocate(Comp1_1Ntheta_1Nrp1(1:Nbthread))
      
    !*** array allocations ***
    !-> real 1D arrays
    do tid = 1,Nbthread
      !--> radial direction
      call glob_allocate(Romp1_0Nr(tid)%val,0,Nr,'Romp1_0Nr')
      call glob_allocate(Romp2_0Nr(tid)%val,0,Nr,'Romp2_0Nr')
      call glob_allocate(Romp3_0Nr(tid)%val,0,Nr,'Romp3_0Nr')
      call glob_allocate(Romp4_0Nr(tid)%val,0,Nr,'Romp4_0Nr')
      call glob_allocate(Romp5_0Nr(tid)%val,0,Nr,'Romp5_0Nr')
      call glob_allocate(Romp6_0Nr(tid)%val,0,Nr,'Romp6_0Nr')
      call glob_allocate(Romp1_1Nrp1(tid)%val,1,Nr+1,'Romp1_1Nrp1')
      call glob_allocate(Romp2_1Nrp1(tid)%val,1,Nr+1,'Romp2_1Nrp1')
      call glob_allocate(Romp3_1Nrp1(tid)%val,1,Nr+1,'Romp3_1Nrp1')
      call glob_allocate(Romp4_1Nrp1(tid)%val,1,Nr+1,'Romp4_1Nrp1')
      !--> theta direction 
      call glob_allocate(Romp1_0Ntheta(tid)%val,0,Ntheta, &
        'Romp1_0Ntheta')
      call glob_allocate(Romp2_0Ntheta(tid)%val,0,Ntheta, &
        'Romp2_0Ntheta')
      call glob_allocate(Romp3_0Ntheta(tid)%val,0,Ntheta, &
        'Romp3_0Ntheta')
      call glob_allocate(Romp4_0Ntheta(tid)%val,0,Ntheta, &
        'Romp4_0Ntheta')
      call glob_allocate(Romp5_0Ntheta(tid)%val,0,Ntheta, &
        'Romp5_0Ntheta')
      call glob_allocate(Romp6_0Ntheta(tid)%val,0,Ntheta, &
        'Romp6_0Ntheta')
      call glob_allocate(Romp1_1Nthetap1(tid)%val,1,Ntheta+1, &
        'Romp1_1Nthetap1')
      call glob_allocate(Romp2_1Nthetap1(tid)%val,1,Ntheta+1, &
        'Romp2_1Nthetap1')
      !--> phi direction 
      call glob_allocate(Romp1_0Nphi(tid)%val,0,Nphi,'Romp1_0Nphi')
      call glob_allocate(Romp2_0Nphi(tid)%val,0,Nphi,'Romp2_0Nphi')
      call glob_allocate(Romp3_0Nphi(tid)%val,0,Nphi,'Romp3_0Nphi')
      call glob_allocate(Romp4_0Nphi(tid)%val,0,Nphi,'Romp4_0Nphi')
      call glob_allocate(Romp5_0Nphi(tid)%val,0,Nphi,'Romp5_0Nphi')
      call glob_allocate(Romp6_0Nphi(tid)%val,0,Nphi,'Romp6_0Nphi')
      !--> vparallel direction 
      call glob_allocate(Romp1_0Nvpar(tid)%val,0,Nvpar, &
        'Romp1_0Nvpar')
      call glob_allocate(Romp2_0Nvpar(tid)%val,0,Nvpar, &
        'Romp2_0Nvpar')
      call glob_allocate(Romp3_0Nvpar(tid)%val,0,Nvpar, &
        'Romp3_0Nvpar')
      call glob_allocate(Romp4_0Nvpar(tid)%val,0,Nvpar, &
        'Romp4_0Nvpar')
      call glob_allocate(Romp5_0Nvpar(tid)%val,0,Nvpar, &
        'Romp5_0Nvpar')
      call glob_allocate(Romp6_0Nvpar(tid)%val,0,Nvpar, &
        'Romp6_0Nvpar')
      call glob_allocate(Romp7_0Nvpar(tid)%val,0,Nvpar, &
        'Romp7_0Nvpar')
      call glob_allocate(Romp8_0Nvpar(tid)%val,0,Nvpar, &
        'Romp8_0Nvpar')
      call glob_allocate(Romp9_0Nvpar(tid)%val,0,Nvpar, &
        'Romp9_0Nvpar')
      call glob_allocate(Romp10_0Nvpar(tid)%val,0,Nvpar, &
        'Romp10_0Nvpar')
      call glob_allocate(Romp11_0Nvpar(tid)%val,0,Nvpar, &
        'Romp11_0Nvpar')
      call glob_allocate(Romp12_0Nvpar(tid)%val,0,Nvpar, &
        'Romp12_0Nvpar')
      !--> used for 2D advections
      call glob_allocate(Romp_fillde_112(tid)%val,1,12, &
        'Romp_fillde_112')
      call glob_allocate(Romp_vect_mpstencil(tid)%val, &
        -stencil,stencil,'Romp_vect_mpstencil')
      !--> used for cubic splines
      call glob_allocate(Romp_deriv_01(tid)%val,0,1, &
        'Romp_deriv_01')
      call glob_allocate(Romp_scoefr(tid)%val,-1,Nr+1, &
        'Romp_scoefr')
      call glob_allocate(Romp_scoeftheta(tid)%val,-1,Ntheta+1, &
        'Romp_scoeftheta')
      call glob_allocate(Romp_scoefphi(tid)%val,-1,Nphi+1, &
        'Romp_scoefphi')
      call glob_allocate(Romp_scoefvpar(tid)%val,-1,Nvpar+1, &
        'Romp_scoefvpar')
    end do
      
    !-> complex 1D arrays
    do tid = 1,Nbthread
      call glob_allocate(Comp_0Nr(tid)%val,0,Nr,'Comp_0Nr')
    end do
      
    !-> real 2D arrays (Nr,Ntheta)
    do tid = 0, MAX(Nmu,Nbthread)
      call glob_allocate(Romp9_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp9_1Nrp1_1Nthetap1')
      call glob_allocate(Romp10_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp10_1Nrp1_1Nthetap1')
      call glob_allocate(Romp11_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp11_1Nrp1_1Nthetap1')
      call glob_allocate(Romp12_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp12_1Nrp1_1Nthetap1')
      call glob_allocate(Romp13_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp13_1Nrp1_1Nthetap1')
      call glob_allocate(Romp14_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp14_1Nrp1_1Nthetap1')
      call glob_allocate(Romp15_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp15_1Nrp1_1Nthetap1')
      call glob_allocate(Romp16_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp16_1Nrp1_1Nthetap1')
      call glob_allocate(Romp17_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp17_1Nrp1_1Nthetap1')
    end do
    do tid = 1,Nbthread
      call glob_allocate(Romp1_0Nr_0Ntheta(tid)%val,0,Nr,0,Ntheta, &
        'Romp1_0Nr_0Ntheta')
      call glob_allocate(Romp2_0Nr_0Ntheta(tid)%val,0,Nr,0,Ntheta, &
        'Romp2_0Nr_0Ntheta')
      call glob_allocate(Romp1_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp1_1Nrp1_1Nthetap1')
      call glob_allocate(Romp2_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp2_1Nrp1_1Nthetap1')
      call glob_allocate(Romp3_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp3_1Nrp1_1Nthetap1')
      call glob_allocate(Romp4_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp4_1Nrp1_1Nthetap1')
      call glob_allocate(Romp5_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp5_1Nrp1_1Nthetap1')
      call glob_allocate(Romp6_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp6_1Nrp1_1Nthetap1')
      call glob_allocate(Romp7_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp7_1Nrp1_1Nthetap1')
      call glob_allocate(Romp8_1Nrp1_1Nthetap1(tid)%val,1,Nr+1, &
        1,Ntheta+1,'Romp8_1Nrp1_1Nthetap1')
      call glob_allocate(Romp1_1Nrpb_1Nthetapb(tid)%val, &
        1,2*(Nr+1),1,Ntheta+1,'Romp1_1Nrpb_1Nthetapb')
      call glob_allocate(Romp1_1Ntheta_1Nrp1(tid)%val,1,Ntheta, &
        1,Nr+1,'Romp1_1Nthetap1_1Nrp1')
      call glob_allocate(Romp1_1nbmom_0Nr(tid)%val,1,nbmoments, &
        0,Nr,'Romp1_1nbmom_0Nr')
    end do
      
    !-> complex 2D arrays (Nr,Ntheta)
    do tid = 1,Nbthread
      call glob_allocate(Comp1_1Nrp1_1Ntheta(tid)%val,1,Nr+1, &
        1,Ntheta,'Comp1_1Nrp1_1Ntheta')
      call glob_allocate(Comp2_1Nrp1_1Ntheta(tid)%val,1,Nr+1, &
        1,Ntheta,'Comp2_1Nrp1_1Ntheta')
      call glob_allocate(Comp1_1Ntheta_1Nrp1(tid)%val, &
        1,Ntheta,1,Nr+1,'Comp1_1Ntheta_1Nrp1')
      call glob_allocate(Comp1_1Nrpb_1Nthetapb(tid)%val, &
        1,2*(Nr+1),1,Ntheta,'Comp1_1Nrpb_1Nthetapb')
      call glob_allocate(Comp1_1Nthetapb_1Nrpb(tid)%val,1,Ntheta,&
        1,2*(Nr+1),'Comp1_1Nthetapb_1Nrpb')
    end do
  end subroutine init_OMP
      
  !-----------------------------------------------------
  ! Destructor
  !-----------------------------------------------------
  subroutine del_OMP
    use globals, only : Nmu, Nbthread
    use mem_alloc_module
      
    integer :: tid
      
    !*** array deallocations ***
    !-> 1D arrays
    do tid = 1,Nbthread
      !--> real 1D arrays
      !---> radial direction
      call glob_deallocate(Romp1_0Nr(tid)%val)
      call glob_deallocate(Romp2_0Nr(tid)%val)
      call glob_deallocate(Romp3_0Nr(tid)%val)
      call glob_deallocate(Romp4_0Nr(tid)%val)
      call glob_deallocate(Romp5_0Nr(tid)%val)
      call glob_deallocate(Romp6_0Nr(tid)%val)
      call glob_deallocate(Romp1_1Nrp1(tid)%val)
      call glob_deallocate(Romp2_1Nrp1(tid)%val)
      call glob_deallocate(Romp3_1Nrp1(tid)%val)
      call glob_deallocate(Romp4_1Nrp1(tid)%val)
      !---> theta direction
      call glob_deallocate(Romp1_0Ntheta(tid)%val)
      call glob_deallocate(Romp2_0Ntheta(tid)%val)
      call glob_deallocate(Romp3_0Ntheta(tid)%val)
      call glob_deallocate(Romp4_0Ntheta(tid)%val)
      call glob_deallocate(Romp5_0Ntheta(tid)%val)
      call glob_deallocate(Romp6_0Ntheta(tid)%val)
      call glob_deallocate(Romp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp2_1Nthetap1(tid)%val)
      !---> phi direction
      call glob_deallocate(Romp1_0Nphi(tid)%val)
      call glob_deallocate(Romp2_0Nphi(tid)%val)
      call glob_deallocate(Romp3_0Nphi(tid)%val)
      call glob_deallocate(Romp4_0Nphi(tid)%val)
      call glob_deallocate(Romp5_0Nphi(tid)%val)
      call glob_deallocate(Romp6_0Nphi(tid)%val)
      !---> vparallel direction
      call glob_deallocate(Romp1_0Nvpar(tid)%val)
      call glob_deallocate(Romp2_0Nvpar(tid)%val)
      call glob_deallocate(Romp3_0Nvpar(tid)%val)
      call glob_deallocate(Romp4_0Nvpar(tid)%val)
      call glob_deallocate(Romp5_0Nvpar(tid)%val)
      call glob_deallocate(Romp6_0Nvpar(tid)%val)
      call glob_deallocate(Romp7_0Nvpar(tid)%val)
      call glob_deallocate(Romp8_0Nvpar(tid)%val)
      call glob_deallocate(Romp9_0Nvpar(tid)%val)
      call glob_deallocate(Romp10_0Nvpar(tid)%val)
      call glob_deallocate(Romp11_0Nvpar(tid)%val)
      call glob_deallocate(Romp12_0Nvpar(tid)%val)
      !-> used for 2D advections
      call glob_deallocate(Romp_vect_mpstencil(tid)%val)
      call glob_deallocate(Romp_fillde_112(tid)%val)
      !--> used for cubic splines
      call glob_deallocate(Romp_deriv_01(tid)%val)
      call glob_deallocate(Romp_scoefr(tid)%val)
      call glob_deallocate(Romp_scoeftheta(tid)%val)
      call glob_deallocate(Romp_scoefphi(tid)%val)
      call glob_deallocate(Romp_scoefvpar(tid)%val)
      !--> complex 1D arrays
      call glob_deallocate(Comp_0Nr(tid)%val)
    end do
      
    !-> 2D arrays (Nr,Ntheta)
    !--> real 2D arrays (Nr,Ntheta)
    do tid = 1,Nbthread
      call glob_deallocate(Romp1_0Nr_0Ntheta(tid)%val)
      call glob_deallocate(Romp2_0Nr_0Ntheta(tid)%val)
      call glob_deallocate(Romp1_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp2_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp3_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp4_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp5_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp6_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp7_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp8_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp1_1Ntheta_1Nrp1(tid)%val)
      call glob_deallocate(Romp1_1Nrpb_1Nthetapb(tid)%val)
      call glob_deallocate(Romp1_1nbmom_0Nr(tid)%val)
    end do
    do tid = 0, MAX(Nmu,Nbthread)
      call glob_deallocate(Romp9_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp10_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp11_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp12_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp13_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp14_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp15_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp16_1Nrp1_1Nthetap1(tid)%val)
      call glob_deallocate(Romp17_1Nrp1_1Nthetap1(tid)%val)
    end do
      
    !-> complex 2D arrays (Nr,Ntheta)
    do tid = 1,Nbthread
      call glob_deallocate(Comp1_1Nrp1_1Ntheta(tid)%val)
      call glob_deallocate(Comp2_1Nrp1_1Ntheta(tid)%val)
      call glob_deallocate(Comp1_1Ntheta_1Nrp1(tid)%val)
      call glob_deallocate(Comp1_1Nrpb_1Nthetapb(tid)%val)
      call glob_deallocate(Comp1_1Nthetapb_1Nrpb(tid)%val)
    end do
      
    !*** structure deallocations ***
    !-> deallocation of structures for real 1D arrays
    !--> radial direction
    deallocate(Romp1_0Nr)
    deallocate(Romp2_0Nr)
    deallocate(Romp3_0Nr)
    deallocate(Romp4_0Nr)
    deallocate(Romp5_0Nr)
    deallocate(Romp6_0Nr)
    deallocate(Romp1_1Nrp1)
    deallocate(Romp2_1Nrp1)
    deallocate(Romp3_1Nrp1)
    deallocate(Romp4_1Nrp1)
    !-> theta direction
    deallocate(Romp1_0Ntheta)
    deallocate(Romp2_0Ntheta)
    deallocate(Romp3_0Ntheta)
    deallocate(Romp4_0Ntheta)
    deallocate(Romp5_0Ntheta)
    deallocate(Romp6_0Ntheta)
    deallocate(Romp1_1Nthetap1)
    deallocate(Romp2_1Nthetap1)
    !-> phi direction
    deallocate(Romp1_0Nphi)
    deallocate(Romp2_0Nphi)
    deallocate(Romp3_0Nphi)
    deallocate(Romp4_0Nphi)
    deallocate(Romp5_0Nphi)
    deallocate(Romp6_0Nphi)
    !-> vparallel direction
    deallocate(Romp1_0Nvpar)
    deallocate(Romp2_0Nvpar)
    deallocate(Romp3_0Nvpar)
    deallocate(Romp4_0Nvpar)
    deallocate(Romp5_0Nvpar)
    deallocate(Romp6_0Nvpar)
    deallocate(Romp7_0Nvpar)
    deallocate(Romp8_0Nvpar)
    deallocate(Romp9_0Nvpar)
    deallocate(Romp10_0Nvpar)
    deallocate(Romp11_0Nvpar)
    deallocate(Romp12_0Nvpar)
    !-> used for 2D advections
    deallocate(Romp_vect_mpstencil)
    deallocate(Romp_fillde_112)
    !--> used for cubic splines
    deallocate(Romp_deriv_01)
    deallocate(Romp_scoefr)
    deallocate(Romp_scoeftheta)
    deallocate(Romp_scoefphi)
    deallocate(Romp_scoefvpar)
    !-> deallocation of structures for complex 1D arrays
    deallocate(Comp_0Nr)
    !-> deallocation of structures for real 2D arrays (Nr,Ntheta)
    deallocate(Romp1_0Nr_0Ntheta)
    deallocate(Romp2_0Nr_0Ntheta)
    deallocate(Romp1_1Nrp1_1Nthetap1)
    deallocate(Romp2_1Nrp1_1Nthetap1)
    deallocate(Romp3_1Nrp1_1Nthetap1)
    deallocate(Romp4_1Nrp1_1Nthetap1)
    deallocate(Romp5_1Nrp1_1Nthetap1)
    deallocate(Romp6_1Nrp1_1Nthetap1)
    deallocate(Romp7_1Nrp1_1Nthetap1)
    deallocate(Romp8_1Nrp1_1Nthetap1)
    deallocate(Romp9_1Nrp1_1Nthetap1)
    deallocate(Romp10_1Nrp1_1Nthetap1)
    deallocate(Romp11_1Nrp1_1Nthetap1)
    deallocate(Romp12_1Nrp1_1Nthetap1)
    deallocate(Romp13_1Nrp1_1Nthetap1)
    deallocate(Romp14_1Nrp1_1Nthetap1)
    deallocate(Romp15_1Nrp1_1Nthetap1)
    deallocate(Romp16_1Nrp1_1Nthetap1)
    deallocate(Romp17_1Nrp1_1Nthetap1)
    deallocate(Romp1_1Ntheta_1Nrp1)
    deallocate(Romp1_1Nrpb_1Nthetapb)
    deallocate(Romp1_1nbmom_0Nr)
      
    !-> deallocation of structures for complex 2D arrays (Nr,Ntheta)
    deallocate(Comp1_1Nrp1_1Ntheta)
    deallocate(Comp2_1Nrp1_1Ntheta)
    deallocate(Comp1_1Ntheta_1Nrp1)
    deallocate(Comp1_1Nrpb_1Nthetapb)
    deallocate(Comp1_1Nthetapb_1Nrpb)
  end subroutine del_OMP
end module OMPutils_module
