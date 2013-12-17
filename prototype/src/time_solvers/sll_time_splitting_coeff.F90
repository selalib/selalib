!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
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
!**************************************************************
!
! Define coefficients for time splitting
!
! contact: Michel Mehrenberger (mehrenbe@math.unistra.fr)
!

module sll_time_splitting_coeff_module


#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  implicit none

  sll_int32, parameter :: SLL_USER_DEFINED         = -3 
  sll_int32, parameter :: SLL_LIE_TV               = -2 
  sll_int32, parameter :: SLL_LIE_VT               = -1 
  sll_int32, parameter :: SLL_STRANG_TVT           = 0 
  sll_int32, parameter :: SLL_STRANG_VTV           = 1 
  sll_int32, parameter :: SLL_TRIPLE_JUMP_TVT      = 2 
  sll_int32, parameter :: SLL_TRIPLE_JUMP_VTV      = 3 
  !sll_int32, parameter :: SLL_ORDER6_TVT           = 4 
  sll_int32, parameter :: SLL_ORDER6_VTV           = 5 
  sll_int32, parameter :: SLL_ORDER6VP_TVT         = 6 
  sll_int32, parameter :: SLL_ORDER6VP_VTV         = 7 
  sll_int32, parameter :: SLL_ORDER6VPnew_TVT      = 8 
  sll_int32, parameter :: SLL_ORDER6VPnew1_VTV     = 9 
  sll_int32, parameter :: SLL_ORDER6VPnew2_VTV     = 10 
  
  type splitting_coeff
    sll_int32 :: split_case
    sll_real64, dimension(:), pointer :: split_step
    sll_int32 :: nb_split_step
    logical :: split_begin_T
  end type splitting_coeff  
  
contains
  function new_time_splitting_coeff( &
    split_case, &
    split_step, &
    nb_split_step, &
    split_begin_T, &
    dt) &
    result(split)  
    type(splitting_coeff), pointer :: split
    sll_int32, intent(in)  :: split_case
    sll_real64, dimension(:), intent(in), optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T
    sll_real64, intent(in), optional :: dt   
    sll_int32 :: ierr
   
    SLL_ALLOCATE(split,ierr)   
    call initialize_time_splitting_coeff( &
      split, &
      split_case, &
      split_step, &
      nb_split_step, &
      split_begin_T, &
      dt)
      
  end function new_time_splitting_coeff  
  
  subroutine initialize_time_splitting_coeff( &
    split, &
    split_case, &
    split_step, &
    nb_split_step, &
    split_begin_T, &
    dt)
    
    type(splitting_coeff) :: split
    sll_int32, intent(in) :: split_case
    sll_real64, dimension(:), intent(in), optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T
    sll_real64, intent(in), optional :: dt   
    sll_int32 :: ierr 

    if(split_case==SLL_USER_DEFINED) then
      if( .not.( (present(split_step)) &
        .and.(present(nb_split_step)) &
        .and.(present(split_begin_T)) )) then
        print *,'#provide split_step, nb_split_step,split_begin_T'
        print *,'#in initialize_time_splitting_coeff'
        stop
      endif  
    else
      if(&
        (present(split_step)) &
        .or.(present(nb_split_step)) &
        .or.(present(split_begin_T)) &
        )then
        print *,'# do not provide split_step, nb_split_step, split_begin_T'
        print *,'#in initialize_time_splitting_coeff'
        stop       
      endif
    endif  
!    if(associated(split%split_step))then
!      print *,'#split%split_step should not be allocated'
!      print *,'#in initialize_time_splitting_coeff'
!      stop
!    endif

    
    select case (split_case)    
      case (SLL_USER_DEFINED)
        if(size(split_step)<nb_split_step) then
          print *,'#bad size for split_step',size(split_step),nb_split_step
          print *,'#in initialize_time_splitting_coeff'
          stop
        else
          SLL_ALLOCATE(split%split_step(nb_split_step),ierr)
          split%split_step(1:nb_split_step)=split_step(1:nb_split_step)
          split%nb_split_step =nb_split_step
          split%split_begin_T = split_begin_T  
        endif
      case (SLL_LIE_TV) 
        split%nb_split_step = 2
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .true.
        split%split_step(1) = 1._f64
        split%split_step(2) = 1._f64
      case (SLL_LIE_VT) 
        split%nb_split_step = 2
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .false.
        split%split_step(1) = 1._f64
        split%split_step(2) = 1._f64
      case (SLL_STRANG_TVT) ! Strang splitting TVT
        split%nb_split_step = 3
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = 1
        split%split_step(1) = 0.5_f64
        split%split_step(2) = 1._f64
        split%split_step(3) = split%split_step(1)
      case (SLL_STRANG_VTV) ! Strang splitting VTV
        split%nb_split_step = 3
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = 0
        split%split_step(1) = 0.5_f64
        split%split_step(2) = 1._f64
        split%split_step(3) = split%split_step(1)
      case (SLL_TRIPLE_JUMP_TVT) ! triple jump TVT
        split%nb_split_step = 7
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = 1
        split%split_step(1) = 0.675603595979829_f64
        split%split_step(2) = 1.351207191959658_f64
        split%split_step(3) = -0.17560359597982855_f64
        split%split_step(4) = -1.702414383919315_f64
        split%split_step(5) = split%split_step(3)
        split%split_step(6) = split%split_step(2)
        split%split_step(7) = split%split_step(1)
      case (SLL_TRIPLE_JUMP_VTV) ! triple jump VTV
        split%nb_split_step = 7
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = 0
        split%split_step(1) = 0.675603595979829_f64
        split%split_step(2) = 1.351207191959658_f64
        split%split_step(3) = -0.17560359597982855_f64
        split%split_step(4) = -1.702414383919315_f64
        split%split_step(5) = split%split_step(3)
        split%split_step(6) = split%split_step(2)
        split%split_step(7) = split%split_step(1)
      case (SLL_ORDER6_VTV) ! Order 6 VTV
        split%nb_split_step = 23
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = 0
        split%split_step(1) = 0.0414649985182624_f64
        split%split_step(2) = 0.123229775946271_f64
        split%split_step(3) = 0.198128671918067_f64
        split%split_step(4) = 0.290553797799558_f64
        split%split_step(5) = -0.0400061921041533_f64
        split%split_step(6) = -0.127049212625417_f64
        split%split_step(7) = 0.0752539843015807_f64          
        split%split_step(8) = -0.246331761062075_f64
        split%split_step(9) = -0.0115113874206879_f64
        split%split_step(10) = 0.357208872795928_f64
        split%split_step(11) = 0.23666992478693111_f64
        split%split_step(12) = 0.20477705429147008_f64
        split%split_step(13) = split%split_step(11)
        split%split_step(14) = split%split_step(10)
        split%split_step(15) = split%split_step(9)          
        split%split_step(16) = split%split_step(8)
        split%split_step(17) = split%split_step(7)
        split%split_step(18) = split%split_step(6)
        split%split_step(19) = split%split_step(5)
        split%split_step(20) = split%split_step(4)
        split%split_step(21) = split%split_step(3)
        split%split_step(22) = split%split_step(2)
        split%split_step(23) = split%split_step(1)          
      case (SLL_ORDER6VP_TVT) ! Order 6 for Vlasov-Poisson TVT 
        if(present(dt))then
          split%nb_split_step = 9
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = 1
          split%split_step(1) = 0.1095115577513980413559540_f64
          split%split_step(2) = 0.268722208204814693684441_f64&
            -2._f64*dt**2*0.000805681667096178271312_f64&
            +4._f64*dt**4*0.000017695766224036466792_f64
          split%split_step(3) = 0.4451715080955340951457244_f64
          split%split_step(4) = 0.2312777917951853063155588_f64&
            -2._f64*dt**2*0.003955911930042478239763_f64&
            +4._f64*dt**4*0.000052384078562246674986_f64
          split%split_step(5) = -0.1093661316938642730033570_f64
          split%split_step(6) = split%split_step(4)
          split%split_step(7) = split%split_step(3)
          split%split_step(8) = split%split_step(2)
          split%split_step(9) = split%split_step(1)
        else
          print *,'#provide dt for use of case SLL_ORDER6VP_TVT=',SLL_ORDER6VP_TVT
          stop
        endif  
      case (SLL_ORDER6VP_VTV) ! Order 6 for Vlasov-Poisson VTV 
        if(present(dt))then        
          split%nb_split_step = 9
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = 0
          split%split_step(1) = 0.359950808794143627485664_f64&
            -2._f64*dt**2*(-0.01359558332625151635_f64)&
            +4._f64*dt**4*(-8.562814848565929e-6_f64)
          split%split_step(2) = 1.079852426382430882456991_f64
          split%split_step(3) = -0.1437147273026540434771131_f64&
            -2._f64*dt**2*(-0.00385637757897273261_f64)&
            +4._f64*dt**4*(0.0004883788785819335822_f64)
          split%split_step(4) = -0.579852426382430882456991_f64
          split%split_step(5) = 0.567527837017020831982899_f64&
            -2._f64*dt**2*(-0.03227361602037480885_f64)&
            +4._f64*dt**4*0.002005141087312622342_f64
          split%split_step(6) = split%split_step(4)
          split%split_step(7) = split%split_step(3)
          split%split_step(8) = split%split_step(2)
          split%split_step(9) = split%split_step(1)
        else
          print *,'#provide dt for use of case SLL_ORDER6VP_VTV=',SLL_ORDER6VP_VTV
          stop
        endif  
      case (SLL_ORDER6VPnew_TVT) ! Order 6 for Vlasov-Poisson TVT (new)
        if(present(dt))then        
          split%nb_split_step = 9
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = 1
          split%split_step(1) = 0.1095115577513980413559540_f64
          split%split_step(2) = 0.268722208204814693684441_f64&
            -2._f64*dt**2*0.000805681667096178271312_f64&
            +4._f64*dt**4*(8.643923349886021963e-6_f64)&
            -8._f64*dt**6*(1.4231479258353431522e-6_f64)
          split%split_step(3) = 0.4451715080955340951457244_f64
          split%split_step(4) = 0.2312777917951853063155588_f64&
            -2._f64*dt**2*0.003955911930042478239763_f64&
            +4._f64*dt**4*(0.000061435921436397119815_f64)
          split%split_step(5) = -0.1093661316938642730033570_f64
          split%split_step(6) = split%split_step(4)
          split%split_step(7) = split%split_step(3)
          split%split_step(8) = split%split_step(2)
          split%split_step(9) = split%split_step(1)
        else
          print *,'#provide dt for use of case SLL_ORDER6VPnew_TVT=',SLL_ORDER6VPnew_TVT
          stop
        endif  
          
      case (SLL_ORDER6VPnew1_VTV) ! Order 6 for Vlasov-Poisson VTV (new1)
        if(present(dt))then        
          split%nb_split_step = 11
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = 0
          split%split_step(1) = 0.0490864609761162454914412_f64&
            -2._f64*dt**2*(0.0000697287150553050840999_f64)
          split%split_step(2) = 0.1687359505634374224481957_f64
          split%split_step(3) = 0.2641776098889767002001462_f64&
            -2._f64*dt**2*(0.000625704827430047189169_f64)&
            +4._f64*dt**4*(-2.91660045768984781644e-6_f64)
          split%split_step(4) = 0.377851589220928303880766_f64
          split%split_step(5) = 0.1867359291349070543084126_f64&
            -2._f64*dt**2*(0.00221308512404532556163_f64)&
            +4._f64*dt**4*(0.0000304848026170003878868_f64)&
            -8._f64*dt**6*(4.98554938787506812159e-7_f64)
          split%split_step(6) = -0.0931750795687314526579244_f64
          split%split_step(7) = split%split_step(5)
          split%split_step(8) = split%split_step(4)
          split%split_step(9) = split%split_step(3)
          split%split_step(10) = split%split_step(2)
          split%split_step(11) = split%split_step(1)
        else
          print *,'#provide dt for use of case SLL_ORDER6VPnew1_VTV=',SLL_ORDER6VPnew1_VTV
          stop
        endif  
      case (SLL_ORDER6VPnew2_VTV) ! Order 6 for Vlasov-Poisson VTV (new2)
        if(present(dt))then        
          split%nb_split_step = 11
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = 0
          split%split_step(1) = 0.083335463273305120964507_f64&
            -2._f64*dt**2*(-0.00015280483587048489661_f64)&
            +4._f64*dt**4*( -0.0017675734111895638156_f64)&
            -8._f64*dt**6*( 0.00021214072262165668039_f64)
          split%split_step(2) = 0.72431592569108212422250_f64
          split%split_step(3) = 0.827694857845135145869413_f64&
            -2._f64*dt**2*(-0.010726848627286273332_f64)&
            +4._f64*dt**4*(0.012324362982853212700_f64)
          split%split_step(4) = -0.4493507217041624582458844_f64
          split%split_step(5) = -0.4110303211184402668339201_f64&
            -2._f64*dt**2*(0.014962337009932798678_f64)
          !  +4._f64*sim%dt**4*(-0.20213028317750837901_f64)
          split%split_step(6) = 0.4500695920261606680467717_f64
          split%split_step(7) = split%split_step(5)
          split%split_step(8) = split%split_step(4)
          split%split_step(9) = split%split_step(3)
          split%split_step(10) = split%split_step(2)
          split%split_step(11) = split%split_step(1)
        else
          print *,'#provide dt for use of case SLL_ORDER6VPnew2_VTV=',SLL_ORDER6VPnew2_VTV
          stop
        endif  
      case default
        print *,'#split_case not defined'
        stop       
    end select
  end subroutine initialize_time_splitting_coeff





end module sll_time_splitting_coeff_module
