#ifndef DOXYGEN_SHOULD_SKIP_THIS

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
!> @author Michel Mehrenberger (mehrenbe@math.unistra.fr)
!

module sll_m_time_splitting_coeff


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_f_new_time_splitting_coeff, &
    sll_p_lie_tv, &
    sll_p_lie_vt, &
    sll_p_order6_tvt, &
    sll_p_order6_vtv, &
    sll_p_order6vp2d_vtv, &
    sll_p_order6vp_tvt, &
    sll_p_order6vp_vtv, &
    sll_p_order6vpnew1_vtv, &
    sll_p_order6vpnew2_vtv, &
    sll_p_order6vpnew_tvt, &
    sll_p_order6vpot_vtv, &
    sll_p_order6vpotnew1_vtv, &
    sll_p_order6vpotnew2_vtv, &
    sll_p_order6vpotnew3_vtv, &
    sll_p_strang_tvt, &
    sll_p_strang_vtv, &
    sll_p_triple_jump_tvt, &
    sll_p_triple_jump_vtv, &
    sll_t_splitting_coeff

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_USER_DEFINED         = -3 
  sll_int32, parameter :: sll_p_lie_tv               = -2 
  sll_int32, parameter :: sll_p_lie_vt               = -1 
  sll_int32, parameter :: sll_p_strang_tvt           = 0 
  sll_int32, parameter :: sll_p_strang_vtv           = 1 
  sll_int32, parameter :: sll_p_triple_jump_tvt      = 2 
  sll_int32, parameter :: sll_p_triple_jump_vtv      = 3 
  !sll_int32, parameter :: sll_p_order6_tvt           = 4 
  sll_int32, parameter :: sll_p_order6_vtv           = 5 
  sll_int32, parameter :: sll_p_order6vp_tvt         = 6 
  sll_int32, parameter :: sll_p_order6vp_vtv         = 7 
  sll_int32, parameter :: sll_p_order6vpnew_tvt      = 8 
  sll_int32, parameter :: sll_p_order6vpnew1_vtv     = 9 
  sll_int32, parameter :: sll_p_order6vpnew2_vtv     = 10 
  sll_int32, parameter :: sll_p_order6vp2d_vtv     = 11 
  sll_int32, parameter :: sll_p_order6vpot_vtv     = 12 
  sll_int32, parameter :: sll_p_order6_tvt           = 13 
  sll_int32, parameter :: sll_p_order6vpotnew1_vtv     = 14 
  sll_int32, parameter :: sll_p_order6vpotnew2_vtv     = 15 
  sll_int32, parameter :: sll_p_order6vpotnew3_vtv     = 16 

  type sll_t_splitting_coeff
    sll_int32 :: split_case
    sll_real64, dimension(:), pointer :: split_step
    sll_int32 :: nb_split_step
    logical :: split_begin_T
    sll_int32 :: dim_split_V
  end type sll_t_splitting_coeff  
  
contains
  function sll_f_new_time_splitting_coeff( &
    split_case, &
    split_step, &
    nb_split_step, &
    split_begin_T, &
    dt) &
    result(split)  
    type(sll_t_splitting_coeff), pointer :: split
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
      
  end function sll_f_new_time_splitting_coeff  
  
  subroutine initialize_time_splitting_coeff( &
    split, &
    split_case, &
    split_step, &
    nb_split_step, &
    split_begin_T, &
    dt)
    
    type(sll_t_splitting_coeff) :: split
    sll_int32, intent(in) :: split_case
    sll_real64, dimension(:), intent(in), optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T
    sll_real64, intent(in), optional :: dt   
    sll_int32 :: ierr 

    split%dim_split_V = 1

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
      case (sll_p_lie_tv) 
        split%nb_split_step = 2
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .true.
        split%split_step(1) = 1._f64
        split%split_step(2) = 1._f64
      case (sll_p_lie_vt) 
        split%nb_split_step = 2
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .false.
        split%split_step(1) = 1._f64
        split%split_step(2) = 1._f64
      case (sll_p_strang_tvt) ! Strang splitting TVT
        split%nb_split_step = 3
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .true.
        split%split_step(1) = 0.5_f64
        split%split_step(2) = 1._f64
        split%split_step(3) = split%split_step(1)
      case (sll_p_strang_vtv) ! Strang splitting VTV
        split%nb_split_step = 3
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .false.
        split%split_step(1) = 0.5_f64
        split%split_step(2) = 1._f64
        split%split_step(3) = split%split_step(1)
      case (sll_p_triple_jump_tvt) ! triple jump TVT
        split%nb_split_step = 7
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .true.
        split%split_step(1) = 0.675603595979829_f64
        split%split_step(2) = 1.351207191959658_f64
        split%split_step(3) = -0.17560359597982855_f64
        split%split_step(4) = -1.702414383919315_f64
        split%split_step(5) = split%split_step(3)
        split%split_step(6) = split%split_step(2)
        split%split_step(7) = split%split_step(1)
      case (sll_p_triple_jump_vtv) ! triple jump VTV
        split%nb_split_step = 7
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .false.
        split%split_step(1) = 0.675603595979829_f64
        split%split_step(2) = 1.351207191959658_f64
        split%split_step(3) = -0.17560359597982855_f64
        split%split_step(4) = -1.702414383919315_f64
        split%split_step(5) = split%split_step(3)
        split%split_step(6) = split%split_step(2)
        split%split_step(7) = split%split_step(1)
      case (sll_p_order6_vtv) ! Order 6 VTV (O6-11 of Blanes)
        split%nb_split_step = 23
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .false.
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
      case (sll_p_order6_tvt) ! Order 6 TVT (O6-14 of Blanes)
        split%nb_split_step = 29
        SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
        split%split_begin_T = .true.
        split%split_step(1) = 0.0378593198406116_f64     
        split%split_step(2) = 0.09171915262446165_f64        
        split%split_step(3) = 0.102635633102435_f64        
        split%split_step(4) = 0.183983170005006_f64        
        split%split_step(5) = -0.0258678882665587_f64      
        split%split_step(6) = -0.05653436583288827_f64      
        split%split_step(7) = 0.314241403071447_f64        
        split%split_step(8) = 0.004914688774712854_f64        
        split%split_step(9) = -0.130144459517415_f64        
        split%split_step(10) = 0.143761127168358_f64        
        split%split_step(11) = 0.106417700369543_f64       
        split%split_step(12) = 0.328567693746804_f64        
        split%split_step(13) = -0.00879424312851058_f64        
        split%split_step(14) = -0.196411466486454234_f64        
        split%split_step(15) = 0.20730506905689536_f64                
        split%split_step(16) = split%split_step(14)
        split%split_step(17) = split%split_step(13)
        split%split_step(18) = split%split_step(12)
        split%split_step(19) = split%split_step(11)
        split%split_step(20) = split%split_step(10)
        split%split_step(21) = split%split_step(9)
        split%split_step(22) = split%split_step(8)
        split%split_step(23) = split%split_step(7)          
        split%split_step(24) = split%split_step(6)          
        split%split_step(25) = split%split_step(5)          
        split%split_step(26) = split%split_step(4)          
        split%split_step(27) = split%split_step(3)          
        split%split_step(28) = split%split_step(2)          
        split%split_step(29) = split%split_step(1)
                  
      case (sll_p_order6vp_tvt) ! Order 6 for Vlasov-Poisson TVT 
        if(present(dt))then
          split%nb_split_step = 9
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = .true.
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
          print *,'#provide dt for use of case sll_p_order6vp_tvt=',sll_p_order6vp_tvt
          stop
        endif  
      case (sll_p_order6vp_vtv) ! Order 6 for Vlasov-Poisson VTV 
        if(present(dt))then        
          split%nb_split_step = 9
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = .false.
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
          print *,'#provide dt for use of case sll_p_order6vp_vtv=',sll_p_order6vp_vtv
          stop
        endif  
      case (sll_p_order6vpnew_tvt) ! Order 6 for Vlasov-Poisson TVT (new)
        if(present(dt))then        
          split%nb_split_step = 9
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = .true.
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
          print *,'#provide dt for use of case sll_p_order6vpnew_tvt=',sll_p_order6vpnew_tvt
          stop
        endif  
          
      case (sll_p_order6vpnew1_vtv) ! Order 6 for Vlasov-Poisson VTV (new1)
        if(present(dt))then        
          split%nb_split_step = 11
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = .false.
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
          print *,'#provide dt for use of case sll_p_order6vpnew1_vtv=',sll_p_order6vpnew1_vtv
          stop
        endif  
      case (sll_p_order6vp2d_vtv) ! Order 6 for Vlasov-Poisson VTV 2D
        if(present(dt))then        
          split%nb_split_step = 11
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = .false.
          split%split_step(1) = 0.0490864609761162454914412_f64&
            +2._f64*dt**2*(0.00166171386175851683711044_f64)
          split%split_step(2) = 0.1687359505634374224481957_f64
          split%split_step(3) = 0.2641776098889767002001462_f64&
            -2._f64*dt**2*(0.00461492847770001641230401_f64)
          split%split_step(4) = 0.377851589220928303880766_f64
          split%split_step(5) = 0.1867359291349070543084126_f64&
            +2._f64*dt**2*(0.0000446959494108217402966857_f64)
          split%split_step(6) = -0.0931750795687314526579244_f64
          split%split_step(7) = split%split_step(5)
          split%split_step(8) = split%split_step(4)
          split%split_step(9) = split%split_step(3)
          split%split_step(10) = split%split_step(2)
          split%split_step(11) = split%split_step(1)
        else
          print *,'#provide dt for use of case sll_p_order6vp2d_vtv=',sll_p_order6vp2d_vtv
          stop
        endif  
      case (sll_p_order6vpot_vtv) ! Order 6 for Vlasov-Poisson VTV 2D with potential modif
        if(present(dt))then        
          split%nb_split_step = 11
          split%dim_split_V = 2
          SLL_ALLOCATE(split%split_step(17),ierr)
          split%split_begin_T = .false.
          split%split_step(1) = 0.0490864609761162454914412_f64&
            +2._f64*dt**2*(0.00166171386175851683711044_f64)
          split%split_step(2) = dt**2*(0.00166171386175851683711044_f64)            
          split%split_step(3) = 0.1687359505634374224481957_f64
          split%split_step(4) = 0.2641776098889767002001462_f64&
            -2._f64*dt**2*(0.00461492847770001641230401_f64)
          split%split_step(5) = -dt**2*(0.00461492847770001641230401_f64)
          split%split_step(6) = 0.377851589220928303880766_f64
          split%split_step(7) = 0.1867359291349070543084126_f64&
            +2._f64*dt**2*(0.0000446959494108217402966857_f64)
          split%split_step(8) = dt**2*(0.0000446959494108217402966857_f64)
          split%split_step(9) = -0.0931750795687314526579244_f64
          split%split_step(10) = split%split_step(7)
          split%split_step(11) = split%split_step(8)          
          split%split_step(12) = split%split_step(6)
          split%split_step(13) = split%split_step(4)
          split%split_step(14) = split%split_step(5)
          split%split_step(15) = split%split_step(3)
          split%split_step(16) = split%split_step(1)
          split%split_step(17) = split%split_step(2)
        else
          print *,'#provide dt for use of case sll_p_order6vpot_vtv=',sll_p_order6vpot_vtv
          stop
        endif  

      case (sll_p_order6vpotnew1_vtv) ! Order 6 for Vlasov-Poisson VTV 2D with potential modif
        !we change also sign here
        if(present(dt))then        
          split%nb_split_step = 9
          split%dim_split_V = 2
          SLL_ALLOCATE(split%split_step(14),ierr)
          split%split_begin_T = .false.
          !b1 a1 b2 a2 b3 a2 b2 a1 b1 
          !b1
          split%split_step(1) = 0.359950808794143627485664_f64&
            +2._f64*dt**2*(0._f64)
          split%split_step(2) = dt**2*(0._f64)            
          !a1
          split%split_step(3) = 1.079852426382430882456991_f64
          !b2
          split%split_step(4) = -0.1437147273026540434771131_f64&
            +2._f64*dt**2*(0.0139652542242388403673_f64)
          split%split_step(5) = dt**2*(0.0139652542242388403673_f64)
          !a2
          split%split_step(6) = -0.579852426382430882456991_f64
          !b3
          split%split_step(7) = 0.567527837017020831982899_f64&
            +2._f64*dt**2*(0.039247029382345626020_f64)
          split%split_step(8) = dt**2*(0.039247029382345626020_f64)
          !a2
          split%split_step(9) = split%split_step(6)
          !b2
          split%split_step(10) = split%split_step(4)
          split%split_step(11) = split%split_step(5)
          !a1
          split%split_step(12) = split%split_step(3)
          !b1
          split%split_step(13) = split%split_step(1)
          split%split_step(14) = split%split_step(2)
        else
          print *,'#provide dt for use of case sll_p_order6vpotnew1_vtv=', &
            sll_p_order6vpotnew1_vtv
          stop
        endif  



      case (sll_p_order6vpotnew2_vtv) ! Order 6 for Vlasov-Poisson VTV 2D with potential modif
        !warning we try with sign change of dt**2
        !this seems to be the right choice
        if(present(dt))then        
          split%nb_split_step = 11
          split%dim_split_V = 2
          SLL_ALLOCATE(split%split_step(17),ierr)
          split%split_begin_T = .false.
          !b1 a1 b2 a2 b3 a3 b3 a2 b2 a1 b1
          !b1
          split%split_step(1) = 0.086971698963920047813358_f64&
            +2._f64*dt**2*(1.98364114652831655458915e-6_f64)
          split%split_step(2) = dt**2*(1.98364114652831655458915e-6_f64)            
          !a1
          split%split_step(3) = 0.303629319055488881944104_f64
          !b2
          split%split_step(4) = 0.560744966588102145251453_f64&
            -2._f64*dt**2*(0.00553752115152236516667268_f64)
          split%split_step(5) = -dt**2*(0.00553752115152236516667268_f64)
          !a2
          split%split_step(6) = 0.303629319055488881944104_f64
          !b3
          split%split_step(7) = -0.1477166655520221930648117_f64&
            -2._f64*dt**2*(0.00284218110811634663914191_f64)
          split%split_step(8) = -dt**2*(0.00284218110811634663914191_f64)
          !a3
          split%split_step(9) = -0.2145172762219555277764167_f64
          !b3
          split%split_step(10) = split%split_step(7)
          split%split_step(11) = split%split_step(8)
          !a2
          split%split_step(12) = split%split_step(6)
          !b2
          split%split_step(13) = split%split_step(4)
          split%split_step(14) = split%split_step(5)
          !a1
          split%split_step(15) = split%split_step(3)
          !b1
          split%split_step(16) = split%split_step(1)
          split%split_step(17) = split%split_step(2)
        else
          print *,'#provide dt for use of case sll_p_order6vpot_vtv=',sll_p_order6vpot_vtv
          stop
        endif  


      case (sll_p_order6vpotnew3_vtv) ! Order 6 for Vlasov-Poisson VTV 2D with potential modif
        !we change also sign
        if(present(dt))then        
          split%nb_split_step = 13
          split%dim_split_V = 2
          SLL_ALLOCATE(split%split_step(20),ierr)
          split%split_begin_T = .false.
          !b1 a1 b2 a2 b3 a3 b4 a3 b3 a2 b2 a1 b1
          !b1
          split%split_step(1) = 0.0482332301753032567427580_f64&
            +2._f64*dt**2*(0.0002566567904012107264_f64)
          split%split_step(2) = dt**2*(0.0002566567904012107264_f64)            
          !a1
          split%split_step(3) = 0.2701015188126056215752542_f64
          !b2
          split%split_step(4) = 0.0482332301753032567427580_f64&
            +2._f64*dt**2*(0.0009439771580927593579_f64)
          split%split_step(5) = dt**2*(0.0009439771580927593579_f64)
          !a2
          split%split_step(6) = -0.108612186368692920020654_f64
          !b3
          split%split_step(7) = 0.2361392603742494444753990_f64&
            -2._f64*dt**2*(0.002494619878121813220_f64)
          split%split_step(8) = -dt**2*(0.002494619878121813220_f64)
          !a3
          split%split_step(9) = 0.3385106675560872984454001_f64
          !b4
          split%split_step(10) = 0.3347885585502880840781703_f64&
            -2._f64*dt**2*(0.002670269183371982607658111_f64)
          split%split_step(11) = -dt**2*(0.002670269183371982607658111_f64)
          !a3
          split%split_step(12) = split%split_step(9)
          !b3
          split%split_step(13) = split%split_step(7)
          split%split_step(14) = split%split_step(8)
          !a2
          split%split_step(15) = split%split_step(6)
          !b2
          split%split_step(16) = split%split_step(4)
          split%split_step(17) = split%split_step(5)
          !a1
          split%split_step(18) = split%split_step(3)
          !b1
          split%split_step(19) = split%split_step(1)
          split%split_step(20) = split%split_step(2)
        else
          print *,'#provide dt for use of case sll_p_order6vpot_vtv=',sll_p_order6vpot_vtv
          stop
        endif  







      case (sll_p_order6vpnew2_vtv) ! Order 6 for Vlasov-Poisson VTV (new2)
        if(present(dt))then        
          split%nb_split_step = 11
          SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
          split%split_begin_T = .false.
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
          print *,'#provide dt for use of case sll_p_order6vpnew2_vtv=',sll_p_order6vpnew2_vtv
          stop
        endif  
      case default
        print *,'#split_case not defined'
        stop       
    end select
  end subroutine initialize_time_splitting_coeff





end module sll_m_time_splitting_coeff

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
