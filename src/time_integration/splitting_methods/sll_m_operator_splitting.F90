!> @ingroup operator_splitting
!> @brief 
!> Base class of operator splitting library. 
!>
!> @details
!> The operator splitting module provides a generic implementation of composition algorithms 
!> of different order for two operators T and V for solving \f$ \frac{dU}{dt} = (T+V) U \f$. 
!> The solution on one time step can be written \f$ U(\Delta t) = \mathcal{S}_{T+V} U(0) \f$. 
!> The composition algorithm consists in successive solutions of the split equations 
!> \f$ \frac{dU}{dt} = T U \f$ and \f$ \frac{dU}{dt} = V U \f$. 
!> Alternating the two reduced solution operators \f$ \mathcal{S}_{T} \f$ 
!> and \f$\mathcal{S}_{V} \f$ with adequately chosen time increments yields arbitrary 
!> order in time for the full solution.
!>
!> The application of an operator splitting method to a concrete problem is done
!> by extending the operator_splitting splitting base class by a new type
!> containing on the one hand the data on which the operators act
!> and a specific implementation of the two operators
!>
module sll_m_operator_splitting
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    do_split_steps, &
    initialize_operator_splitting, &
    operator_splitting, &
    sll_lie_tv, &
    sll_lie_vt, &
    sll_order6_tvt, &
    sll_order6_vtv, &
    sll_strang_tvt, &
    sll_strang_vtv, &
    sll_triple_jump_tvt, &
    sll_triple_jump_vtv

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! global variables 
  sll_int32, parameter :: SLL_USER_DEFINED      = 0 !< user defined splitting
  sll_int32, parameter :: SLL_LIE_TV            = 1 !< Lie splitting T first (order 1)
  sll_int32, parameter :: SLL_LIE_VT            = 2 !< Lie splitting V first x (order 1)
  sll_int32, parameter :: SLL_STRANG_TVT        = 3 !< Strang splitting T first (order 2)
  sll_int32, parameter :: SLL_STRANG_VTV        = 4 !< Strang splitting V first (order 2)
  sll_int32, parameter :: SLL_TRIPLE_JUMP_TVT   = 5 !< Triple jump splitting T first (order 4)
  sll_int32, parameter :: SLL_TRIPLE_JUMP_VTV   = 6 !< Triple jump splitting V first (order 4)
  sll_int32, parameter :: SLL_ORDER6_TVT        = 7 !< Order 6 splitting T first (order 6)
  sll_int32, parameter :: SLL_ORDER6_VTV        = 8 !< Order 6 splitting V first (order 6)
  sll_int32, parameter :: SLL_ORDER6VP_TVT      = 9 !< Specific Vlasov-Poisson splitting T first (order 6)
  sll_int32, parameter :: SLL_ORDER6VP_VTV      = 10 !< Specific Vlasov-Poisson splitting V first (order 6)
  sll_int32, parameter :: SLL_ORDER6VPNEW_TVT   = 11 !< Specific Vlasov-Poisson splitting T first (order 6)
  sll_int32, parameter :: SLL_ORDER6VPNEW1_VTV  = 12 !< Specific Vlasov-Poisson splitting V first (order 6)
  sll_int32, parameter :: SLL_ORDER6VPNEW2_VTV  = 13 !< Specific Vlasov-Poisson splitting V first (order 6)
  sll_int32, parameter :: SLL_ORDER6VP2D_VTV    = 14 !< Specific Vlasov-Poisson splitting V first (order 6)
  sll_int32, parameter :: SLL_ORDER6VPOT_VTV    = 15 !< Specific Vlasov-Poisson splitting V first (order 6)

 !> operator splitting object
  type :: operator_splitting
     !> current time to be incremented 
     sll_real64 :: current_time = 0.0_f64  
     !> defines the splitting method to be chosen from those defined as global variables in 
     !> sll_m_operator_splitting module
     sll_int32 :: split_case  
     !> number of split steps in the method 
     sll_int32 :: nb_split_step
     !> array containing the coefficients of the each split step 
     sll_real64, dimension(:), pointer :: split_step
     !> Start with operatorT if true and with operatorV if false. 
     logical :: split_begin_T
     !> Used for specific Vlasov-Poisson splitting.
     sll_int32 :: dim_split_v
   contains 
     !> pointer on routine defined first operator
     procedure, pass(this) :: operatorT => operator
     !> pointer on routine defined second operator
     procedure, pass(this) :: operatorV => operator
  end type operator_splitting

contains
  !> @brief dummy implementation of an operator_splitting needed to provide the interface. 
  !> The class cannot be abstract because the splitting coefficients are initialised here 
  !> as they are generic for any type of operator.
  subroutine operator (this, dt)
    class(operator_splitting), intent(inout)  :: this !< operator_splitting object
    sll_real64, intent(in)  :: dt !< time increment on which operator is applied

    print*, 'This is a dummy implementation for providint the interface  &
    &    in sll_operator_splitting_base.F90.                             &
    &    The class operator_splitting needs to be extended               &
    &    providing the data to be evolved by the operator for an actual  &
    &    implementation.'

    return
    SLL_ASSERT(storage_size(this)>0)
    SLL_ASSERT(dt>0.0_f64)
  end subroutine operator

  !> @brief Returns a pointer to a heap-allocated operator_splitting object.
  !> @param[in]	split_case	splitting method to be chosen among
  !> SLL_USER_DEFINED : user provides coefficients <br>
  !> SLL_LIE_TV : first order Lie with T first <br>
  !> SLL_LIE_VT : first order Lie with V first <br>
  !> SLL_STRANG_TVT : second order Strang with T first <br>
  !> SLL_STRANG_VTV : second order Strang with V first <br>
  !> SLL_TRIPLE_JUMP_TVT : fourth order triple jump with T first <br>
  !> SLL_TRIPLE_JUMP_VTV : fourth order triple jump with V first <br>
  !> SLL_ORDER6_TVT : sixth order with T first <br>
  !> SLL_ORDER6_VTV : sixth order with V first <br>
  !> SLL_ORDER6VP_TVT : sixth order optimized for Vlasov-Poisson with T first <br>
  !> SLL_ORDER6VP_VTV : sixth order optimized for Vlasov-Poisson with V first <br>
  !> SLL_ORDER6VPnew_TVT <br>
  !> SLL_ORDER6VPnew1_VTV <br>
  !> SLL_ORDER6VPnew2_VTV <br>
  !> SLL_ORDER6VP2D_VTV <br>
  !> SLL_ORDER6VPOT_VTV
  !> @param[in]	split_step	array containing splitting coefficients
  !> @param[in]	nb_split_step	number of split steps in the method
  !> @param[in]	split_begin_T	logical True if T first false else
  !> @param[in]	dt	time step
  !> @return split a pointer to a heap-allocated operator_splitting object.
  function new_operator_splitting( &
       split_case, &
       split_step, &
       nb_split_step, &
       split_begin_T, &
       dt) &
       result(split)  
    class(operator_splitting), pointer :: split
    sll_int32, intent(in)  :: split_case
    sll_real64, dimension(:), intent(in), optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T
    sll_real64, intent(in), optional :: dt   
    sll_int32 :: ierr

    SLL_ALLOCATE(split,ierr)   
    call initialize_operator_splitting( &
         split, &
         split_case, &
         split_step, &
         nb_split_step, &
         split_begin_T, &
         dt)

  end function new_operator_splitting

  !> @brief Initialises data for operator splitting.
  !> @param[in]	split	operator splitting object to be initialised
  !> @param[in]	split_case	splitting method to be chosen among <br>
  !> SLL_USER_DEFINED : user provides coefficients <br>
  !> SLL_LIE_TV : first order Lie with T first <br>
  !> SLL_LIE_VT : first order Lie with V first <br>
  !> SLL_STRANG_TVT : second order Strang with T first <br>
  !> SLL_STRANG_VTV : second order Strang with V first <br>
  !> SLL_TRIPLE_JUMP_TVT : fourth order triple jump with T first <br>
  !> SLL_TRIPLE_JUMP_VTV : fourth order triple jump with V first <br>
  !> SLL_ORDER6_TVT : sixth order with T first <br>
  !> SLL_ORDER6_VTV : sixth order with V first <br>
  !> SLL_ORDER6VP_TVT : sixth order optimized for Vlasov-Poisson with T first <br>
  !> SLL_ORDER6VP_VTV : sixth order optimized for Vlasov-Poisson with V first <br>
  !> SLL_ORDER6VPnew_TVT <br>
  !> SLL_ORDER6VPnew1_VTV  <br>
  !> SLL_ORDER6VPnew2_VTV  <br>
  !> SLL_ORDER6VP2D_VTV  <br>
  !> SLL_ORDER6VPOT_VTV <br>
  !> @param[in]	split_step	array containing splitting coefficients
  !> @param[in]	nb_split_step	number of split steps in the method
  !> @param[in]	split_begin_T	logical True if T first false else
  !> @param[in]	dt	time step
  subroutine initialize_operator_splitting( &
       split, &
       split_case, &
       split_step, &
       nb_split_step, &
       split_begin_T, &
       dt)

    class(operator_splitting) :: split
    sll_int32, intent(in) :: split_case
    sll_real64, dimension(:), intent(in), optional :: split_step
    sll_int32, intent(in), optional :: nb_split_step
    logical, intent(in), optional :: split_begin_T
    sll_real64, intent(in), optional :: dt   
    sll_int32 :: ierr 

    split%dim_split_V = 1
    split%split_case = split_case

    if(split_case==sll_user_defined) then
       if( .not.( (present(split_step)) &
            .and.(present(nb_split_step)) &
            .and.(present(split_begin_T)) )) then
          print *,'#provide split_step, nb_split_step,split_begin_T'
          print *,'#in initialize_operator_splitting_coeff'
          stop
       endif
    else
       if(&
            (present(split_step)) &
            .or.(present(nb_split_step)) &
            .or.(present(split_begin_T)) &
            )then
          print *,'# do not provide split_step, nb_split_step, split_begin_T'
          print *,'#in initialize_operator_splitting_coeff'
          stop       
       endif
    endif

    select case (split_case)    
    case (sll_user_defined)
       if(size(split_step)<nb_split_step) then
          print *,'#bad size for split_step',size(split_step),nb_split_step
          print *,'#in initialize_operator_splitting_coeff'
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
       split%split_begin_T = .true.
       split%split_step(1) = 0.5_f64
       split%split_step(2) = 1._f64
       split%split_step(3) = split%split_step(1)
    case (SLL_STRANG_VTV) ! Strang splitting VTV
       split%nb_split_step = 3
       SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
       split%split_begin_T = .false.
       split%split_step(1) = 0.5_f64
       split%split_step(2) = 1._f64
       split%split_step(3) = split%split_step(1)
    case (SLL_TRIPLE_JUMP_TVT) ! triple jump TVT
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
    case (SLL_TRIPLE_JUMP_VTV) ! triple jump VTV
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
    case (SLL_ORDER6_TVT) ! Order 6 TVT
       split%nb_split_step = 23
       SLL_ALLOCATE(split%split_step(split%nb_split_step),ierr)
       split%split_begin_T = .true.
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
    case (SLL_ORDER6_VTV) ! Order 6 VTV
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
    case (SLL_ORDER6VP_TVT) ! Order 6 for Vlasov-Poisson TVT 
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
          print *,'#provide dt for use of case SLL_ORDER6VP_TVT=',sll_order6vp_tvt
          stop
       endif
    case (SLL_ORDER6VP_VTV) ! Order 6 for Vlasov-Poisson VTV 
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
          print *,'#provide dt for use of case SLL_ORDER6VP_VTV=',sll_order6vp_vtv
          stop
       endif
    case (SLL_ORDER6VPNEW_TVT) ! Order 6 for Vlasov-Poisson TVT (new)
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
          print *,'#provide dt for use of case SLL_ORDER6VPnew_TVT=',sll_order6vpnew_tvt
          stop
       endif

    case (SLL_ORDER6VPNEW1_VTV) ! Order 6 for Vlasov-Poisson VTV (new1)
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
          print *,'#provide dt for use of case SLL_ORDER6VPnew1_VTV=',sll_order6vpnew1_vtv
          stop
       endif
    case (SLL_ORDER6VP2D_VTV) ! Order 6 for Vlasov-Poisson VTV 2D
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
          print *,'#provide dt for use of case SLL_ORDER6VP2D_VTV=',sll_order6vp2d_vtv
          stop
       endif
    case (SLL_ORDER6VPOT_VTV) ! Order 6 for Vlasov-Poisson VTV 2D with potential modif
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
          print *,'#provide dt for use of case SLL_ORDER6VPOT_VTV=',sll_order6vpot_vtv
          stop
       endif
    case (SLL_ORDER6VPNEW2_VTV) ! Order 6 for Vlasov-Poisson VTV (new2)
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
          print *,'#provide dt for use of case SLL_ORDER6VPnew2_VTV=',sll_order6vpnew2_vtv
          stop
       endif
    case default
       print *,'#split_case not defined'
       stop       
    end select
  end subroutine initialize_operator_splitting

  !> @brief Apply the composition method for given number of times steps.
  !> @param[inout]	split	operator_splitting object
  !> @param[in]	dt	: time step
  !> @param[in]	number_time_steps	: number of time steps to be performed
  subroutine do_split_steps(split, dt, number_time_steps)
    class(operator_splitting)   :: split
    sll_real64, intent(in)  :: dt
    sll_int32, intent(in)   :: number_time_steps
    ! local variables
    sll_int32  :: i, istep

    ! Straight implementation of Lie splitting in both direction
    if (split%split_case == SLL_LIE_TV) then
       do i = 1, number_time_steps
          call split%operatorT(dt)
          call split%operatorV(dt)
          ! Increment current_time
          split%current_time = split%current_time + dt
       end do
    else if (split%split_case == SLL_LIE_VT) then
       do i = 1, number_time_steps
          call split%operatorV(dt)
          call split%operatorT(dt)
          ! Increment current_time
          split%current_time = split%current_time + dt
       end do
    else
       ! Combine first and last step of symmetric operators
       ! if more than one time step if performed
       !---------------------------------------------------
       if (split%split_begin_T) then
          ! Apply T first
          ! First split step
          istep = 1 
          call split%operatorT(split%split_step(istep)*dt)
          do i = 1, number_time_steps - 1
             ! Start with second split step as first already accounted
             ! for in last time step
             istep = 2
             do while (istep < split%nb_split_step - 2)
                call split%operatorV(split%split_step(istep)*dt)
                istep = istep + 1
                call split%operatorT(split%split_step(istep)*dt)
                istep = istep + 1
             end do
             ! Last two steps with T applied on twice the step as 
             ! it is combine with the first push of the next time step
             call split%operatorV(split%split_step(istep)*dt)
             istep = istep + 1
             call split%operatorT(2*split%split_step(istep)*dt)
             istep = istep + 1
          enddo
          ! Last time step of the sequence. Start with second split step
          ! as first already accounted for in last time step
          istep = 2
          do while (istep < split%nb_split_step)
             call split%operatorV(split%split_step(istep)*dt)
             istep = istep + 1
             call split%operatorT(split%split_step(istep)*dt)
             istep = istep + 1
          end do
          ! Increment current_time
          split%current_time = split%current_time + dt
       else     
          ! Apply V first
          ! First split step
          istep = 1 
          call split%operatorV(split%split_step(istep)*dt)
          do i = 1, number_time_steps - 1
             ! Start with second split step as first already accounted
             ! for in last time step
             istep = 2
             do while (istep < split%nb_split_step - 2)
                call split%operatorT(split%split_step(istep)*dt)
                istep = istep + 1
                call split%operatorV(split%split_step(istep)*dt)
                istep = istep + 1
             end do
             ! Last two steps with V applied on twice the step as 
             ! it is combine with the first push of the next time step
             call split%operatorT(split%split_step(istep)*dt)
             istep = istep + 1
             call split%operatorV(2*split%split_step(istep)*dt)
             istep = istep + 1
          enddo
          ! Last time step of the sequence. Start with second split step
          ! as first already accounted for in last time step
          istep = 2
          do while (istep < split%nb_split_step) 
             call split%operatorT(split%split_step(istep)*dt)
             istep = istep + 1
             call split%operatorV(split%split_step(istep)*dt)
             istep = istep + 1
          end do
          ! Increment current_time
          split%current_time = split%current_time + dt
       endif
    endif
  end subroutine do_split_steps


end module sll_m_operator_splitting
