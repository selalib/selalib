! Define coefficients for composition methods of symmetric and first order methods
!
!> @author Jakob Ameres (jakob.ameres@tum.de)
!

module sll_time_composition


#include "sll_working_precision.h"
#include "sll_memory.h"

  implicit none

  !> Symmetric Composition of Symmetric Methods
  type comp_coeff_sym_sym
    sll_int32 :: order
    sll_real64, dimension(:), allocatable :: gamma  !coefficients, already symmetric
    sll_int32 :: stages
    contains
    procedure, pass(this) :: init => time_composition_sym_sym
!     procedure, pass(this) :: coeff => comp_coeff_
  end type comp_coeff_sym_sym  
  
contains

!> returns available stages for composition of desired order
! function time_composition_sym_sym_stages( order)
! end function time_composition_sym_sym_stages


!> Symmetric Composition of Symmetric Methods
!> relaxing the number of stages can lead to better methods of high order
subroutine time_composition_sym_sym(this, order, stages) 
 sll_int32, intent(in) :: order, stages
 class(comp_coeff_sym_sym), intent(inout) :: this
 sll_int32 :: idx
 
 this%order=order
 this%stages=stages
 
 allocate(this%gamma(1:stages))
SELECT CASE( order)
   CASE(1)
    this%gamma=1.0_f64/stages
   
   !########################################
   CASE(2)
    SELECT CASE (stages)
    CASE(3)
     this%gamma(1)=1.3512071919596576340476878089715
     this%gamma(2)=-1.7024143839193152680953756179429
    END SELECT
   
   
   !########################################
   CASE(6)
    SELECT CASE (stages)
     CASE(7)
      this%gamma(1)=0.78451361047755726381949763_f64
      this%gamma(2)=0.23557321335935813368479318_f64
      this%gamma(3)=-1.17767998417887100694641568
      this%gamma(4)=-1.31518632068391121888424973
     CASE(9)
      this%gamma(1) = 0.39216144400731413927925056
      this%gamma(2)= 0.33259913678935943859974864
      this%gamma(3)= -0.70624617255763935980996482
      this%gamma(4)= 0.08221359629355080023149045
      this%gamma(5)= 0.79854399093482996339895035
    END SELECT
  
  !###################################################
  CASE (8)
     SELECT CASE (stages)
     CASE(15)
       this%gamma(1)= 0.74167036435061295344822780
       this%gamma(2)=-0.40910082580003159399730010
       this%gamma(3)= 0.19075471029623837995387626
       this%gamma(4)=-0.57386247111608226665638773
       this%gamma(5)= 0.29906418130365592384446354
       this%gamma(6)= 0.33462491824529818378495798
       this%gamma(7)= 0.31529309239676659663205666
       this%gamma(8)=-0.79688793935291635401978884
     CASE(17)
       this%gamma(1)= 0.13020248308889008087881763
       this%gamma(2)= 0.56116298177510838456196441
       this%gamma(3)=-0.38947496264484728640807860
       this%gamma(4)= 0.15884190655515560089621075
       this%gamma(5)=-0.39590389413323757733623154
       this%gamma(6)= 0.18453964097831570709183254
       this%gamma(7)= 0.25837438768632204729397911
       this%gamma(8)= 0.29501172360931029887096624
       this%gamma(9)=-0.60550853383003451169892108
     END SELECT
  !###################################################
  CASE(10)
     SELECT CASE (stages)
     CASE(35)
       this%gamma(1) = 0.07879572252168641926390768
       this%gamma(2) = 0.31309610341510852776481247
       this%gamma(3) = 0.02791838323507806610952027
       this%gamma(4) =-0.22959284159390709415121340
       this%gamma(5) = 0.13096206107716486317465686
       this%gamma(6) =-0.26973340565451071434460973
       this%gamma(7) = 0.07497334315589143566613711
       this%gamma(8) = 0.11199342399981020488957508
       this%gamma(9) = 0.36613344954622675119314812
       this%gamma(10)=-0.39910563013603589787862981
       this%gamma(11)= 0.10308739852747107731580277
       this%gamma(12)= 0.41143087395589023782070412
       this%gamma(13)=-0.00486636058313526176219566
       this%gamma(14)=-0.39203335370863990644808194
       this%gamma(15)= 0.05194250296244964703718290
       this%gamma(16)= 0.05066509075992449633587434
       this%gamma(17)= 0.04967437063972987905456880
       this%gamma(18)= 0.04931773575959453791768001
      END SELECT
END SELECT

!symmetry
do idx=1,stages/2  !fortran intrinsic ceiling
 this%gamma(stages-(idx-1))=this%gamma(idx)
end do

end subroutine time_composition_sym_sym


end module sll_time_composition

