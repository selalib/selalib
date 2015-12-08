!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

module sll_m_pif_fieldsolver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_i1, &
    sll_pi

  implicit none

  public :: &
    diag_dot_matrix_real64, &
    pif_fieldsolver

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    type :: pif_fieldsolver


    sll_int32 :: dimx !spatial dimensions
    sll_real64, allocatable, dimension(:) :: unitmode !normalized to domain fourier mode
    sll_int32, allocatable, dimension(:,:) :: allmodes
    sll_comp64, allocatable, dimension(:) :: rhs_one ! 1 
    contains
    procedure,  pass(this) :: init =>sll_pif_fieldsolver_init
    !sets the domain to be a d-dimensional cube with constant edge length
    procedure,  pass(this) :: set_box_len=>sll_pif_fieldsolver_set_box_len
    !sets user defined edge length in every direction
    procedure,  pass(this) :: set_box_lens=>sll_pif_fieldsolver_set_box_lens
    procedure, pass(this) :: problemsize=>sll_pif_fieldsolver_get_problemsize
    
    procedure,  pass(this) :: get_fourier_modes=>get_fourier_modes
    procedure,  pass(this) :: get_fourier_modes_chunk=>get_fourier_modes_chunk
    procedure,  pass(this) :: calc_fourier_modes=>calc_fourier_modes
    
    
    procedure,  pass(this) :: get_fourier_modes2=>get_fourier_modes2
    procedure, pass(this) :: calc_fourier_modes2=>calc_fourier_modes2
    procedure,  pass(this) :: get_fourier_modes2_chunk=>get_fourier_modes2_chunk
    procedure, pass(this) :: calc_fourier_modes2_chunk=>calc_fourier_modes2_chunk
    
    procedure, pass(this) :: solve_poisson=>sll_pif_fieldsolver_solve_poisson
    procedure, pass(this) :: solve_mass=>sll_pif_fieldsolver_solve_mass
    procedure, pass(this) :: solve_quasineutral=>sll_pif_fieldsolver_solve_quasineutral
    procedure, pass(this) ::  solve_qn_rho_wo_zonalflow=>sll_pif_fieldsolver_solve_qn_rho_wo_zonalflow
    procedure, pass(this) ::  eval_gradient=>sll_pif_fieldsolver_eval_gradient
    procedure, pass(this) ::  eval_solution=>sll_pif_fieldsolver_eval_solution
    procedure, pass(this) :: get_rhs_particle=>get_fourier_modes
    procedure, pass(this) :: visu_info=>visu_info_sll_pif_fieldsolver
     
     procedure, pass(this) :: l2norm=>l2norm_sll_pif_fieldsolver
 end type pif_fieldsolver

contains
pure function sll_pif_fieldsolver_get_problemsize(this) result(sz)
 class(pif_fieldsolver), intent(in) :: this
 sll_int32 :: sz
 if (allocated(this%allmodes)) then
 sz=size(this%allmodes,2)
 else
 sz=0
 endif
end function sll_pif_fieldsolver_get_problemsize

!>Returns the l2norm for a coefficient vector of a solution
function l2norm_sll_pif_fieldsolver(this, solution) result(l2norm)
 class(pif_fieldsolver), intent(in) :: this
 sll_comp64, dimension(:), intent(in) :: solution
!  sll_int32 :: idx
 sll_real64 :: l2norm
 
 SLL_ASSERT(size(solution)==this%problemsize())
 
 l2norm=real(sqrt(dot_product(solution, solution)),f64)
 
end function l2norm_sll_pif_fieldsolver



subroutine sll_pif_fieldsolver_init(this,maxmode)
 class(pif_fieldsolver), intent(inout) :: this
 sll_int32, intent(in) :: maxmode
 sll_int32, allocatable, dimension(:) :: maxmodes,minmodes
 sll_int32 :: ierr,idx
  
 SLL_ALLOCATE(maxmodes(1:this%dimx),ierr)
 SLL_ALLOCATE(minmodes(1:this%dimx),ierr)
 
 maxmodes=maxmode
 minmodes=-maxmode
 minmodes(1)=0
  SLL_ALLOCATE(this%allmodes(1:this%dimx,1:product(maxmodes-minmodes+1)),ierr)
  this%allmodes=generate_exponents(minmodes,maxmodes)
  
 !Define the one  
 SLL_ALLOCATE(this%rhs_one(1:this%problemsize()),ierr)
 do idx=1,size(this%allmodes,2)
   if (sum(abs(this%allmodes(:,idx)))==0) then
      this%rhs_one(idx)=cmplx(product((2*sll_pi/this%unitmode)),0.,f64)
   endif
 end do
end subroutine  sll_pif_fieldsolver_init


subroutine visu_info_sll_pif_fieldsolver(this)
 class(pif_fieldsolver), intent(inout) :: this

print *,"Spatial Dimensions (x)", this%dimx
print *,"Number of Fourier modes: ", this%problemsize()
print *,"Domain length", 1.0/this%unitmode*sll_pi*2.0

end subroutine


subroutine sll_pif_fieldsolver_set_box_len(this, length)
 class(pif_fieldsolver), intent(inout) :: this
 sll_real64, intent(in) :: length
 sll_int32 :: ierr 

 if (.not. allocated(this%unitmode)) then
 SLL_ALLOCATE(this%unitmode(this%dimx),ierr)
 endif
 
 this%unitmode=2*sll_pi/length
end subroutine sll_pif_fieldsolver_set_box_len
 
 
 
subroutine sll_pif_fieldsolver_set_box_lens(this, lengths)
 class(pif_fieldsolver), intent(inout) :: this
 sll_real64, intent(in),dimension(:) :: lengths
 sll_int32 :: ierr 
 
 SLL_ASSERT(size(lengths)==this%dimx)
 
 if (.not. allocated(this%unitmode)) then
 SLL_ALLOCATE(this%unitmode(this%dimx),ierr)
 endif
 
 this%unitmode(:)=2*sll_pi/lengths(:)
end subroutine sll_pif_fieldsolver_set_box_lens




! sll_comp64 function get_fourier_mode( particle, mode)
!  sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,v,weight)
!  sll_real64, dimension(dimx), intent(in) :: mode !normalized to domain fourier mode
!  !warning the dot_product(x,y)=conjg(x)*y
!  get_fourier_mode=dot_product(particle(dimx+dimv+1,:),exp(-sll_i1*matmul(mode,particle(1:dimx,:)))) 
! end function get_fourier_mode
 

! sll_pif_fieldsolver_get_rhs
 

 function get_fourier_modes2_chunk(this, particle,chunksize) result(fouriermodes)
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
   sll_int32, intent(in) :: chunksize
 sll_comp64, dimension(:), allocatable :: fouriermodes
 sll_int32 :: ierr
 
 SLL_ALLOCATE(fouriermodes(1:this%problemsize()),ierr)
 call this%calc_fourier_modes2_chunk(particle,fouriermodes,chunksize)
end function get_fourier_modes2_chunk
!  
 
 
 subroutine calc_fourier_modes2_chunk(this, particle,fouriermodes,chunksize)
  class(pif_fieldsolver), intent(in) :: this
  sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
  sll_comp64, dimension(:),intent(out) :: fouriermodes
  sll_int32, intent(in) :: chunksize
  sll_int32 :: num, chunk
  sll_comp64,  dimension(size(fouriermodes)) :: fmodechunk

  
  num=size(particle,2)
  fouriermodes=(0.0_f64,0.0_f64)
do chunk=1,ceiling(real(num/chunksize,8))
  call this%calc_fourier_modes2&
      (particle(:, (chunk-1)*chunksize+1:min(chunk*chunksize,num)), fmodechunk)
  fouriermodes=fouriermodes+fmodechunk;
enddo
 end subroutine calc_fourier_modes2_chunk
 
 
 
 
 
 !please apply this chunked
subroutine calc_fourier_modes2(this, particle,fouriermodes)
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
 !sll_real64, dimension(size(particle,2)) :: weight
 sll_comp64, dimension(this%dimx,size(particle,2)) :: unitmodes
 sll_comp64, dimension(:),intent(out) :: fouriermodes
 sll_int32 :: idx!,ierr
 !Calculate unit fourier modes first
 unitmodes=exp(cmplx(0.0,-diag_dot_matrix_real64(this%unitmode, particle(1:this%dimx,:)),f64));
 !Extract weights from particle array
 !weight=particle(this%dimx+1,:)
 
 do idx=1, this%problemsize()
  fouriermodes(idx) =  &
  sum(product(array_exponent_comp64(unitmodes,this%allmodes(:,idx)),1)*cmplx(particle(this%dimx+1,:),0.0,f64))
 end do
 
end subroutine calc_fourier_modes2

 
!please apply this chuncked
function get_fourier_modes2(this, particle) result(fouriermodes)
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
 
 sll_comp64, dimension(:), allocatable :: fouriermodes
 sll_int32 :: ierr
 
 SLL_ALLOCATE(fouriermodes(1:this%problemsize()),ierr)
 call this%calc_fourier_modes2(particle,fouriermodes)
end function get_fourier_modes2
!  

!>Get phi from the right hand side
!> - Phi''= pho
function sll_pif_fieldsolver_solve_poisson(this, rhs) result(solution)
 class(pif_fieldsolver), intent(in) :: this
 sll_comp64, dimension(:), intent(in) :: rhs
 sll_comp64, dimension(size(rhs)) :: solution
 sll_int32 :: idx
 
 SLL_ASSERT(size(rhs)==this%problemsize())
 do idx=1,size(rhs)
 !intermit constant mode
  if (sum(abs(this%allmodes(:,idx)))/=0) then
   solution(idx)=rhs(idx)                               &
     / cmplx(sum((this%allmodes(:,idx)*this%unitmode(:))**2)  &
     *       product(this%unitmode/sll_pi/2.0_f64),0.0,f64)
  else
       solution(idx)=(0.0_f64,0.0_f64)
  endif
 end do
 
  do idx=1,this%problemsize()
  if ( .not. this%allmodes(1,idx)==0) then
   solution(idx)=solution(idx)*2
  endif
 end do
end function sll_pif_fieldsolver_solve_poisson




function sll_pif_fieldsolver_solve_qn_rho_wo_zonalflow(this, rhs) result(solution)
 class(pif_fieldsolver), intent(in) :: this
 sll_comp64, dimension(:), intent(in) :: rhs
 sll_comp64, dimension(size(rhs)) :: solution
 sll_int32 :: idx
 sll_int32 :: zonaldim=3
 
 SLL_ASSERT(size(rhs)==this%problemsize())
 do idx=1,size(rhs)
 !intermit constant mode
  if (sum(abs(this%allmodes(:,idx)))/=0) then
   solution(idx)=rhs(idx)*cmplx(product(this%unitmode/sll_pi/2),0.0,f64)
  else
       solution(idx)=(0.0_f64,0.0_f64)
  endif
  
  if (this%dimx==3) then
  !remove Zonal flow in zonaldim
  if (this%allmodes(zonaldim,idx)==0) then
    solution(idx)=(0.0_f64,0.0_f64)
  endif
  endif
  
 end do
 
  do idx=1,this%problemsize()
  if (this%allmodes(1,idx)/=0) then
   solution(idx)=solution(idx)*2
  endif
 end do
 
end function sll_pif_fieldsolver_solve_qn_rho_wo_zonalflow






!Get rho from the right hand side
function sll_pif_fieldsolver_solve_mass(this, rhs) result(solution)
 class(pif_fieldsolver), intent(in) :: this
 sll_comp64, dimension(:), intent(in) :: rhs
 sll_comp64, dimension(size(rhs)) :: solution
 sll_int32 :: idx
 
 SLL_ASSERT(size(rhs)==this%problemsize())
 
 solution=rhs*cmplx(product(this%unitmode/sll_pi/2),0.0_f64,f64)
 do idx=1,this%problemsize()
  if ( .not. this%allmodes(1,idx)==0) then
   solution(idx)=solution(idx)*2
  endif
 end do
end function sll_pif_fieldsolver_solve_mass

!Get rho from the right hand side
function sll_pif_fieldsolver_solve_quasineutral(this, rhs) result(solution)
 class(pif_fieldsolver), intent(in) :: this
 sll_comp64, dimension(:), intent(in) :: rhs
 sll_comp64, dimension(size(rhs)) :: solution
 sll_int32 :: idx
 
 SLL_ASSERT(size(rhs)==this%problemsize())
 
 solution=rhs*cmplx(product(this%unitmode/sll_pi/2),0.0,f64)
 do idx=1,this%problemsize()
  if ( .not. this%allmodes(1,idx)==0) then
   solution(idx)=solution(idx)*2
  endif
  
  !Remove constant mode, could also be a dimensional average
  if (all(this%allmodes(:,idx)==0)) then
    solution(idx)=(0.0_f64,0.0_f64)
  endif
 end do
end function sll_pif_fieldsolver_solve_quasineutral

function sll_pif_fieldsolver_eval_gradient(this, pos,fouriermodes) result(gradient)
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: pos !pos vector (x), no weights
 sll_comp64, dimension(:), intent(in) :: fouriermodes
 sll_real64, dimension(size(pos,1),size(pos,2)) ::  gradient
 !sll_comp64 :: partmode
 sll_comp64, dimension(size(pos,2)) :: partmode
 sll_int32 :: idx, jdx
 gradient=0.0_f64
 
!   do jdx=1,size(pos,2)
!   
!    do idx=1,this%problemsize()
!        partmode=sll_i1*exp(sll_i1*(dot_product(this%allmodes(:,idx)*this%unitmode(:),pos(1:this%dimx,jdx))  ))
!  
!        gradient(:,jdx)=gradient(:,jdx)+(real(partmode)*real(fouriermodes(idx))-aimag(partmode)*aimag(fouriermodes(idx)))*&
!                               (this%allmodes(:,idx)*this%unitmode(:))
!    end do
!  end do
 
 do idx=1,this%problemsize()
 partmode=sll_i1*exp(cmplx(0.0_f64,matmul(this%allmodes(:,idx)*this%unitmode(:), pos(1:this%dimx,:)),f64))
 do jdx=1,size(partmode)
      gradient(:,jdx)=gradient(:,jdx)+(real(partmode(jdx))*real(fouriermodes(idx))-aimag(partmode(jdx))*aimag(fouriermodes(idx)))*&
                             (this%allmodes(:,idx)*this%unitmode(:))
  end do
 end do
end function sll_pif_fieldsolver_eval_gradient


function sll_pif_fieldsolver_eval_solution(this, pos,fouriermodes) result(fun)
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: pos !pos vector (x), no weights
 sll_comp64, dimension(:), intent(in) :: fouriermodes
 sll_real64, dimension(size(pos,2)) ::  fun
 sll_comp64, dimension(size(pos,2)) :: partmode
 sll_int32 :: idx!, jdx
 fun=0.0_f64
 do idx=1,this%problemsize()
   partmode=exp(cmplx(0.0_f64,matmul(this%allmodes(:,idx)*this%unitmode(:), pos(1:this%dimx,:)),f64))
   fun(:)=fun(:)+(real(partmode(:))*real(fouriermodes(idx))-aimag(partmode(:))*aimag(fouriermodes(idx)))
 end do
end function sll_pif_fieldsolver_eval_solution



function get_fourier_modes_chunk(this, particle, chunksize) result(fouriermodes)
  class(pif_fieldsolver), intent(in) :: this
  sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
  sll_comp64, dimension(:), allocatable  :: fouriermodes
  sll_int32, intent(in) :: chunksize
  sll_int32 :: num, ierr, chunk
  sll_comp64,  dimension(this%problemsize()) :: fmodechunk
  
   SLL_ALLOCATE(fouriermodes(1:this%problemsize()),ierr)

  num=size(particle,2)
  fouriermodes=(0.0_f64,0.0_f64)
do chunk=1,ceiling(real(num/chunksize,8))
  call this%calc_fourier_modes&
      (particle(:, (chunk-1)*chunksize+1:min(chunk*chunksize,num)), fmodechunk)
  fouriermodes=fouriermodes+fmodechunk;
enddo
end function get_fourier_modes_chunk



subroutine calc_fourier_modes(this, particle, fouriermodes) 
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
 sll_comp64, dimension(:), intent(out) :: fouriermodes
 sll_int32 :: idx

do idx=1, this%problemsize()
  fouriermodes(idx)=sum(exp(-cmplx(0.0_f64,matmul(this%allmodes(:,idx)*this%unitmode,particle(1:this%dimx,:)),f64))*cmplx(particle(this%dimx+1,:),0.0,f64))
 end do
 
end subroutine calc_fourier_modes

function kahan_sum_comp64(summands) result(sum)
 sll_comp64, dimension(:), intent(in) :: summands
 sll_comp64 :: sum, c,t,y
 sll_int32 :: idx
 sum=(0.0_f64,0.0_f64)
 c=(0.0_f64,0.0_f64)
 t=(0.0_f64,0.0_f64)
    do idx = 1,size(summands)
        y = summands(idx) - c  
        t = sum + y         
        c = (t - sum) - y 
        sum = t 
  end do
end function


function kahan_sum_real64(summands) result(sum)
 sll_real64, dimension(:), intent(in) :: summands
 sll_real64 :: sum, c,t,y
 sll_int32 :: idx
 sum=0.0_f64
 c=0.0_f64
 t=0.0_f64
    do idx = 1,size(summands)
        y = summands(idx) - c  
        t = sum + y         
        c = (t - sum) - y 
        sum = t 
  end do
end function


function get_fourier_modes(this, particle) result(fouriermodes)
 class(pif_fieldsolver), intent(in) :: this
 sll_real64, dimension(:,:), intent(in)  :: particle !particle vector (x,weight)
 sll_comp64, dimension(:), allocatable :: fouriermodes
 sll_int32 :: ierr
 
 SLL_ALLOCATE(fouriermodes(1:this%problemsize()),ierr)
 call this%calc_fourier_modes(particle,fouriermodes)
end function get_fourier_modes
 
 
recursive function generate_exponents(min_exponents, max_exponents) result(list)
 sll_int32, dimension(:),intent(in) :: min_exponents,max_exponents
 sll_int32, dimension(:,:), allocatable :: list
 sll_int32, dimension(:,:), allocatable :: sublist
 sll_int32 :: dim, idx, sz, dz, num,subsz
sll_int32 :: ierr
 dim=size(min_exponents);
 
 !determine problemsize
 sz=product(max_exponents-min_exponents+1)
 dz=max_exponents(1)-min_exponents(1)+1
 subsz=sz/dz
 
 !allocate memory
 SLL_ALLOCATE(list(dim,1:sz),ierr)
 SLL_ALLOCATE(sublist(dim-1, subsz),ierr)
 
 if (dim>1) then
 sublist=generate_exponents(min_exponents(2:dim),max_exponents(2:dim))

 num=min_exponents(1);
 do idx=1,dz 
 list(1,(1+(idx-1)*subsz):idx*subsz)=num
 list(2:dim,(1+(idx-1)*subsz):idx*subsz)=sublist
 num=num+1
 end do
 else
 num=min_exponents(1);
 do idx=1,dz
 list(1,idx)=num
 num=num+1
 end do
 endif
end function generate_exponents
 
 
! function powers_comp64( basis, min_exponent, max_exponent)
!  sll_comp64, dimension(:
! 
!  fouriermode(idx)=fouriermode(idx-1)
! 
! end
 
 
function array_exponent_comp64(basis, exponent)
  sll_comp64, dimension(:,:),intent(in) :: basis
  sll_int32, dimension(:),intent(in) :: exponent
  sll_comp64, dimension(size(basis,1),size(basis,2)) :: array_exponent_comp64
  sll_int32 :: idx
  
  do idx=1,size(basis,2)
  !If exponent is negative take complex conjugate, compiler should know...
  array_exponent_comp64(:,idx)=basis(:,idx)**exponent(:)
  end do
end function array_exponent_comp64
 
 
! function weighted_array_power(array, weight, powers)
!  sll_real64, dimension(:) :: array
! 
! end function
 
 
function  diag_dot_matrix_real64( diagonal,matrix)
 sll_real64, dimension(:,:), intent(in) :: matrix !full matrix
 sll_real64, dimension(:), intent(in) :: diagonal !Entries of a diagonal matrix
  sll_real64 , dimension(size(matrix,1),size(matrix,2)) :: diag_dot_matrix_real64
  sll_int32 :: idx,sz
 
 SLL_ASSERT(size(matrix,1)==size(diagonal))
 sz=size(matrix,2)
 do idx=1,sz
 diag_dot_matrix_real64(:,idx)=diagonal(:)*matrix(:,idx)
 end do
end function diag_dot_matrix_real64
 
 
end module sll_m_pif_fieldsolver
