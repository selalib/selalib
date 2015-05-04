
module contrib_rho_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use cubic_non_uniform_splines
  !use utils
  implicit none
contains  

subroutine compute_rho_mapped_mesh&
  (rho,f,integration_points,rho_case,nc_eta1,nc_eta2,geom,jac_array,spl_per_x1)
  sll_int :: rho_case
  sll_real64,dimension(:,:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_val
  sll_real64,dimension(:), pointer :: rho
  sll_real64,dimension(:), pointer :: node_positions_x1
  sll_real64,dimension(:), pointer :: new_node_positions
  sll_real64,dimension(:), pointer :: buf_1d
  sll_real64,dimension(:,:), pointer :: jac_array
  sll_real64,dimension(:,:), pointer :: f
  sll_real64,intent(in) :: geom(2,2)
  sll_int,intent(in) :: nc_eta1,nc_eta2
  sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max
  sll_int :: i1,i2
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1
  sll_int :: err

  eta1_min = geom(1,1)
  eta1_max = geom(2,1)
  eta2_min = geom(1,2)
  eta2_max = geom(2,2)
  
  SLL_ALLOCATE(node_positions_x1(nc_eta1+1),err)
  SLL_ALLOCATE(new_node_positions(nc_eta1),err)
  SLL_ALLOCATE(buf_1d(nc_eta1+1),err)
  SLL_ALLOCATE(integration_points_val(2,nc_eta1+1),err)
  do i1=1,nc_eta1+1
    !node_positions_x1(i1) = eta1_min+(real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
    node_positions_x1(i1) = eta1_min+(real(i1,f64)-1._f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
  enddo
    do i2=1,nc_eta2
      do i1 = 1,nc_eta1
        new_node_positions(i1) = integration_points(1,i1,i2)
        if((new_node_positions(i1)>eta1_max).or.(new_node_positions(i1)<eta1_min) )then
          print *,'problem of new_node_position:',new_node_positions(i1),eta1_min,eta1_max
          stop
        endif
        if(new_node_positions(i1)<node_positions_x1(1))then
          new_node_positions(i1)=new_node_positions(i1)+eta1_max-eta1_min
        endif      
        !buf_1d(i1) = sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)!-f_equil(i1,i2)
        buf_1d(i1) = f(i1,i2)/jac_array(i1,i2)
        
      enddo
      buf_1d(nc_eta1+1) = buf_1d(1)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      call interpolate_array_value_nonunif( new_node_positions, buf_1d(1:nc_eta1), nc_eta1, spl_per_x1)
      do i1 = 1,nc_eta1
        integration_points(3,i1,i2) =  buf_1d(i1)
      enddo
    enddo
   
   

  !Compute rho0 and E0
    do i1 = 1, nc_eta1
      do i2=1,nc_eta2
        integration_points_val(1,i2) = integration_points(2,i1,i2)
        integration_points_val(2,i2) = integration_points(3,i1,i2)
      enddo
      rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2,rho_case)
!      if(rho_case==1)then
!        rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2)
!      endif
!      if(rho_case==2)then
!        rho(i1)=compute_non_unif_integral_spline(integration_points_val,nc_eta2)
!      endif
!      if(rho_case==3)then
!        rho(i1)=compute_non_unif_integral_gaussian(integration_points_val,nc_eta2)
!      endif      
!      if(rho_case==4)then
!        rho(i1)=compute_non_unif_integral_gaussian_sym(integration_points_val,nc_eta2)
!      endif      
!      !if(test_case==4)then      
!      !  rho(i1) = rho(i1)+1._f64
!      !endif  
   enddo
   
   rho(nc_eta1+1)=rho(1)
 
     
  SLL_DEALLOCATE(node_positions_x1,err)
  SLL_DEALLOCATE(new_node_positions,err)
  SLL_DEALLOCATE(buf_1d,err)
  SLL_DEALLOCATE(integration_points_val,err)



end subroutine compute_rho_mapped_mesh


subroutine compute_rho_mapped_mesh2&
  (rho,f,integration_points,rho_case,nc_eta1,nc_eta2,geom,jac_array,spl_per_x1)
  sll_int :: rho_case
  sll_real64,dimension(:,:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_val
  sll_real64,dimension(:), pointer :: rho
  sll_real64,dimension(:), pointer :: node_positions_x1
  sll_real64,dimension(:), pointer :: new_node_positions
  sll_real64,dimension(:), pointer :: buf_1d
  sll_real64,dimension(:,:), pointer :: jac_array
  sll_real64,dimension(:,:), pointer :: f
  sll_real64,intent(in) :: geom(2,2)
  sll_int,intent(in) :: nc_eta1,nc_eta2
  sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max
  sll_int :: i1,i2
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1
  sll_int :: err

  eta1_min = geom(1,1)
  eta1_max = geom(2,1)
  eta2_min = geom(1,2)
  eta2_max = geom(2,2)
  
  SLL_ALLOCATE(node_positions_x1(nc_eta1+1),err)
  SLL_ALLOCATE(new_node_positions(nc_eta1),err)
  SLL_ALLOCATE(buf_1d(nc_eta1+1),err)
  SLL_ALLOCATE(integration_points_val(2,nc_eta1+1),err)
  do i1=1,nc_eta1+1
    node_positions_x1(i1) = eta1_min+(real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
    !node_positions_x1(i1) = eta1_min+(real(i1,f64)-1._f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
  enddo
    do i2=1,nc_eta2
      do i1 = 1,nc_eta1
        new_node_positions(i1) = integration_points(1,i1,i2)
        if((new_node_positions(i1)>eta1_max).or.(new_node_positions(i1)<eta1_min) )then
          print *,'problem of new_node_position:',new_node_positions(i1),eta1_min,eta1_max
          stop
        endif
        if(new_node_positions(i1)<node_positions_x1(1))then
          new_node_positions(i1)=new_node_positions(i1)+eta1_max-eta1_min
        endif      
        !buf_1d(i1) = sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)!-f_equil(i1,i2)
        buf_1d(i1) = f(i1,i2)/jac_array(i1,i2)
        
      enddo
      buf_1d(nc_eta1+1) = buf_1d(1)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      call interpolate_array_value_nonunif( new_node_positions, buf_1d(1:nc_eta1), nc_eta1, spl_per_x1)
      do i1 = 1,nc_eta1
        integration_points(3,i1,i2) =  buf_1d(i1)
      enddo
    enddo
   
   

  !Compute rho0 and E0
    do i1 = 1, nc_eta1
      do i2=1,nc_eta2
        integration_points_val(1,i2) = integration_points(2,i1,i2)
        integration_points_val(2,i2) = integration_points(3,i1,i2)
      enddo
      rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2,rho_case)
!      if(rho_case==1)then
!        rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2)
!      endif
!      if(rho_case==2)then
!        rho(i1)=compute_non_unif_integral_spline(integration_points_val,nc_eta2)
!      endif
!      if(rho_case==3)then
!        rho(i1)=compute_non_unif_integral_gaussian(integration_points_val,nc_eta2)
!      endif      
!      if(rho_case==4)then
!        rho(i1)=compute_non_unif_integral_gaussian_sym(integration_points_val,nc_eta2)
!      endif      
!      !if(test_case==4)then      
!      !  rho(i1) = rho(i1)+1._f64
!      !endif  
   enddo
   
   rho(nc_eta1+1)=rho(1)
 
     
  SLL_DEALLOCATE(node_positions_x1,err)
  SLL_DEALLOCATE(new_node_positions,err)
  SLL_DEALLOCATE(buf_1d,err)
  SLL_DEALLOCATE(integration_points_val,err)



end subroutine compute_rho_mapped_mesh2



function compute_non_unif_integral(integration_points,N_points,rho_case)
  sll_real64 :: compute_non_unif_integral
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_int,intent(in) :: N_points,rho_case
  if(rho_case==1)then
    compute_non_unif_integral= compute_non_unif_integral_trapezoid(integration_points,N_points)
  endif
  if(rho_case==2)then
    compute_non_unif_integral=compute_non_unif_integral_spline(integration_points,N_points)
  endif
  if(rho_case==3)then
    compute_non_unif_integral=compute_non_unif_integral_gaussian(integration_points,N_points)
  endif      
  if(rho_case==4)then
    compute_non_unif_integral=compute_non_unif_integral_gaussian_sym(integration_points,N_points)
  endif        
  if(rho_case==5)then
    compute_non_unif_integral=compute_non_unif_integral_spline_per(integration_points,N_points)
  endif
  
end  function compute_non_unif_integral



function compute_non_unif_integral_trapezoid(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_trapezoid
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_int,intent(in) :: N_points
  sll_int :: i
  sll_real64 :: tmp,x1,x2,fval1,fval2
  compute_non_unif_integral_trapezoid = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    if(x2<x1)then
      print *,i,'bad integration points x1=',x1,'x2=',x2
      stop
    endif
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    tmp = 0.5_f64*(fval1+fval2)*(x2-x1)
    compute_non_unif_integral_trapezoid=compute_non_unif_integral_trapezoid+tmp
  enddo
  
  
end  function compute_non_unif_integral_trapezoid


function compute_non_unif_integral_spline_old(integration_points,N_points,Nb)
  sll_real64 :: compute_non_unif_integral_spline_old
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_fine
  sll_int,intent(in) :: N_points,Nb
  sll_int :: i,N_points_fine,ierr,j
  sll_real64 :: x1,x2
  type(cubic_nonunif_spline_1D), pointer :: spl
  compute_non_unif_integral_spline_old = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  N_points_fine = (N_points-1)*Nb+1
  spl =>  new_cubic_nonunif_spline_1D( N_points-1, SLL_HERMITE)
  SLL_ALLOCATE(integration_points_fine(2,N_points_fine),ierr)
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    if(x2<x1)then
      print *,i,'(spl) bad integration points x1=',x1,'x2=',x2
      stop
    endif
    do j=1,Nb
      integration_points_fine(1,(i-1)*Nb+j)=x1+real(j-1,f64)/real(Nb,f64)*(x2-x1)
    enddo
  enddo  
  integration_points_fine(1,N_points_fine)=integration_points(1,N_points)
  !stop  
  call compute_spline_nonunif( integration_points(2,1:N_points), spl, integration_points(1,1:N_points),0._f64,0._f64)
  call interpolate_array_value_nonunif( integration_points_fine(1,1:N_points_fine), &
  &integration_points_fine(2,1:N_points_fine),N_points_fine-1, spl)
  call delete_cubic_nonunif_spline_1D( spl, ierr)

  do i=1,N_points_fine-1
     x1 = integration_points_fine(1,i)
     x2 = integration_points_fine(1,i+1)
     if(x2<x1)then
       print *,i,'(spl2) bad integration points x1=',x1,'x2=',x2
       stop
     endif
     !print *,i,integration_points_fine(1,i)
  enddo
  
  
  compute_non_unif_integral_spline_old = compute_non_unif_integral_trapezoid(integration_points_fine,N_points_fine)
  SLL_DEALLOCATE_ARRAY(integration_points_fine,ierr)
  
end  function compute_non_unif_integral_spline_old

function compute_non_unif_integral_spline(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_spline
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_middle
  sll_int,intent(in) :: N_points
  sll_int :: i,ierr
  sll_real64 :: tmp,x1,x2,fval1,fval2,fvalm
  type(cubic_nonunif_spline_1D), pointer :: spl
  compute_non_unif_integral_spline = 0._f64

  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  spl =>  new_cubic_nonunif_spline_1D( N_points-1, SLL_HERMITE)
  SLL_ALLOCATE(integration_points_middle(2,N_points-1),ierr)
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    integration_points_middle(1,i)=0.5_f64*(x1+x2)
  enddo  
  call compute_spline_nonunif( integration_points(2,1:N_points), spl, integration_points(1,1:N_points),0._f64,0._f64)
  call interpolate_array_value_nonunif( integration_points_middle(1,1:N_points-1), &
  &integration_points_middle(2,1:N_points-1),N_points-1, spl)
  call delete_cubic_nonunif_spline_1D( spl, ierr)
  
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    fvalm = integration_points_middle(2,i)
    tmp = (fval1+4._f64*fvalm+fval2)*(x2-x1)/6._f64
    compute_non_unif_integral_spline=compute_non_unif_integral_spline+tmp
  enddo
  
  SLL_DEALLOCATE_ARRAY(integration_points_middle,ierr)
  
end  function compute_non_unif_integral_spline


function compute_non_unif_integral_spline_per(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_spline_per
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_middle
  sll_int,intent(in) :: N_points
  sll_int :: i,ierr
  sll_real64 :: tmp,x1,x2,fval1,fval2,fvalm
  type(cubic_nonunif_spline_1D), pointer :: spl
  compute_non_unif_integral_spline_per = 0._f64

  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  spl =>  new_cubic_nonunif_spline_1D( N_points-1, SLL_PERIODIC)
  SLL_ALLOCATE(integration_points_middle(2,N_points-1),ierr)
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    integration_points_middle(1,i)=0.5_f64*(x1+x2)
  enddo  
  call compute_spline_nonunif( integration_points(2,1:N_points), spl, integration_points(1,1:N_points))
  call interpolate_array_value_nonunif( integration_points_middle(1,1:N_points-1), &
  &integration_points_middle(2,1:N_points-1),N_points-1, spl)
  call delete_cubic_nonunif_spline_1D( spl, ierr)
  
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    fvalm = integration_points_middle(2,i)
    tmp = (fval1+4._f64*fvalm+fval2)*(x2-x1)/6._f64
    compute_non_unif_integral_spline_per=compute_non_unif_integral_spline_per+tmp
  enddo
  
  SLL_DEALLOCATE_ARRAY(integration_points_middle,ierr)
  
end  function compute_non_unif_integral_spline_per



function compute_non_unif_integral_gaussian(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_gaussian
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_new
  sll_int,intent(in) :: N_points
  sll_int::i,ierr
  compute_non_unif_integral_gaussian = 0.5_f64*compute_non_unif_integral_gaussian_sym(integration_points,N_points)
  
  SLL_ALLOCATE(integration_points_new(2,N_points),ierr)
  
  
  integration_points_new(1,1:N_points)=integration_points(1,1:N_points)
  do i=1,N_points
    integration_points_new(2,i)=integration_points(2,N_points-i+1)
  enddo
  
  compute_non_unif_integral_gaussian = compute_non_unif_integral_gaussian+&
  0.5_f64*compute_non_unif_integral_gaussian_sym(integration_points_new,N_points)
  
  SLL_DEALLOCATE_ARRAY(integration_points_new,ierr)
  
end  function compute_non_unif_integral_gaussian


function compute_non_unif_integral_gaussian_sym(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_gaussian_sym
  sll_real64,dimension(:,:),pointer :: integration_points
  !sll_real64,dimension(:,:),pointer :: integration_points_new
  sll_real64,dimension(:,:),allocatable :: integration_points_new
  sll_int,intent(in) :: N_points
  sll_int :: i,ierr,j,is_center_point,N_points_new
  sll_real64 :: tmp,x1,x2,x3,fval1,fval2,fval3,x4,fval4,dx_int
  sll_int :: N_int,d_gauss,j_gauss
  !sll_real64,dimension(:,:),pointer :: gauss_points
  sll_real64,dimension(:,:),allocatable :: gauss_points
  compute_non_unif_integral_gaussian_sym = 0._f64
  N_int = 2
  d_gauss = 10
  
  SLL_ALLOCATE(gauss_points(2,d_gauss+1),ierr)
  
  if(d_gauss==0)then
    gauss_points(1,1) = 0.5_f64
    gauss_points(2,1) = 1._f64
  endif

  if(d_gauss==2)then
    gauss_points(1,1) = 0.11270166537925831148_f64 
    gauss_points(1,2) = 0.50000000000000000000_f64
    gauss_points(1,3) = 0.88729833462074168852_f64
    gauss_points(2,1) = 0.27777777777777777775_f64
    gauss_points(2,2) = 0.44444444444444444445_f64
    gauss_points(2,3) = 0.27777777777777777778_f64
  endif

  if(d_gauss==4)then
    gauss_points(1,1) = 0.046910077030668003601_f64 
    gauss_points(1,2) = 0.23076534494715845448_f64
    gauss_points(1,3) = 0.50000000000000000000_f64
    gauss_points(1,4) = 0.76923465505284154552_f64
    gauss_points(1,5) = 0.95308992296933199640_f64
    gauss_points(2,1) = 0.11846344252809454382_f64
    gauss_points(2,2) = 0.23931433524968323302_f64
    gauss_points(2,3) = 0.28444444444444444382_f64
    gauss_points(2,4) = 0.23931433524968323408_f64
    gauss_points(2,5) = 0.11846344252809454332_f64
  endif

  if(d_gauss==6)then
    gauss_points(1,1) = 0.025446043828620737737_f64 
    gauss_points(1,2) = 0.12923440720030278007_f64
    gauss_points(1,3) = 0.29707742431130141655_f64
    gauss_points(1,4) = 0.50000000000000000000_f64
    gauss_points(1,5) = 0.70292257568869858345_f64
    gauss_points(1,6) = 0.87076559279969721993_f64
    gauss_points(1,7) = 0.97455395617137926226_f64
    gauss_points(2,1) = 0.064742483084434846538_f64
    gauss_points(2,2) = 0.13985269574463833578_f64
    gauss_points(2,3) = 0.19091502525255949935_f64
    gauss_points(2,4) = 0.20897959183673468797_f64
    gauss_points(2,5) = 0.19091502525255948899_f64
    gauss_points(2,6) = 0.13985269574463833022_f64
    gauss_points(2,7) = 0.064742483084434844832_f64
  endif

  if(d_gauss==8)then
    gauss_points(1,1) = 0.015919880246186955082_f64 
    gauss_points(1,2) = 0.081984446336682102850_f64
    gauss_points(1,3) = 0.19331428364970480135_f64
    gauss_points(1,4) = 0.33787328829809553548_f64
    gauss_points(1,5) = 0.50000000000000000000_f64
    gauss_points(1,6) = 0.66212671170190446452_f64
    gauss_points(1,7) = 0.80668571635029519865_f64
    gauss_points(1,8) = 0.91801555366331789715_f64
    gauss_points(1,9) = 0.98408011975381304492_f64
    gauss_points(2,1) = 0.040637194180787175822_f64
    gauss_points(2,2) = 0.090324080347428468792_f64
    gauss_points(2,3) = 0.13030534820146771147_f64
    gauss_points(2,4) = 0.15617353852000146170_f64
    gauss_points(2,5) = 0.16511967750063012022_f64
    gauss_points(2,6) = 0.15617353852000133681_f64
    gauss_points(2,7) = 0.13030534820146773892_f64
    gauss_points(2,8) = 0.090324080347428458790_f64
    gauss_points(2,9) = 0.040637194180787192773_f64
  endif

  if(d_gauss==10)then
    gauss_points(1,1) = 0.010885670926971503598_f64 
    gauss_points(1,2) = 0.056468700115952350462_f64
    gauss_points(1,3) = 0.13492399721297533795_f64
    gauss_points(1,4) = 0.24045193539659409204_f64
    gauss_points(1,5) = 0.36522842202382751383_f64
    gauss_points(1,6) = 0.50000000000000000000_f64
    gauss_points(1,7) = 0.63477157797617248617_f64
    gauss_points(1,8) = 0.75954806460340590796_f64
    gauss_points(1,9) = 0.86507600278702466205_f64
    gauss_points(1,10) = 0.94353129988404764954_f64
    gauss_points(1,11) = 0.98911432907302849640_f64
    gauss_points(2,1) = 0.027834283558086630522_f64
    gauss_points(2,2) = 0.062790184732454913642_f64
    gauss_points(2,3) = 0.093145105463867582027_f64
    gauss_points(2,4) = 0.11659688229597882542_f64
    gauss_points(2,5) = 0.13140227225513898990_f64
    gauss_points(2,6) = 0.13646254338895855961_f64
    gauss_points(2,7) = 0.13140227225511460324_f64
    gauss_points(2,8) = 0.11659688229599488815_f64
    gauss_points(2,9) = 0.093145105463868674217_f64
    gauss_points(2,10) = 0.062790184732451937905_f64
    gauss_points(2,11) = 0.027834283558088336992_f64
  endif


  
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  i=N_points/2
  is_center_point=0
  if(2*i/=N_points)then
    is_center_point=1
  endif
  !check that the integration points are increasing
  do i=1,N_points-1
    if(integration_points(1,i+1)<=integration_points(1,i))then
      print *,'order problem for integration_points',i,integration_points(1,i),integration_points(1,i+1)
      do j=1,N_points
        print *,j,integration_points(1,j)
      enddo
      stop
    endif
  enddo
  !check for the symmetry  
  tmp=0._f64
  if(is_center_point==0)then
    do i=1,N_points/2
      if(abs(integration_points(1,N_points/2+i)+integration_points(1,N_points/2-i+1))>tmp)then
        tmp=abs(integration_points(1,N_points/2+i)+integration_points(1,N_points/2-i+1))
      endif
    enddo
    if(tmp>1.e-13)then
      print *,'integration_points are not symmetric',tmp
      do j=1,N_points
        print *,j,integration_points(1,j)
      enddo
      stop
    endif
  endif
  
  if(is_center_point==1)then
    tmp=0._f64
    tmp=abs(integration_points(1,(N_points+1)/2))
    do i=1,(N_points-1)/2
      if(abs(integration_points(1,(N_points+1)/2+i)+integration_points(1,(N_points+1)/2-i))>tmp)then
        tmp=abs(integration_points(1,(N_points+1)/2+i)+integration_points(1,(N_points+1)/2-i))
      endif
    enddo
    if(tmp>1.e-14)then
      print *,'integration_points are not symmetric',tmp
      do j=1,N_points
        print *,j,integration_points(1,j)
      enddo
      stop
    endif
  endif
  
  !we will store in a new tab so that we are in the even case
  N_points_new = N_points
  if(is_center_point==0)then
    N_points_new= N_points+1
  endif
  SLL_ALLOCATE(integration_points_new(2,N_points_new),ierr)
  if(is_center_point==1)then
    integration_points_new = integration_points
  endif
  if(is_center_point==0)then
    if(N_points/2+3>N_points)then
      print *,'N_points is too small',N_points
    endif
    integration_points_new(1:2,1:N_points/2) = integration_points(1:2,1:N_points/2)
    integration_points_new(1:2,N_points/2+2:N_points+1) = integration_points(1:2,N_points/2+1:N_points)
    !we have to predict the value of the center
    integration_points_new(1,N_points/2+1)=0._f64
    x1=integration_points(1,N_points/2+1)
    x2=integration_points(1,N_points/2+2)
    x3=integration_points(1,N_points/2+3)
    fval1=integration_points(2,N_points/2+1)*exp(0.5_f64*x1*x1)
    fval2=integration_points(2,N_points/2+2)*exp(0.5_f64*x2*x2)
    fval3=integration_points(2,N_points/2+3)*exp(0.5_f64*x3*x3)
    !print *,x1,x2,x3,fval1,fval2,fval3
    x4=0._f64
    fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))+fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
    fval4=fval4*exp(-0.5_f64*x4*x4)
    integration_points_new(2,N_points/2+1)=fval4
  endif
  
  do i=1,(N_points_new-1)/4
    x1 = integration_points_new(1,(N_points_new-1)/2+2*i-1)
    x2 = integration_points_new(1,(N_points_new-1)/2+2*i)
    x3 = integration_points_new(1,(N_points_new-1)/2+2*i+1)
    fval1 = integration_points_new(2,(N_points_new-1)/2+2*i-1)*exp(0.5_f64*x1*x1)
    fval2 = integration_points_new(2,(N_points_new-1)/2+2*i)*exp(0.5_f64*x2*x2)
    fval3 = integration_points_new(2,(N_points_new-1)/2+2*i+1)*exp(0.5_f64*x3*x3)
    tmp=0._f64
    dx_int=(x3-x1)/real(N_int,f64)
    do j=1,N_int
      do j_gauss=1,d_gauss+1
        !x4 =x1+(real(j,f64)-0.5_f64)*dx_int
        x4 =x1+(real(j-1,f64)+gauss_points(1,j_gauss))*dx_int
        fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))+fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
        fval4=fval4*exp(-0.5_f64*x4*x4)
        tmp=tmp+fval4*gauss_points(2,j_gauss)
      enddo
    enddo  
    tmp=tmp*dx_int
    !print *,i,x1,x2,x3
    compute_non_unif_integral_gaussian_sym =compute_non_unif_integral_gaussian_sym+tmp
    !integration_points_middle(1,i)=0.5_f64*(x1+x2)
  enddo
  j=(N_points_new-1)/4
  if(2*j/=(N_points_new-1)/2)then
    !print *,2*j,(N_points_new-1)/2
    if(2*j+1/=(N_points_new-1)/2)then
      print *,'Problem concerning N_points',2*j+1,(N_points_new-1)/2
    endif
    x1 = integration_points_new(1,N_points_new-2)
    x2 = integration_points_new(1,N_points_new-1)
    x3 = integration_points_new(1,N_points_new)
    fval1 = integration_points_new(2,N_points_new-2)*exp(0.5_f64*x1*x1)
    fval2 = integration_points_new(2,N_points_new-1)*exp(0.5_f64*x2*x2)
    fval3 = integration_points_new(2,N_points_new)*exp(0.5_f64*x3*x3)
    x4 = 2._f64*x3-x2
    fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))+fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
    fval4=fval4*exp(-0.5_f64*x4*x4)


    
    x1 = integration_points_new(1,N_points_new-1)
    x2 = integration_points_new(1,N_points_new)
    x3 = 2._f64*x2-x1
    fval1 = integration_points_new(2,N_points_new-1)*exp(0.5_f64*x1*x1)
    fval2 = integration_points_new(2,N_points_new)*exp(0.5_f64*x2*x2)
    fval3 = fval4*exp(0.5_f64*x3*x3)
    
    !print *,x1,x2,x3
    !print *,fval1,fval2,fval3
    !print *,x1,x2,x3
    !print *,fval1,fval2,fval3
    
    !stop
    
    tmp=0._f64
    N_int =N_int*2
    dx_int=(x3-x1)/real(N_int,f64)
    do j=1,N_int
      do j_gauss=1,d_gauss+1
        !x4 =x1+(real(j,f64)-0.5_f64)*dx_int
        x4 =x1+(real(j-1,f64)+gauss_points(1,j_gauss))*dx_int
        fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))&
        +fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
        fval4=fval4*exp(-0.5_f64*x4*x4)
        tmp=tmp+fval4*gauss_points(2,j_gauss)
      enddo
    enddo
    tmp=tmp*dx_int
    !print *,i,x1,x2,x3
    compute_non_unif_integral_gaussian_sym =compute_non_unif_integral_gaussian_sym+tmp
  endif
  
  
  
  
  compute_non_unif_integral_gaussian_sym = 2._f64*compute_non_unif_integral_gaussian_sym
  
  SLL_DEALLOCATE_ARRAY(integration_points_new,ierr)
  SLL_DEALLOCATE_ARRAY(gauss_points,ierr)

  !SLL_DEALLOCATE(integration_points_new,ierr)
  !SLL_DEALLOCATE(gauss_points,ierr)
  
end  function compute_non_unif_integral_gaussian_sym



function compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib,index_contrib,value_contrib,Ntot,N1,N2)
  sll_int :: compute_contrib_rho
  sll_int,intent(in)         ::  Ntot,N1,N2,N_int1,N_int2
  sll_real64,dimension(2,N_int2+1,N_int1+1)  ::  integration_points
  sll_int,dimension(-1:N2+1,N_int1+1)    ::  size_contrib
  sll_int,dimension(Ntot)         ::  index_contrib
  sll_real64,dimension(Ntot)         :: value_contrib
  sll_int                    ::  i,j,ii,jj,i_loc,j_loc
  sll_real64                  ::x1,x2,w(1:2,-1:2),xx
  sll_int,dimension(:,:,:),pointer :: point_flag
  sll_int,dimension(:,:),pointer :: line_flag
  sll_real64,dimension(:,:),pointer ::  point_val
  sll_int :: N_size,s,err
  SLL_ALLOCATE(point_flag(1:2,-1:N1+1,-1:N2+1),err)
  SLL_ALLOCATE(line_flag(1:2,-1:N2+1),err)
  SLL_ALLOCATE(point_val(-1:N1+1,-1:N2+1),err)
  N_size = 0
  line_flag=0
  point_flag=0
  s=0
  size_contrib=0
  do i=1,N_int1+1
    !buf=0
    do j=1,N_int2+1
      x1 = integration_points(1,j,i)
      x2 = integration_points(2,j,i)
      !treatment of boundary conditions
      ii=floor(x1*real(N1,f64))
      if(x1>=1._f64)then
        x1=1._f64
        ii=N1-1
      endif
      if(x1<=0._f64)then
        x1=0._f64
        ii=0
      endif
      jj=floor(x2*real(N2,f64))
      if(x2>=1._f64)then
        x2=1._f64
        jj=N2-1
      endif
      if(x2<=0._f64)then
        x2=0._f64
        jj=0
      endif
      
      xx=(x1*real(N1,f64))-real(ii,f64)
      
      if((xx>1_f64).or.(xx<0._f64))then
        print *,'x1=',xx
        stop
      endif
      
      w(1,-1)=(1.0_f64-xx)*(1.0_f64-xx)*(1.0_f64-xx)/6._f64    
      w(1,0)=2._f64/3._f64+xx*xx*xx/2._f64-xx*xx
      w(1,1)=xx*xx/2._f64+xx/2._f64+1._f64/6._f64-xx*xx*xx/2._f64
      w(1,2)=xx*xx*xx/6._f64
      xx=(x2*real(N2,f64))-real(jj,f64)      

      if((xx>1_f64).or.(xx<0._f64))then
        print *,'x2=',xx
        stop
      endif


      w(2,-1)=(1.0_f64-xx)*(1.0_f64-xx)*(1.0_f64-xx)/6._f64    
      w(2,0)=2._f64/3._f64+xx*xx*xx/2._f64-xx*xx
      w(2,1)=xx*xx/2._f64+xx/2._f64+1._f64/6._f64-xx*xx*xx/2._f64
      w(2,2)=xx*xx*xx/6._f64
      
      do j_loc=-1,2
        do i_loc=-1,2
          if(line_flag(1,jj+j_loc)/=i)then
            line_flag(1,jj+j_loc)=i
            line_flag(2,jj+j_loc)=-1
         endif
          if(point_flag(1,ii+i_loc,jj+j_loc)/=i)then
            size_contrib(jj+j_loc,i)=size_contrib(jj+j_loc,i)+1
            N_size  = N_size +1
            point_flag(1,ii+i_loc,jj+j_loc)=i
            point_val(ii+i_loc,jj+j_loc) = 0._f64
            point_flag(2,ii+i_loc,jj+j_loc) = line_flag(2,jj+j_loc)
            line_flag(2,jj+j_loc) = ii+i_loc
          endif          
          point_val(ii+i_loc,jj+j_loc)=point_val(ii+i_loc,jj+j_loc)+w(1,i_loc)*w(2,j_loc)
        enddo
      enddo
      
    enddo
    
    !fill the contrib arrays
    !do j=-1,N2+1
    !  !print *,i,j,size_contrib(j,i),line_flag(2,j)
    !  ii=line_flag(2,j)
    !  do jj=1,size_contrib(j,i)
    !    !print *,"##",ii,point_val(ii,j)
    !    ii=point_flag(2,ii,j) 
    !  enddo      
    !enddo
    !stop
    
    do j = -1,N2+1
      if(size_contrib(j,i)>=1)then
        ii=line_flag(2,j)
      endif
      do jj=1,size_contrib(j,i)
        s=s+1
        if(s<=Ntot)then          
          index_contrib(s) = ii
          value_contrib(s) = point_val(ii,j)
        endif
        ii=point_flag(2,ii,j) 
      enddo
    enddo
    
  enddo
  
  SLL_DEALLOCATE(point_flag,err)
  SLL_DEALLOCATE(line_flag,err)
  SLL_DEALLOCATE(point_val,err)
  !deallocate(point_flag,line_flag,point_val)
  
  if(s/=N_size)then
    print *,'s=',s,'N_size=',N_size
    !stop
  endif
  
  compute_contrib_rho = s
  
end function compute_contrib_rho




end module contrib_rho_module

