module sll_reduction_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_remapper
  use sll_logical_meshes

contains

  !--------------------------------------------------
  ! Generic function for computing charge density
  ! we should also add a choice for the integration 
  ! also should be generalized for more complicated data
  ! as here array of values f_0,\dots,f_N 
  !---------------------------------------------------  
  
  subroutine compute_reduction_4d_to_3d(&
      m_x4, &
      data_4d, &
      data_3d, &
      layout4d, &
      m4d)
    type(sll_logical_mesh_1d), pointer    :: m_x4
    sll_real64, dimension(:,:,:,:), intent(in)    :: data_4d
    sll_real64, dimension(:,:,:)  , intent(out) :: data_3d
    type(layout_4D), pointer, optional    :: layout4d
    type(sll_logical_mesh_4d), pointer, optional    :: m4d
    sll_int32  :: np_x1_loc
    sll_int32  :: np_x2_loc
    sll_int32  :: np_x3
    sll_int32  :: np_x4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_real64 :: delta_x4, tmp
    sll_int32  :: loc4d_sz_x1
    sll_int32  :: loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3
    sll_int32  :: loc4d_sz_x4

    if(present(layout4d))then
      call compute_local_sizes_4d( layout4d, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )
    else if(present(m4d)) then          
      loc4d_sz_x1 = m4d%num_cells1+1
      loc4d_sz_x2 = m4d%num_cells2+1
      loc4d_sz_x3 = m4d%num_cells3+1
      loc4d_sz_x4 = m4d%num_cells4+1      
    else
      print *,'#Problem in compute_reduction_4d_to_3d'
      print *,'#logical_mesh4d or layout_4D are not provided'
      stop 
    endif
    
    
    if(loc4d_sz_x1>size(data_4d,1))then
      print *,'#Problem for size1 in compute_reduction_4d_to_3d'
      stop
    endif
    if(loc4d_sz_x2>size(data_4d,2))then
      print *,'#Problem for size2 in compute_reduction_4d_to_3d'
      stop
    endif
    if(loc4d_sz_x3>size(data_4d,3))then
      print *,'#Problem for size3 in compute_reduction_4d_to_3d'
      stop
    endif
    if(loc4d_sz_x4>size(data_4d,4))then
      print *,'#Problem for size3 in compute_reduction_4d_to_3d'
      stop
    endif

    if((loc4d_sz_x4).ne.(m_x4%num_cells+1))then
      print *,'#Problem for size in compute_reduction_4d_to_3d'
    endif

    do iloc3 = 1,loc4d_sz_x3
      do iloc2 = 1,loc4d_sz_x2
        do iloc1 = 1,loc4d_sz_x1
          tmp = 0.5_f64*(data_4d(iloc1,iloc2,iloc3,1)&
            +data_4d(iloc1,iloc2,iloc3,loc4d_sz_x4))
          do iloc4 = 2,loc4d_sz_x4-1
            tmp = tmp + data_4d(iloc1,iloc2,iloc3,iloc4)
          end do
          data_3d(iloc1,iloc2,iloc3) = tmp*m_x4%delta_eta
        end do
      end do
    end do
  end subroutine compute_reduction_4d_to_3d

end module sll_reduction_module