program VP_test

  use sll_nu_cart_mesh
  use gausslobatto_mod

  implicit none

  type(gausslobatto1d) :: gl
  integer :: i
  
  do i=4,4

     print*,i
     call init_gausslobatto_1d(i,gl)
     print*,gl%node
     print*,gl%weight
     call delete_gausslobatto_1d(gl)
     print*,' '

  end do

end program VP_test
