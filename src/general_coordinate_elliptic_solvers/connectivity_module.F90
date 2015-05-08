module connectivity_module
#include "sll_working_precision.h"

implicit none

sll_int32, parameter :: CONNECT_PERIODIC  = 0
sll_int32, parameter :: CONNECT_DIRICHLET = 1
sll_int32, parameter :: CONNECT_NEUMANN   = 2

contains

subroutine initconnectivity( num_cells1,                    &
                             num_cells2,                    &
                             spline_degree1,                &
                             spline_degree2,                &
                             bc1_min,                       &
                             bc1_max,                       &
                             bc2_min,                       &
                             bc2_max,                       &
                             local_spline_indices,          &
                             global_spline_indices,         &
                             local_to_global_spline_indices )
  
sll_int32, intent(in) :: num_cells1
sll_int32, intent(in) :: num_cells2
sll_int32, intent(in) :: spline_degree1
sll_int32, intent(in) :: spline_degree2
sll_int32, intent(in) :: bc1_min
sll_int32, intent(in) :: bc1_max
sll_int32, intent(in) :: bc2_min
sll_int32, intent(in) :: bc2_max

sll_int32, dimension(:,:), intent(out) :: local_spline_indices
sll_int32, dimension(:)  , intent(out) :: global_spline_indices
sll_int32, dimension(:,:), intent(out) :: local_to_global_spline_indices

sll_int32 :: nb_spl_x
sll_int32 :: nb_spl_y
sll_int32 :: i, j
sll_int32 :: ii, jj
sll_int32 :: Bloc, B, maille    
  
do j = 1, num_cells2    
  do i=1, num_cells1  
    maille = num_cells1*(j-1) + i
    do jj = 0 , spline_degree2
      do ii = 0, spline_degree1
        B = (j-1)*(num_cells1+spline_degree1) + i + &
                  jj*(num_cells1+spline_degree1) + ii 
        Bloc = jj * (spline_degree1 + 1) + ii + 1
        local_spline_indices(Bloc, maille) = B
      end do
    end do
  end do
end do

nb_spl_x= num_cells1 + spline_degree1
nb_spl_y= num_cells2 + spline_degree2 
  
if( bc1_min == CONNECT_PERIODIC .and. bc1_max == CONNECT_PERIODIC .and.&
    bc2_min == CONNECT_PERIODIC .and. bc2_max == CONNECT_PERIODIC ) then

  call xi_PER_eta_PER_init( nb_spl_x,            &
                            nb_spl_y,            &
                            spline_degree1,      &
                            spline_degree2,      &
                            global_spline_indices)
      
else if (bc1_min == CONNECT_PERIODIC  .and. bc1_max == CONNECT_PERIODIC .and. &
         bc2_min == CONNECT_DIRICHLET .and. bc2_max == CONNECT_DIRICHLET) then
     
  call xi_PER_eta_DIR_init( nb_spl_x,            &
                            nb_spl_y,            &
                            spline_degree1,      &
                            global_spline_indices)

else if(bc1_min == CONNECT_DIRICHLET .and. bc1_max == CONNECT_DIRICHLET .and. &
        bc2_min == CONNECT_PERIODIC  .and. bc2_max == CONNECT_PERIODIC) then

  call xi_DIR_eta_PER_init( nb_spl_x,            &
                            nb_spl_y,            &
                            spline_degree2,      &
                            global_spline_indices)

elseif(bc1_min == CONNECT_NEUMANN  .and. bc1_max == CONNECT_DIRICHLET .and. &
       bc2_min == CONNECT_PERIODIC .and. bc2_max == CONNECT_PERIODIC ) then

  call xi_NEU_DIR_eta_PER_init( nb_spl_x,            &
                                nb_spl_y,            &
                                spline_degree2,      &
                                global_spline_indices)       


else if(bc1_min == CONNECT_DIRICHLET .and. bc1_max==CONNECT_DIRICHLET .and. &
        bc2_min == CONNECT_DIRICHLET .and. bc2_max==CONNECT_DIRICHLET) then

  call xi_DIR_eta_DIR_init(nb_spl_x,nb_spl_y,global_spline_indices)

end if
  
call initLM( num_cells1,                   &
             num_cells2,                   &
             spline_degree1,               &
             spline_degree2,               &
             local_spline_indices,         &
             global_spline_indices,        &
             local_to_global_spline_indices)

end subroutine initconnectivity


subroutine xi_DIR_eta_PER_init( num_cells1,          &
                                num_cells2,          &
                                spline_degree2,      &
                                global_spline_indices)
  
sll_int32,               intent(in)  :: num_cells1
sll_int32,               intent(in)  :: num_cells2
sll_int32,               intent(in)  :: spline_degree2
sll_int32, dimension(:), intent(out) :: global_spline_indices

sll_int32 :: d
sll_int32 :: i
sll_int32 :: j
sll_int32 :: A
sll_int32 :: L

d = 0

do j = 1, num_cells2
do i = 1, num_cells1

  A = i + num_cells1*(j-1)
  if ((i == 1) .OR. (i == num_cells1)) then
    global_spline_indices(A) = 0
  else
    if (j /= num_cells2) then
      if (global_spline_indices(A) == 0) then
        d = d + 1
        global_spline_indices(A) = d
      end if
      if ( 1<=j .AND. j<=spline_degree2 ) then
        L = (num_cells2 - spline_degree2) * num_cells1
        global_spline_indices(A+L) = global_spline_indices(A)
      end if
    end if
  end if

end do
end do

end subroutine xi_DIR_eta_PER_init

subroutine xi_NEU_DIR_eta_PER_init( num_cells1,          &
                                    num_cells2,          &
                                    spline_degree2,      &
                                    global_spline_indices)
  
sll_int32              , intent(in)  :: num_cells1
sll_int32              , intent(in)  :: num_cells2
sll_int32              , intent(in)  :: spline_degree2
sll_int32, dimension(:), intent(out) :: global_spline_indices

sll_int32 :: d
sll_int32 :: i
sll_int32 :: j
sll_int32 :: A
sll_int32 :: L

d = 0
do j = 1, num_cells2
do i = 1, num_cells1
  A = i + num_cells1*(j-1)
  if (i == num_cells1) then
    global_spline_indices(A) = 0
  else
    if (j /= num_cells2) then
      if (global_spline_indices(A) == 0) then
        d = d + 1
        global_spline_indices(A) = d
      end if
      if ( (1 <= j ) .AND. ( j <= spline_degree2 ) ) then
        L = (num_cells2 - spline_degree2) * num_cells1
        global_spline_indices(A+L) = global_spline_indices(A)
      end if
    end if
  end if
end do
end do

end subroutine xi_NEU_DIR_eta_PER_init

subroutine xi_PER_eta_DIR_init( num_cells1,          &
                                num_cells2,          &
                                spline_degree1,      &
                                global_spline_indices)
  
sll_int32,               intent(in)  :: num_cells1
sll_int32,               intent(in)  :: num_cells2
sll_int32,               intent(in)  :: spline_degree1
sll_int32, dimension(:), intent(out) :: global_spline_indices

sll_int32 :: d
sll_int32 :: i
sll_int32 :: j
sll_int32 :: A

d = 0
do j = 1, num_cells2
do i = 1, num_cells1
  A = i + num_cells1*(j-1)
  if ((j == 1) .OR. (j == num_cells2)) then
    global_spline_indices(A) = 0
  else
    if (i /= num_cells1) then
      if (global_spline_indices(A) == 0) then
        d = d + 1
        global_spline_indices(A) = d
      end if
      if ( 1 <= i .AND. i <= spline_degree1 ) then
        global_spline_indices(A+num_cells1-spline_degree1) = global_spline_indices(A)
      end if
    end if
  end if
end do
end do

end subroutine xi_PER_eta_DIR_init

subroutine xi_DIR_eta_DIR_init( num_cells1,          &
                                num_cells2,          &
                                global_spline_indices)

sll_int32,               intent(in)  :: num_cells1
sll_int32,               intent(in)  :: num_cells2
sll_int32, dimension(:), intent(out) :: global_spline_indices

sll_int32 :: d
sll_int32 :: i
sll_int32 :: j
sll_int32 :: A

d = 0
do j = 2, num_cells2-1
do i = 2, num_cells1-1
  A = i + num_cells1*(j-1)
  d = d + 1
  global_spline_indices(A) = d
end do
end do

end subroutine xi_DIR_eta_DIR_init

subroutine xi_eta_init( num_cells1,          &
                        num_cells2,          &
                        global_spline_indices)

sll_int32,               intent(in)  :: num_cells1
sll_int32,               intent(in)  :: num_cells2
sll_int32, dimension(:), intent(out) :: global_spline_indices

sll_int32 :: d
sll_int32 :: i
sll_int32 :: j
sll_int32 :: A


d = 0
do j = 1, num_cells2
do i = 1, num_cells1
  A = i + num_cells1*(j-1)
  d = d + 1
  global_spline_indices(A) = d
end do
end do

end subroutine xi_eta_init

subroutine xi_PER_eta_PER_init( num_cells1,          &
                                num_cells2,          &
                                spline_degree1,      &
                                spline_degree2,      &
                                global_spline_indices)

sll_int32              , intent(in)  :: num_cells1
sll_int32              , intent(in)  :: num_cells2
sll_int32              , intent(in)  :: spline_degree1
sll_int32              , intent(in)  :: spline_degree2
sll_int32, dimension(:), intent(out) :: global_spline_indices

sll_int32               :: d
sll_int32               :: i
sll_int32               :: j
sll_int32               :: A
sll_int32               :: L

d = 0
do j = 1, num_cells2
do i = 1, num_cells1

  A = i + num_cells1*(j-1)
  if (i /= num_cells1 .and. j/= num_cells2 ) then
    if (global_spline_indices(A) == 0) then
      d = d+1
      global_spline_indices(A) = d
    end if
    if ( (1 <= i ) .AND. ( i <= spline_degree1 ) ) then
      global_spline_indices(A+num_cells1-spline_degree1) = global_spline_indices(A)
    end if
    if ( (1 <= j ) .AND. ( j <= spline_degree2 ) ) then
      L = (num_cells2 - spline_degree2) * num_cells1
      global_spline_indices(A+L) = global_spline_indices(A)
    end if
  end if

end do
end do

global_spline_indices(num_cells2*num_cells1) = &
     global_spline_indices(num_cells1*(num_cells2-1) + spline_degree2)

end subroutine xi_PER_eta_PER_init

subroutine initLM( num_cells1,                   &
                   num_cells2,                   &
                   spline_degree1,               &
                   spline_degree2,               &
                   local_spline_indices,         &
                   global_spline_indices,        &
                   local_to_global_spline_indices)

sll_int32                , intent(in)  :: num_cells1
sll_int32                , intent(in)  :: num_cells2
sll_int32                , intent(in)  :: spline_degree1
sll_int32                , intent(in)  :: spline_degree2
sll_int32, dimension(:)  , intent(out) :: global_spline_indices
sll_int32, dimension(:,:), intent(out) :: local_spline_indices
sll_int32, dimension(:,:), intent(out) :: local_to_global_spline_indices

sll_int32                 :: e
sll_int32                 :: b

do e = 1, num_cells1*num_cells2
  do b = 1, (spline_degree1+1)*(spline_degree2+1)
    local_to_global_spline_indices(b,e) = global_spline_indices(local_spline_indices(b,e))
  end do
end do

end subroutine initLM

end module connectivity_module

!subroutine initconnectivity_bis( num_cells1, &
!                                 num_cells2, &
!                                 spline_degree1, &
!                                 spline_degree2, &
!                                 local_spline_indices, &
!                                 global_spline_indices, &
!                                 local_to_global_spline_indices )
!  
!  sll_int32 :: num_cells1, num_cells2, spline_degree1, spline_degree2
!  sll_int32 :: nb_spl_x,nb_spl_y
!  sll_int32, dimension(:,:) :: local_spline_indices
!  sll_int32, dimension(:)   :: global_spline_indices
!  sll_int32, dimension(:,:) :: local_to_global_spline_indices
!  
! ! sll_int32 :: BC  ! 0 if periodic-Dirichlet and 1 if Dirichlet-Dirichlet 
!  sll_int32 :: i, j!,t1,t2
!  sll_int32 :: ii, jj, Bloc, B,maille    
!  
!  !print*, 'local_spline_indices=', local_spline_indices(1,1)
!
!  ! loop over elements
! do j = 1, num_cells2    
!     do i=1, num_cells1  
!        maille = num_cells1*(j-1) + i
!        do jj = 0 , spline_degree2
!           do ii = 0, spline_degree1
!              B = (j -1)*(num_cells1+spline_degree1) + i + &
!                        jj*(num_cells1+spline_degree1) + ii 
!              Bloc = jj * (spline_degree1 + 1) + ii + 1
!              local_spline_indices(Bloc, maille) = B
!           end do
!        end do
!     end do
!  end do
!
! ! print*, local_spline_indices
!  ! INITIALIIZING THE ID ARRAY
!  !select case ( ai_TYPE_PBC)
!
!  !case ( BC_XI_DIR_ETA_PER )
!  nb_spl_x= num_cells1 + spline_degree1
!  nb_spl_y= num_cells2 + spline_degree2
!  
!  call xi_eta_init(&
!       nb_spl_x,&
!       nb_spl_y,&
!       global_spline_indices)
!  
!  ! INITIALIIZING THE LOCATION MATRIX LM ARRAY
!  call initLM(&
!       num_cells1,&
!       num_cells2,&
!       spline_degree1,&
!       spline_degree2,&
!       local_spline_indices,&
!       global_spline_indices,&
!       local_to_global_spline_indices)
!
!end subroutine initconnectivity_bis
