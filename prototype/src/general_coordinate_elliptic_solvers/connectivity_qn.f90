module connectivity_module
  implicit none

  integer, parameter :: CONNECT_PERIODIC = 0, CONNECT_DIRICHLET = 1

contains
!----------------------------------------------------------------------------					
  subroutine initconnectivity( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    local_spline_indices, &
    global_spline_indices, &
    local_to_global_spline_indices )
    
    integer :: num_cells1, num_cells2, spline_degree1, spline_degree2
    integer :: bc_left, bc_right, bc_bottom, bc_top
    integer :: nb_spl_x,nb_spl_y
    integer, dimension(:,:) :: local_spline_indices
    integer, dimension(:)   :: global_spline_indices
    integer, dimension(:,:) :: local_to_global_spline_indices

    integer :: BC  ! 0 if periodic-Dirichlet and 1 if Dirichlet-Dirichlet 
     integer :: li_e,li_i_, li_j_,t1,t2
    integer :: li_i, li_j, li_iloc, li_jloc, li_Bloc, li_B,maille    
    
    !print*, 'local_spline_indices=', local_spline_indices(1,1)

    ! loop over elements
    do li_j_ = 1, num_cells2    
       do li_i_=1, num_cells1  
          maille = num_cells1*(li_j_-1) + li_i_
          do li_jloc = 0 , spline_degree2
             do li_iloc = 0, spline_degree1
                li_B = (li_j_ -1)*(num_cells1+spline_degree1) + li_i_ + &
                          li_jloc*(num_cells1+spline_degree1) + li_iloc 
                li_Bloc = li_jloc * (spline_degree1 + 1) + li_iloc + 1
                local_spline_indices(li_Bloc, maille) = li_B
             end do
          end do
       end do
    end do

    ! INITIALIIZING THE ID ARRAY
    !select case ( ai_TYPE_PBC)

    !case ( BC_XI_DIR_ETA_PER )
    nb_spl_x= num_cells1 + spline_degree1
    nb_spl_y= num_cells2 + spline_degree2 
    
   if( (bc_left == CONNECT_PERIODIC) .and. (bc_right == CONNECT_PERIODIC).and.&
       (bc_bottom == CONNECT_PERIODIC) .and. (bc_top == CONNECT_PERIODIC) ) then
      call xi_PER_eta_PER_init( &
           nb_spl_x,&
           nb_spl_y,&
           spline_degree1,&
           spline_degree2,&
           global_spline_indices)
        
    else if ((bc_left==CONNECT_PERIODIC).and.(bc_right==CONNECT_PERIODIC).and. &
             (bc_bottom==CONNECT_DIRICHLET).and.(bc_top==CONNECT_DIRICHLET))then
       
       call xi_PER_eta_DIR_init(&
            nb_spl_x,&
            nb_spl_y,&
            spline_degree1,&
            global_spline_indices)

    else if((bc_left == CONNECT_DIRICHLET).and.(bc_right==CONNECT_DIRICHLET).and. &
             (bc_bottom==CONNECT_PERIODIC).and.(bc_top==CONNECT_PERIODIC) )then

       call xi_DIR_eta_PER_init(&
            nb_spl_x,&
            nb_spl_y,&
            spline_degree2,&
            global_spline_indices)

    else if((bc_left == CONNECT_DIRICHLET).and.(bc_right==CONNECT_DIRICHLET).and. &
            (bc_bottom==CONNECT_DIRICHLET).and.(bc_top==CONNECT_DIRICHLET))then

       call xi_DIR_eta_DIR_init(nb_spl_x,nb_spl_y,global_spline_indices)
    end if
    
    ! INITIALIIZING THE LOCATION MATRIX LM ARRAY
    call initLM(&
         num_cells1,&
         num_cells2,&
         spline_degree1,&
         spline_degree2,&
         local_spline_indices,&
         global_spline_indices,&
         local_to_global_spline_indices)

  end subroutine initconnectivity

  !--------------------------------------------------------------------------
  subroutine xi_DIR_eta_PER_init( &
       num_cells1,&
       num_cells2,&
       spline_degree2,&
       global_spline_indices)
    
    integer :: li_d,num_cells1,num_cells2,spline_degree2
    integer :: li_i, li_j, li_A
    integer :: li_dof,ai_sizePB
    integer :: li_L
    integer, dimension(:) :: global_spline_indices
    
    li_d = 0
    
    do li_j = 1, num_cells2
       do li_i = 1, num_cells1
          li_A = li_i + num_cells1*(li_j-1)
          
          if ((li_i == 1) .OR. (li_i == num_cells1)) then
             global_spline_indices(li_A) = 0
             
          else
             
             if (li_j /= num_cells2) then
                
                if (global_spline_indices(li_A) == 0) then
                   
                   li_d = li_d + 1
                   
                   global_spline_indices(li_A) = li_d
                   
                end if
                
                if ( (1 <= li_j ) .AND. ( li_j <= spline_degree2 ) ) then
                   
                   li_L = (num_cells2 - spline_degree2) * num_cells1
                   global_spline_indices(li_A + li_L) =&
                        global_spline_indices(li_A)

                end if
             end if
          end if
       end do
    end do

  end subroutine xi_DIR_eta_PER_init

  !-------------------------------------------------------------------------
  subroutine xi_PER_eta_DIR_init( &
       num_cells1,&
       num_cells2,&
       spline_degree1,&
       global_spline_indices)
    
    integer :: li_d,num_cells1,num_cells2,spline_degree1
    integer :: li_i, li_j, li_A
    integer :: li_dof,ai_sizePB
    integer, dimension(:) :: global_spline_indices

    li_d = 0
    do li_j = 1, num_cells2
       do li_i = 1, num_cells1
          li_A = li_i + num_cells1*(li_j-1)

          if ((li_j == 1) .OR. (li_j == num_cells2)) then

             global_spline_indices(li_A) = 0
             
          else

             if (li_i /= num_cells1) then

                if (global_spline_indices(li_A) == 0) then

                   li_d = li_d + 1
                   
                   global_spline_indices(li_A) = li_d

                end if

                if ( (1 <= li_i ) .AND. ( li_i <= spline_degree1 ) ) then

                   global_spline_indices(li_A + num_cells1 - spline_degree1) = &
                        global_spline_indices(li_A)

                end if
             end if
          end if
       end do
    end do
  end subroutine xi_PER_eta_DIR_init

!-----------------------------------------------------------------------

  subroutine xi_DIR_eta_DIR_init(&
       num_cells1,&
       num_cells2,&
       global_spline_indices)

    integer :: li_d,num_cells1,num_cells2
    integer :: li_i, li_j, li_A
    integer :: li_dof,ai_sizePB
    integer, dimension(:) :: global_spline_indices
    
    li_d = 0

    do li_j = 2, num_cells2-1
        do li_i = 2, num_cells1-1
            li_A = li_i + num_cells1*(li_j-1)
            li_d = li_d + 1
            global_spline_indices(li_A) = li_d
        end do
    end do
  end subroutine xi_DIR_eta_DIR_init


  subroutine xi_PER_eta_PER_init(&
       num_cells1,&
       num_cells2,&
       spline_degree1,&
       spline_degree2,&
       global_spline_indices)

    integer :: li_d,num_cells1,num_cells2,spline_degree1,spline_degree2
    integer :: li_i, li_j, li_A, li_L
    integer :: li_dof,ai_sizePB
    integer, dimension(:) :: global_spline_indices
    
    li_d = 0

    do li_j = 1, num_cells2
       
       do li_i = 1, num_cells1
          
          li_A = li_i + num_cells1*(li_j-1)
          
          !print*, li_A
          
          
          if (li_i /= num_cells1 .and. li_j/= num_cells2 ) then

             if (global_spline_indices(li_A) == 0) then
                
                li_d = li_d + 1
                
                global_spline_indices(li_A) = li_d
                
             end if

             if ( (1 <= li_i ) .AND. ( li_i <= spline_degree1 ) ) then

                global_spline_indices(li_A + num_cells1 - spline_degree1) = &
                     global_spline_indices(li_A)
             end if

             if ( (1 <= li_j ) .AND. ( li_j <= spline_degree2 ) ) then

                li_L = (num_cells2 - spline_degree2) * num_cells1
                global_spline_indices(li_A  + li_L) = &
                     global_spline_indices(li_A)
                !print*, 'hey',li_A + li_L, li_L
             end if
          end if
       end do
    end do

    global_spline_indices(num_cells2*num_cells1) = &
         global_spline_indices(num_cells1*(num_cells2-1) + spline_degree2)

  end subroutine xi_PER_eta_PER_init


!------------------------------------------------------------------------
  subroutine initLM( &
       num_cells1,&
       num_cells2,&
       spline_degree1,&
       spline_degree2,&
       local_spline_indices,&
       global_spline_indices,&
       local_to_global_spline_indices)

    integer :: li_e, li_b,num_cells1,num_cells2,spline_degree1,spline_degree2
    integer, dimension(:) :: global_spline_indices
    integer, dimension(:,:) :: local_spline_indices
    integer,dimension(:,:) :: local_to_global_spline_indices
    
    do li_e = 1, num_cells1*num_cells2
       
       do li_b = 1, (spline_degree1+1)*(spline_degree2+1)
          
          ! print*, global_spline_indices(local_spline_indices(li_b, li_e))
          local_to_global_spline_indices(li_b, li_e) = &
               global_spline_indices(local_spline_indices(li_b, li_e))
       end do
    end do
  end subroutine initLM
      
end module connectivity_module
