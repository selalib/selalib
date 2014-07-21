module connectivity_Module
  implicit none

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
    integer :: BC  ! 0 if periodic-Dirichlet and 1 if Dirichlet-Dirichlet 
    integer, dimension(:,:) :: local_spline_indices
    integer, dimension(:)   :: global_spline_indices
    integer, dimension(:,:) :: local_to_global_spline_indices
    !LOCAL VARIABLES
    integer :: li_e,li_i_, li_j_,t1,t2
    integer :: li_i, li_j, li_iloc, li_jloc, li_Bloc, li_B,maille    
    
    !print*, 'api_IEN=', api_IEN(1,1)

    ! loop over elements
    do li_j_ = 1, Ny    
        do li_i_=1, Nx  
            
            maille = Nx*(li_j_-1) + li_i_
            do li_jloc = 0 , spline_degree

                do li_iloc = 0, spline_degree
                    

                    li_B = (li_j_ -1)*(Nx+spline_degree) + li_i_ + li_jloc*(Nx+spline_degree) + li_iloc 
                    

                    li_Bloc = li_jloc * (spline_degree + 1) + li_iloc + 1
                    
                    api_IEN(li_Bloc, maille) = li_B

                end do

            end do

        end do
    end do           

    ! INITIALIIZING THE ID ARRAY
    !select case ( ai_TYPE_PBC)

    !case ( BC_XI_DIR_ETA_PER )

    nb_spl_x= Nx + spline_degree
    nb_spl_y= Ny + spline_degree 
    
    
    if ( BC == 0 ) then 
        
        call xi_DIR_eta_PER_init(nb_spl_x,nb_spl_y,spline_degree,api_ID)
        
    else if ( BC == 1 ) then 
        
        call xi_DIR_eta_DIR_init(nb_spl_x,nb_spl_y,spline_degree,api_ID)
    
    end if 
    
!    call xi_PER_eta_DIR_init(nb_spl_x,nb_spl_y,spline_degree,api_ID)
    

    ! INITIALIIZING THE LOCATION MATRIX LM ARRAY
    call initLM(Nx,Ny,spline_degree,api_IEN,api_ID,api_LM)

end subroutine initconnectivity

!----------------------------------------------------------------------------------------------
subroutine xi_PER_eta_DIR_init(Nx,Ny,spline_degree,api_ID)
    implicit none
    !LOCAL VARIABLES
    integer :: li_d,Nx,Ny,spline_degree
    integer :: li_i, li_j, li_A
    integer :: li_dof,ai_sizePB
    integer :: li_L
    integer, dimension(:) :: api_ID

    li_d = 0

    do li_j = 1, Ny

        do li_i = 1, Nx
            
            li_A = li_i + Nx*(li_j-1)

            if ((li_i == 1) .OR. (li_i == Nx)) then

                api_ID(li_A) = 0

            else

                if (li_j /= Ny) then

                    if (api_ID(li_A) == 0) then

                        li_d = li_d + 1

                        api_ID(li_A) = li_d

                    end if

                    if ( (1 <= li_j ) .AND. ( li_j <= spline_degree ) ) then

                        li_L = (Ny - 1) * Nx
                        api_ID(li_A + li_L) = api_ID(li_A)

                    end if

                end if

            end if

        end do

    end do


end subroutine xi_PER_eta_DIR_init
!----------------------------------------------------------------------------------------------
subroutine xi_DIR_eta_PER_init(Nx,Ny,spline_degree,api_ID)
    implicit none
    !LOCAL VARIABLES
    integer :: li_d,Nx,Ny,spline_degree
    integer :: li_i, li_j, li_A
    integer :: li_dof,ai_sizePB
    integer, dimension(:) :: api_ID

    li_d = 0

    do li_j = 1, Ny

        do li_i = 1, Nx

            li_A = li_i + Nx*(li_j-1)

            if ((li_j == 1) .OR. (li_j == Ny)) then

                api_ID(li_A) = 0

            else

                if (li_i /= Nx) then

                    if (api_ID(li_A) == 0) then

                        li_d = li_d + 1

                        api_ID(li_A) = li_d

                    end if

                    if ( (1 <= li_i ) .AND. ( li_i <= spline_degree ) ) then

                        api_ID(li_A + Nx - spline_degree) = api_ID(li_A)

                    end if

                end if

            end if

        end do

    end do


end subroutine xi_DIR_eta_PER_init

!-----------------------------------------------------------------------------------

subroutine xi_DIR_eta_DIR_init(Nx,Ny,spline_degree,api_ID)
    implicit none
    !LOCAL VARIABLES
    integer :: li_d,Nx,Ny,spline_degree
    integer :: li_i, li_j, li_A
    integer :: li_dof,ai_sizePB
    integer, dimension(:) :: api_ID
    
    
    li_d = 0

    do li_j = 2, Ny-1

        do li_i = 2, Nx-1

            li_A = li_i + Nx*(li_j-1)

            li_d = li_d + 1

            api_ID(li_A) = li_d
                

        end do

    end do




end subroutine xi_DIR_eta_DIR_init

!----------------------------------------------------------------------------------------------
subroutine initLM(Nx,Ny,spline_degree,api_IEN,api_ID,api_LM)
    implicit none
    !LOCAL VARIABLES
    integer :: li_e, li_b,Nx,Ny,spline_degree
    integer, dimension(:) :: api_ID
    integer, dimension(:,:) :: api_IEN
    integer,dimension(:,:) :: api_LM
    

    do li_e = 1, Nx*Ny

        do li_b = 1, (spline_degree+1)**2
            
           ! print*, api_ID(api_IEN(li_b, li_e))
            api_LM(li_b, li_e) = api_ID(api_IEN(li_b, li_e))

        end do

    end do

end subroutine initLM
      
end Module connectivity_Module
