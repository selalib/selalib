!     
! File:   connectivity.F90
! Author: root
!
! Created on December 20, 2011, 3:36 PM
!

module connectivity_module
    use tracelog_module
    use connectivities_def
    implicit none

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !---------------------------------------------------------------
    subroutine create_connectivity(self, ai_npatch, ai_maxnen, ai_nelt, ai_N &
        , api_IEN, api_ID, api_LM, api_nen)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_npatch
        integer :: ai_maxnen
        integer :: ai_nelt
        integer :: ai_N
        integer, dimension(:,:,:), optional :: api_IEN
        integer, dimension(:), optional :: api_ID
        integer, dimension(:,:,:), optional :: api_LM
        integer, dimension(:,:,:), optional :: api_nen
        ! LOCAL

        CALL printlog("create_connectivity : Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_IEN(ai_npatch, ai_maxnen, ai_nelt))
        ALLOCATE(self % opi_ID(ai_N))
        ALLOCATE(self % opi_LM(ai_npatch, ai_maxnen, ai_nelt))
        ALLOCATE(self % opi_nen(ai_npatch, ai_nelt))

        ! opi_nen is always initialized as : opi_nen = maxnen
        self % opi_nen = ai_maxnen

#ifdef _DEBUG
call concatmsg("ai_npatch = ", ai_dtllevel = mi_dtllevel_base + 1)
call concatmsg(ai_npatch, ai_dtllevel = mi_dtllevel_base + 1)

call concatmsg("ai_nelt = ", ai_dtllevel = mi_dtllevel_base + 1)
call concatmsg(ai_nelt, ai_dtllevel = mi_dtllevel_base + 1)

call concatmsg("ai_maxnen = ", ai_dtllevel = mi_dtllevel_base + 1)
call concatmsg(ai_maxnen, ai_dtllevel = mi_dtllevel_base + 1)

call printmsg(ai_dtllevel = mi_dtllevel_base + 1)

call concatmsg("self % opi_nen = ", ai_dtllevel = mi_dtllevel_base + 1)
call concatmsg(self % opi_nen, ai_dtllevel = mi_dtllevel_base + 1)

call printmsg(ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (present(api_IEN)) then
            self % opi_IEN = api_IEN
        end if

        if (present(api_ID)) then
            self % opi_ID = api_ID
        end if

        if (present(api_LM)) then
            self % opi_LM = api_LM
        end if

        if (present(api_nen)) then
            print *, 'create_connectivity : Not yet implemented for api_nen'
            self % opi_nen = ai_maxnen
!            self % opi_nen = api_nen
        end if

        CALL printlog("create_connectivity : End", ai_dtllevel = 1)
    end subroutine create_connectivity
    !---------------------------------------------------------------
    subroutine create_real_elts(self, ai_npatch, ai_nelt, api_real_elts)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_npatch
        integer :: ai_nelt
        integer, dimension(:,:), optional :: api_real_elts
        ! LOCAL

        CALL printlog("create_real_elts : Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_real_elts(0:ai_npatch-1, ai_nelt))

        if (present(api_real_elts)) then
            self % opi_real_elts = api_real_elts
        end if

        CALL printlog("create_real_elts : End", ai_dtllevel = 1)
    end subroutine create_real_elts
    !---------------------------------------------------------------
    subroutine free_connectivity(self)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self

        CALL printlog("free_connectivity : Start", ai_dtllevel = 1)

        DEALLOCATE(self % opi_IEN)
        DEALLOCATE(self % opi_ID)
        DEALLOCATE(self % opi_LM)
        DEALLOCATE(self % opi_nen)
        DEALLOCATE(self % opi_real_elts)
        
        CALL printlog("free_connectivity : End", ai_dtllevel = 1)
    end subroutine free_connectivity


end module connectivity_module
