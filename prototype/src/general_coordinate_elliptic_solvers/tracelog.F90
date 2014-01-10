Module tracelog_Module
use used_precision
Implicit None

    !> ml_stdoutput, IF THE OUTPUT IS THE SCREEN, DEFAULT VALUE IS 0
    logical :: ml_stdoutput
    !> DETAIL LEVEL, DEFAULT VALUE IS 0
    integer :: mi_dtllevel
!    character(LEN=len("/Runs/log/trace.log")), parameter	:: ms_fileoutput = "/Runs/log/trace.log"
    integer, parameter :: mi_traceflux = 1000
#ifdef _DEBUG
    integer, parameter :: MAXMSG = 80000
#else
    integer, parameter :: MAXMSG = 1000
#endif
    character(len=MAXMSG) :: ms_msg

    CHARACTER(LEN=10), PARAMETER :: FMT1_REAL = "(F8.3)"

    integer, parameter :: MAXSTAMP = 16
    character(len=MAXSTAMP) :: ms_stamp


    !> CURRENT TIME CPU
    real    ( kind = 8 ) mr_time

    interface concatmsg
            module procedure concatmsg_string           &
                   , concatmsg_int              &
                   , concatmsg_realwp           &
                   , concatmsg_array1realwp     &
                   , concatmsg_array2realwp     &
                   , concatmsg_array3realwp     &
                   , concatmsg_array4realwp     &
                   , concatmsg_array1int        &
                   , concatmsg_array2int        &
                   , concatmsg_array1complexwp
    end interface

contains
!-------------------------------------------------------------------------------------------
    subroutine opentracelog ( al_stdoutput, ai_dtllevel )
    implicit none
        !> DETAIL LEVEL
        logical, optional :: al_stdoutput
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        integer  :: li_ios

        mi_dtllevel  = 0
        if ( present ( ai_dtllevel )  ) then
            mi_dtllevel = ai_dtllevel
        end if

        ml_stdoutput = .TRUE.
        if ( present ( al_stdoutput )  ) then
            ml_stdoutput = al_stdoutput
        end if

        ms_msg = ""
        ms_stamp = ""

!        call printcputime ( )

        ! IF USING FILE FOR OUTPUT
        ! WE MUST OPEN THE LOG FILE
        if ( .not. ml_stdoutput ) then

            open(unit=mi_traceflux, file='trace.log',action="write", iostat=li_ios)
            if (li_ios/=0) STOP "erreur d'ouverture du fichier du fichier trace.log"

        end if

    end subroutine opentracelog
!-------------------------------------------------------------------------------------------
    subroutine closetracelog (  )
    implicit none

!        call printcputime ( )

!        call concatmsg_string ( "Total CPU time :" )
!
!		call cpu_time ( mr_time )
!
!        call concatmsg_realwp ( mr_time )
!        call printmsg ( ai_dtllevel = 0 )

        ! IF USING FILE FOR OUTPUT
        ! WE MUST CLOSE THE LOG FILE
        if ( .not. ml_stdoutput ) then

            close(mi_traceflux)

        end if

    end subroutine closetracelog
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE PRINTS THE MESSAGE as_message IF THE al_condition IS TRUE
    subroutine printlog ( as_message, al_condition, ai_dtllevel )
    implicit none
        !> DETAIL LEVEL
        character(len=*), intent(in) :: as_message
        logical, optional :: al_condition
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        logical :: ll_condition
        integer :: li_dtllevel

        ll_condition = .TRUE.
        if ( present ( al_condition ) ) then
            ll_condition = al_condition
        end if

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! IF al_condition IS NOT TRUE, THEN EXIT
        if ( .not. ll_condition ) then
            return
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

        call setstamp ( li_dtllevel )

        ! STANDARD OUTPUT
        if ( ml_stdoutput ) then

            print*,TRIM ( ms_stamp ) // as_message

        ! FILE OUTPUT
        else

            write(mi_traceflux, *)TRIM ( ms_stamp ) // as_message

        end if

        ms_stamp = ""

    end subroutine printlog
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE PRINTS THE MESSAGE as_message IF THE al_condition IS TRUE
    subroutine printmsg ( al_condition, ai_dtllevel )
    implicit none
        !> DETAIL LEVEL
        logical, optional :: al_condition
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        logical :: ll_condition
        integer :: li_dtllevel

        ll_condition = .TRUE.
        if ( present ( al_condition ) ) then
            ll_condition = al_condition
        end if

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! IF al_condition IS NOT TRUE, THEN EXIT
        if ( .not. ll_condition ) then
            ms_msg = ""
            ms_stamp = ""
            return
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            ms_msg = ""
            ms_stamp = ""
            return
        end if

        call setstamp ( li_dtllevel )

        ms_msg = TRIM ( ms_stamp ) // TRIM ( ms_msg )

        ! STANDARD OUTPUT
        if ( ml_stdoutput ) then

            print*,TRIM ( ms_msg )

        ! FILE OUTPUT
        else

            write(mi_traceflux, *)TRIM ( ms_msg )

        end if

        ! reinitilize the message
        ms_msg = ""
        ms_stamp = ""

    end subroutine printmsg
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_string ( as_message, ai_dtllevel )
    implicit none
        character(len=*), intent(in) :: as_message
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        integer :: li_dtllevel

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( as_message )

    end subroutine concatmsg_string
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_int ( ai_value, ai_dtllevel )
    implicit none
        integer  :: ai_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

        write(ls_msg,*) ai_value
        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_int
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_realwp ( ar_value, ai_dtllevel )
    implicit none
        real(wp)  :: ar_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

!        write(ls_msg,*) ar_value
        write(ls_msg,FMT1_REAL) ar_value
        
        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_realwp
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array1realwp ( apr_value, ai_dtllevel )
    implicit none
        real(wp), dimension(:)  :: apr_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel
        integer :: li_i

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

!        write(ls_msg,*) apr_value

        DO li_i = 1, SIZE(apr_value,1)
            write(ls_msg,FMT1_REAL) apr_value(li_i)
        END DO

        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array1realwp
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array2realwp ( apr_value, ai_dtllevel )
    implicit none
        real(wp), dimension(:,:)  :: apr_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel
        integer :: li_i
        integer :: li_j

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

!        write(ls_msg,*) apr_value
        
        DO li_j = 1, SIZE(apr_value,2)
            DO li_i = 1, SIZE(apr_value,1)
                write(ls_msg,FMT1_REAL) apr_value(li_i, li_j)
            END DO
        END DO

        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array2realwp
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array3realwp ( apr_value, ai_dtllevel )
    implicit none
        real(wp), dimension(:,:,:)  :: apr_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel
        integer :: li_i
        integer :: li_j
        integer :: li_k
        
        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

!        write(ls_msg,*) apr_value

        DO li_k = 1, SIZE(apr_value,3)
            DO li_j = 1, SIZE(apr_value,2)
                DO li_i = 1, SIZE(apr_value,1)
                    write(ls_msg,FMT1_REAL) apr_value(li_i, li_j, li_k)
                END DO
            END DO
        END DO

        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array3realwp
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array4realwp ( apr_value, ai_dtllevel )
    implicit none
        real(wp), dimension(:,:,:,:)  :: apr_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel
        integer :: li_i
        integer :: li_j
        integer :: li_k
        integer :: li_p
        
        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

!        write(ls_msg,*) apr_value

        DO li_p = 1, SIZE(apr_value,3)
            DO li_k = 1, SIZE(apr_value,3)
                DO li_j = 1, SIZE(apr_value,2)
                    DO li_i = 1, SIZE(apr_value,1)
                        write(ls_msg,FMT1_REAL) apr_value(li_i, li_j, li_k, li_p)
                    END DO
                END DO
            END DO
        END DO
        
        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array4realwp
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array1int ( api_value, ai_dtllevel )
    implicit none
        integer, dimension(:)  :: api_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

        write(ls_msg,*) api_value
        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array1int
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array2int ( api_value, ai_dtllevel )
    implicit none
        integer, dimension(:,:)  :: api_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

        write(ls_msg,*) api_value
        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array2int
!-------------------------------------------------------------------------------------------
    !> THIS ROUTINE COPIES as_message INTO ms_msg
    subroutine concatmsg_array1complexwp ( apc_value, ai_dtllevel )
    implicit none
        complex(wp), dimension(:)  :: apc_value
        !> DETAIL LEVEL
        integer, optional :: ai_dtllevel
        ! LOCAL VARIABLES
        character(len=MAXMSG) :: ls_msg
        ! LOCAL VARIABLES
        integer :: li_dtllevel

        li_dtllevel = mi_dtllevel
        if ( present ( ai_dtllevel ) ) then
            li_dtllevel = ai_dtllevel
        end if

        ! WE PRINT THE MESSAGE ONLY IF THE DETAIL LEVEL IS LESS THAN THE DESIRED ONE
        if ( li_dtllevel > mi_dtllevel ) then
            return
        end if

        write(ls_msg,*) apc_value
        ms_msg = TRIM ( ms_msg ) // " " // TRIM ( ADJUSTL ( ls_msg ) )

    end subroutine concatmsg_array1complexwp
!-------------------------------------------------------------------------------------------
    subroutine printcputime ( ai_detail )
    implicit none
    integer, optional, intent ( in ) :: ai_detail
        !LOCAL VARIABLES
    integer :: li_detail
    real    ( kind = 8 ) lr_time
    
    li_detail = 2
    if ( present ( ai_detail ) ) then
       li_detail = ai_detail
    end if
    
    call concatmsg_string ( "CPU time :" )
    
    call cpu_time ( lr_time )
    
    call concatmsg_realwp ( lr_time - mr_time )
    call printmsg ( ai_dtllevel = li_detail )
    
    mr_time = lr_time
    
  end subroutine printcputime
  !-----------------------------------------------------------------------
  subroutine setstamp ( ai_dtllevel )
    implicit none
    integer  :: ai_dtllevel
    !LOCAL VARIABLES
    integer  :: li_i
    
    do li_i = 1, ai_dtllevel
       
       ms_stamp = TRIM ( ADJUSTL ( ADJUSTR ( ms_stamp ) ) ) // "---"
       
    end do
    
    ms_stamp = TRIM ( ADJUSTL ( ADJUSTR ( ms_stamp ) ) ) // "> "
    
  end subroutine setstamp
  !-----------------------------------------------------------------------
  integer function print_evolution( ai_i, ai_N, ai_k, ai_q, ai_err )
    implicit none
    integer  :: ai_i, ai_N, ai_k, ai_err, ai_q
    ! LOCAL VARIABLES
    real(wp) :: lr_d
    integer  :: li_q
    
    if ( ai_err == 0 ) then
       print_evolution = 0
       return
    end if
    
    lr_d = float ( ai_N ) / float ( ai_k )
    li_q = int ( ai_i / lr_d )
    if ( li_q > ai_q ) then
       ai_q = li_q
       write(6,'(a)')'.' ! originally the format specifier was (a$)
    end if
    
    if ( ai_i == ai_N ) then
       write(6,'(a)')' done'
       !			print*,' '
    end if
    
    print_evolution = 1
    
  end function print_evolution
End Module tracelog_Module



