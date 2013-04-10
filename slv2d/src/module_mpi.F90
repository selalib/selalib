module module_mpi 
!!$ =====================================
!!$    
!!$    File:          module_mpi
!!$    Project:       vlasov
!!$    Author(s):     
!!$    Creation:      01.05.1998
!!$    Last modified: 05.03.1999
!!$    
!!$ =============================================================
!!$ =============================================================
!
! Les variables globales 
!     Nombre_Processeurs
!     Numero_Processeur
!     Root
!
! Les procédures
!     initialise_moduleMPI   : initialise les variables MPI
!     termine_moduleMPI      : arrete proprement le programme
!     erreur_mpi             : arrete proprement le programme
!
!
  use used_precision 
!  implicit none
!
use mpi
integer, parameter :: mpi_realtype = mpi_double_precision

integer :: Nombre_Processeurs, Numero_Processeur, Root
contains
!
! initialise les variables MPI
!
  subroutine initialise_moduleMPI
    integer info,nstatus
!
    call mpi_init(info)
    if(info/=0) call erreur_mpi(info)
    nstatus=mpi_status_size
    call mpi_comm_rank(mpi_comm_world, Numero_Processeur, info)
    if(info/=0) call erreur_mpi(info)
    call mpi_comm_size(mpi_comm_world, Nombre_Processeurs, info)
    if(info/=0) call erreur_mpi(info)
    Root = 0
!    write(*,*) 'Nombre de processeurs : ',Nombre_Processeurs
!    Write(*,*) 'Numero du processeurs : ',Numero_Processeur

  end subroutine initialise_moduleMPI
!
! termine le programme mpi
!
  subroutine termine_moduleMPI
    integer info
!  
    call mpi_finalize(info)       
    if(info/=0) call erreur_mpi(info)
!
  end subroutine termine_moduleMPI
!
! Pour gerer les erreurs MPI
!
  subroutine erreur_mpi(info)
    integer, intent(in) :: info
!
    write(*,*) 'Erreur MPI numero : ',info
    write(*,*) 'Arret du programme.'
    call  termine_moduleMPI()
    stop
  end subroutine erreur_mpi
end module module_mpi

MODULE Clock

  IMPLICIT NONE
  INTEGER :: rate

CONTAINS

  SUBROUTINE clck_temps(P_temps)
    INTEGER, INTENT(OUT) :: P_temps
    call system_clock(COUNT=P_temps)
  END SUBROUTINE clck_temps

  SUBROUTINE clck_diff(P_temps1,P_temps2,P_diffTemps)
    INTEGER, INTENT(IN) :: P_temps1, P_temps2
    REAL*8, INTENT(OUT) :: P_diffTemps
    call system_clock(COUNT_RATE=rate)
    P_diffTemps=float(P_temps2-P_temps1)/rate
  END SUBROUTINE clck_diff

END MODULE Clock

