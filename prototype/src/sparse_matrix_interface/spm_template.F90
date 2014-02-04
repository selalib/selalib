MODULE SPM_MODULE
    IMPLICIT NONE 

    PRIVATE 


CONTAINS

! .............................................
SUBROUTINE SPM_INITIALIZE(ai_nMatrices, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_INITIALIZE
! .............................................

! .............................................
SUBROUTINE SPM_SetDefaultOptions(, ierr)
!(id, zero, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_SetDefaultOptions
! .............................................

! .............................................
SUBROUTINE SPM_SetOptionINT(, ierr)
!(id, IPARM_VERBOSE, API_VERBOSE_NOT, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_SetOptionINT
! .............................................

! .............................................
SUBROUTINE SPM_SetOptionREAL(, ierr)
!(id, MURGE_RPARAM_EPSILON_ERROR, PREC, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_SetOptionREAL
! .............................................

! .............................................
SUBROUTINE SPM_GRAPHBEGIN(, ierr)
!(id, n, edgenbr, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_GRAPHBEGIN
! .............................................

! .............................................
SUBROUTINE SPM_GRAPHEDGE(, ierr)
!(id, i, j, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_GRAPHEDGE
! .............................................

! .............................................
SUBROUTINE SPM_GRAPHEND(, ierr)
!(id, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_GRAPHEND
! .............................................

! .............................................
SUBROUTINE SPM_CLEAN(id, ierr)
!(id, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_CLEAN
! .............................................

! .............................................
SUBROUTINE SPM_CLEANALL(ierr)
!(ierr)
IMPLICIT NONE

END SUBROUTINE SPM_CLEANALL
! .............................................

! .............................................
SUBROUTINE SPM_FINALIZE(ierr)
!(ierr)
IMPLICIT NONE

END SUBROUTINE SPM_FINALIZE
! .............................................

! .............................................
SUBROUTINE SPM_MATRIXRESET(id, ierr)
!(id, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_MATRIXRESET
! .............................................

! .............................................
SUBROUTINE SPM_ASSEMBLYBEGIN()
!(id, ncount, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_FOOL, MURGE_BOOLEAN_FALSE, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_ASSEMBLYBEGIN
! .............................................

! .............................................
SUBROUTINE SPM_ASSEMBLYEND(id, ierr)
!(id, ierr)
IMPLICIT NONE

END SUBROUTINE SPM_ASSEMBLYEND
! .............................................

END MODULE SPM_MODULE
