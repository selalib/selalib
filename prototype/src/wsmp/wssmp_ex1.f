C------------------------------------------------------------------------------
C       Example program to show the use of the "wssmp" routine.
C------------------------------------------------------------------------------
C.. This program can be obtained from:
C
C   http://www.cs.umn.edu/~agupta/wsmp
C
C   (C) IBM Corporation, 1999, 2013.
C
c	Acceptance and use of this program constitutes the user's understanding
c	that he/she will have no recourse to IBM for any actual or consequential
c	damages, including, but not limited to, lost profits or savings, arising
c	out of the use or inability to use this program. 
C------------------------------------------------------------------------------

        program wssmp_ex1
        implicit none

        integer n, ldb, nrhs, naux
        integer iparm(64)
        integer ia(10)
        integer ja(29)
        integer perm(9), invp(9)
        double precision dparm(64)
        double precision avals(29)
        double precision b(9)

        integer mrp                     ! just a placeholder in this program
        double precision aux, diag      ! just placeholders in this program

        integer i
        double precision waltime, wsmprtc

C.. Fill all arrays containing matrix data.

        data n /9/, ldb /9/, nrhs /1/, naux /0/

        data ia /1,5,9,13,17,21,25,27,29,30/

        data ja 
     1        /1,          3,                      7,    8,
     2	             2,    3,                            8,    9,
     3                     3,                      7,    8,    9,
     4                           4,          6,    7,    8,
     5                                 5,    6,          8,    9,
     6                                       6,    7,    8,    9,
     7                                             7,    8,
     8                                                   8,    9,
     9                                                         9/
        data avals 
     1      /14.d0,      -1.d0,                  -1.d0,-3.d0,
     2             14.d0,-1.d0,                        -3.d0,-1.d0,
     3                   16.d0,                  -2.d0,-4.d0,-2.d0,
     4                         14.d0,      -1.d0,-1.d0,-3.d0,
     5                               14.d0,-1.d0,      -3.d0,-1.d0,
     6                                     16.d0,-2.d0,-4.d0,-2.d0,
     7                                           16.d0,-4.d0,
     8                                                 71.d0,-4.d0,
     9                                                       16.d0/

C.. Executable part of the program starts here.

        call wsetmaxthrds(2)    ! Use 2 threads

        waltime = wsmprtc()

C.. Fill 'iparm' and 'dparm' arrays with default values.

C.. As an alternative to this step, the values in 'iparm' and 'dparm' can be
c   filled with values suitable for the application either manually or by 
c   using a "data" statement according to the description in the User's guide.

	call wsmp_initialize
        iparm(1) = 0
        iparm(2) = 0
        iparm(3) = 0
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Initialization complete in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

C.. Ordering.

        waltime = wsmprtc()
        iparm(2) = 1
        iparm(3) = 1
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Ordering complete in time - ',wsmprtc()-waltime

        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

C.. Symbolic Factorization.

        waltime = wsmprtc()
        iparm(2) = 2
        iparm(3) = 2
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Symbolic complete in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

        print *,'Number of nonzeros in factor L = 1000 X ',iparm(24)
        print *,'Number of FLOPS in factorization = ',dparm(23)
        print *,'Double words needed to factor = 1000 X ',iparm(23) 

C.. Cholesky Factorization.

        waltime = wsmprtc()
        iparm(2) = 3
        iparm(3) = 3
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        waltime = wsmprtc() - waltime
        print *,'Cholesky complete in time - ',waltime
        print *,'Factorization MegaFlops = ',
     +           (dparm(23) * 1.d-6) / waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)

C.... If Cholesky factorization failed due to non-positive-definite
c     matrix, then try LDL' factorization.

          if (iparm(64) .gt. 0) go to 1000
          stop
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do

        waltime = wsmprtc()
        iparm(2) = 4
        iparm(3) = 4
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Back substitution done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

C.. Iterative refinement.

        waltime = wsmprtc()
        iparm(2) = 5
        iparm(3) = 5
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Iterative ref. done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

        print *,'The solution of the system is as follows:'
        do i = 1, n
          print *,i,' : ',b(i)
        end do

C.. Solve the same system using LDL' factorization.

1000	continue

        waltime = wsmprtc()
        iparm(2) = 3
        iparm(3) = 3
        iparm(31) = 1
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        waltime = wsmprtc()-waltime
        print *,'LDL^T factorization complete in time - ',waltime
        print *,'Factorization MegaFlops = ',
     +           (dparm(23) * 1.d-6) / waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do

        waltime = wsmprtc()
        iparm(2) = 4
        iparm(3) = 4
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Back substitution done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if

C.. Iterative refinement.

        waltime = wsmprtc()
        iparm(2) = 5
        iparm(3) = 5
        call wssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +              aux, naux, mrp, iparm, dparm)

        print *,'Iterative ref. done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          stop
        end if
        print *,'The solution of the system is as follows:'
        do i = 1, n
          print *,i,' : ',b(i)
        end do

        stop
        end

