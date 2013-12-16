C------------------------------------------------------------------------------
C       Example parallel program to show the use of the "pwssmp" routine 
c       in the peer mode. Note: Use at least 3 MPI processes.
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

        program pwssmp_ex1
        implicit none
        include 'mpif.h'
        integer mperr, myid

        call mpi_init(mperr)
        if (mperr .ne. 0) then
          print *,'MPI initialization error'
          stop
        end if
        call mpi_comm_rank(MPI_COMM_WORLD,myid,mperr)

C..  This program must be called on *at least* 3 nodes.

        call wsetmaxthrds(2)    ! Use 2 threads on each process

        if (myid .eq. 0) then
          call master0
        else if (myid .eq. 1) then
          call master1
        else if (myid .eq. 2) then
          call master2
        else
          call slave (myid)
        end if

        call mpi_finalize(mperr)
        stop
        end
C------------------------------------------------------------------------------
        subroutine master0 
        implicit none

        integer n, ldb, nrhs, naux
        integer iparm(64)
        integer ia(4)
        integer ja(12)
        integer perm(9), invp(9)
        double precision dparm(64)
        double precision avals(12)
        double precision b(3)

        integer mrp                     ! just a placeholder in this program
        double precision aux, diag      ! just placeholders in this program

        integer i
        double precision waltime, wsmprtc

C.. Fill all arrays containing matrix data.

        data n /3/, ldb /3/, nrhs /1/, naux /0/

        data ia /1,5,9,13/

        data ja 
     1        /1,          3,                      7,    8,
     2	             2,    3,                            8,    9,
     3                     3,                      7,    8,    9/
        data avals 
     1      /14.d0,      -1.d0,                  -1.d0,-3.d0,
     2             14.d0,-1.d0,                        -3.d0,-1.d0,
     3                   16.d0,                  -2.d0,-4.d0,-2.d0/

C.. Executable part of the program starts here.

        waltime = wsmprtc()

C.. Fill 'iparm' and 'dparm' arrays with default values.

C.. As an alternative to this step, the values in 'iparm' and 'dparm' can be
c   filled with values suitable for the application either manually or by 
c   using a "data" statement according to the description in the User's guide.

        iparm(1) = 0
        iparm(2) = 0
        iparm(3) = 0
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        print *,'Initialization complete in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
        end if

C.. Ordering.

        waltime = wsmprtc()
        iparm(2) = 1
        iparm(3) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        print *,'Ordering complete in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
        end if

C.. Symbolic Factorization.

        waltime = wsmprtc()
        iparm(2) = 2
        iparm(3) = 2
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        print *,'Symbolic complete in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
        end if
        print *,'Number of nonzeros in factor L = 1000 X ',iparm(24)
        print *,'Number of FLOPS in factorization = ',dparm(23)
        print *,'Double words needed to factor on 0 = 1000 X ',
     +		 iparm(23) 

C.. Cholesky Factorization.

        waltime = wsmprtc()
        iparm(2) = 3
        iparm(3) = 3
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        waltime = wsmprtc() - waltime
        print *,'Cholesky complete in time - ',waltime
        print *,'Factorization MegaFlops = ',
     +           (dparm(23) * 1.d-6) / waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)

C.... If Cholesky factorization failed due to non-positive-definite
c     matrix, then try LDL' factorization.

          if (iparm(64) .gt. 0) go to 1000
          return
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        waltime = wsmprtc()
        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        print *,'Back substitution done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
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
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        waltime = wsmprtc()-waltime
        print *,'LDL^T factorization complete in time - ',waltime
        print *,'Factorization MegaFlops = ',
     +           (dparm(23) * 1.d-6) / waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        waltime = wsmprtc()
        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        print *,'Back substitution done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
        end if

C.. Iterative refinement.

        waltime = wsmprtc()
        iparm(2) = 5
        iparm(3) = 5
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        print *,'Iterative ref. done in time - ',wsmprtc()-waltime
        if (iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',iparm(64)
          return
        end if
        print *,'The solution of the system is as follows:'
        do i = 1, n
          print *,i,' : ',b(i)
        end do

        return
        end
C------------------------------------------------------------------------------
        subroutine master1
        implicit none

        integer n, ldb, nrhs, naux
        integer iparm(64)
        integer ia(3)
        integer ja(8)
        integer perm(9), invp(9)
        double precision dparm(64)
        double precision avals(8)
        double precision b(2)

        integer mrp                     ! just a placeholder in this program
        double precision aux, diag      ! just placeholders in this program

        integer i
        double precision waltime, wsmprtc

C.. Fill all arrays containing matrix data.

        data n /2/, ldb /2/, nrhs /1/, naux /0/

        data ia /1,5,9/

        data ja
     4      /                    4,          6,    7,    8,
     5                                 5,    6,          8,    9/
        data avals
     4      /                  14.d0,      -1.d0,-1.d0,-3.d0,
     5                               14.d0,-1.d0,      -3.d0,-1.d0/

C.. Fill 'iparm' and 'dparm' arrays with default values.

        iparm(1) = 0
        iparm(2) = 0
        iparm(3) = 0
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Ordering.

        iparm(2) = 1
        iparm(3) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Symbolic Factorization.

        iparm(2) = 2
        iparm(3) = 2
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

        print *,'Double words needed to factor on 1 = 1000 X ',
     +		 iparm(23)

C.. Cholesky Factorization.

        iparm(2) = 3
        iparm(3) = 3
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        if (iparm(64) .ne. 0) then
          if (iparm(64) .gt. 0) go to 1000
          return
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

        do i = 1, n
          print *,i+3,' : ',b(i)
        end do

C.. Solve the same system using LDL' factorization.

1000    continue

        iparm(2) = 3
        iparm(3) = 3
        iparm(31) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return 

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Iterative refinement.

        iparm(2) = 5
        iparm(3) = 5
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

        do i = 1, n
          print *,i+3,' : ',b(i)
        end do

        return
        end
C------------------------------------------------------------------------------
        subroutine master2
        implicit none

        integer n, ldb, nrhs, naux
        integer iparm(64)
        integer ia(5)
        integer ja(9)
        integer perm(9), invp(9)
        double precision dparm(64)
        double precision avals(9)
        double precision b(4)

        integer mrp                     ! just a placeholder in this program
        double precision aux, diag      ! just placeholders in this program

        integer i
        double precision waltime, wsmprtc

C.. Fill all arrays containing matrix data.

        data n /4/, ldb /4/, nrhs /1/, naux /0/

        data ia /1,5,7,9,10/

        data ja
     6      /                                6,    7,    8,    9,
     7                                             7,    8,
     8                                                   8,    9,
     9                                                         9/
        data avals
     6      /                              16.d0,-2.d0,-4.d0,-2.d0,
     7                                           16.d0,-4.d0,
     8                                                 71.d0,-4.d0,
     9                                                       16.d0/

C.. Fill 'iparm' and 'dparm' arrays with default values.

        iparm(1) = 0
        iparm(2) = 0
        iparm(3) = 0
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Ordering.

        iparm(2) = 1
        iparm(3) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Symbolic Factorization.

        iparm(2) = 2
        iparm(3) = 2
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

        print *,'Double words needed to factor on 2 = 1000 X ',
     +		 iparm(23)

C.. Cholesky Factorization.

        iparm(2) = 3
        iparm(3) = 3
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        if (iparm(64) .ne. 0) then
          if (iparm(64) .gt. 0) go to 1000
          return
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

        do i = 1, n
          print *,i+5,' : ',b(i)
        end do

C.. Solve the same system using LDL' factorization.

1000    continue

        iparm(2) = 3
        iparm(3) = 3
        iparm(31) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Iterative refinement.

        iparm(2) = 5
        iparm(3) = 5
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

        do i = 1, n
          print *,i+5,' : ',b(i)
        end do

        return
        end
C------------------------------------------------------------------------------
        subroutine slave (myid)
        implicit none
        integer myid

        integer n, ldb, nrhs, naux
        integer iparm(64)
        integer ia, ja 
        integer perm(9), invp(9)
        double precision dparm(64)
        double precision avals, b       ! just placeholders in this routine

        integer mrp                     ! just a placeholder in this program
        double precision aux, diag      ! just placeholders in this program

        data n /0/, ldb /1/, nrhs /1/, naux /0/

C.. Fill 'iparm' and 'dparm' arrays with default values.

        iparm(1) = 0
        iparm(2) = 0
        iparm(3) = 0
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Ordering.

        iparm(2) = 1
        iparm(3) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Symbolic Factorization.

        iparm(2) = 2
        iparm(3) = 2
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        print *,'Double words needed to factor on',myid,'= 1000 X ',
     +		 iparm(23) 
        if (iparm(64) .ne. 0) return

C.. Cholesky Factorization.

        iparm(2) = 3
        iparm(3) = 3
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) then
          if (iparm(64) .gt. 0) go to 1000
          return
        end if

C.. Back substitution.

        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Iterative refinement.

        iparm(2) = 5
        iparm(3) = 5
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Solve the same system using LDL' factorization.

1000    continue

        iparm(2) = 3
        iparm(3) = 3
        iparm(31) = 1
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Back substitution.

        iparm(2) = 4
        iparm(3) = 4
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)
        if (iparm(64) .ne. 0) return

C.. Iterative refinement.

        iparm(2) = 5
        iparm(3) = 5
        call pwssmp (n, ia, ja, avals, diag, perm, invp, b, ldb, nrhs,
     +               aux, naux, mrp, iparm, dparm)

        return
        end

