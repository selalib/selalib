C------------------------------------------------------------------------------
C       Example program to show the use of the simple interface for
c       solving symmetric sparse systems in the WSMP library.
C------------------------------------------------------------------------------
C.. This program can be obtained from:
C
C   http://www.cs.umn.edu/~agupta/wsmp
C
C   (C) IBM Corporation, 1999.
C
c	Acceptance and use of this program constitutes the user's understanding
c	that he/she will have no recourse to IBM for any actual or consequential
c	damages, including, but not limited to, lost profits or savings, arising
c	out of the use or inability to use this program. 
C------------------------------------------------------------------------------

        program wssmp_ex2
        implicit none

        integer n, ldb, nrhs, naux, nnzl, wspace, niter, info
        integer options(5)
        integer ia(10)
        integer ja(29)
        integer perm(9), invp(9)
        real*8 avals(29)
        real*8 b(9)
        real*8 rcond
        integer aux            ! just a placeholder in this program

        integer i
        real*8 waltime, wsmprtc

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

        data options /5 * 0/

C.. Ordering and Symbolic.

        call wsetmaxthrds(2)    ! Use 2 threads

        waltime = wsmprtc()
        call wscalz (n, ia, ja, options, perm, invp, nnzl, wspace, 
     +               aux, naux, info)

        print *,'Analyze step complete in time - ',wsmprtc()-waltime
        if (info .ne. 0) then
          print *,'The following ERROR was detected: ',info
          stop
        end if
        print *,'Number of nonzeros in factor L = ',nnzl
        print *,'Double words needed to factor = 1000 X ',
     +		 wspace

C.. Cholesky Factorization.

        waltime = wsmprtc()
        call wscchf (n, ia, ja, avals, perm, invp, aux, naux, info)

        print *,'Cholesky complete in time - ', wsmprtc() - waltime
        if (info .ne. 0) then
          print *,'The following ERROR was detected: ',info

C.... If Cholesky factorization failed due to non-positive-definite
c     matrix, then try LDL' factorization.

          if (info .gt. 0) go to 1000
          stop
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        waltime = wsmprtc()
        niter = 1
        call wsslv (n, perm, invp, b, ldb, nrhs, niter, aux, naux)

        print *,'Back substitution done in time - ',wsmprtc()-waltime
        print *,'The solution of the system is as follows:'
        do i = 1, n
          print *,i,' : ',b(i)
        end do

C.. Solve the same system using LDL' factorization.

1000	continue

        waltime = wsmprtc()
        call wscldl (n, ia, ja, avals, perm, invp, aux, naux, 
     +               info)      

        print *,'LDL^T factorization done in time - ',wsmprtc()-waltime
        if (info .ne. 0) then
          print *,'The following ERROR was detected: ',info
          stop
        end if

C.. Back substitution.

        do i = 1, n
          b(i) = 1.d0
        end do
        waltime = wsmprtc()
        niter = 1
        call wsslv (n, perm, invp, b, ldb, nrhs, niter, aux, naux)

        print *,'Back substitution done in time - ',wsmprtc()-waltime
        print *,'The solution of the system is as follows:'
        do i = 1, n
          print *,i,' : ',b(i)
        end do

C.. Solve the same system using the expert driver.

        waltime = wsmprtc()
        do i = 1, n
          b(i) = 1.d0
        end do
        call wscsvx (n, ia, ja, avals, perm, invp, b, ldb, nrhs,
     +               aux, naux, rcond, info)

        print *,'Run time for expert driver = ',wsmprtc() - waltime
        if (info .ne. 0) then
          print *,'The following ERROR was detected: ',info
          stop
        end if
        print *,'Condition number estimate of matrix = ',1.d0/rcond
        print *,'The solution of the system is as follows:'
        do i = 1, n
          print *,i,' : ',b(i)
        end do

        stop
        end

