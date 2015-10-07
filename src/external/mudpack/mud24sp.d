c
c     file mud24sp.d
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 5.0                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c     For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c
c
c ... file mud24sp.d
c
C     contains documentation for subroutine mud24sp(work,phi,ierror)
c     A sample fortran driver is file "tmud24sp.f".
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213
c
c ... required MUDPACK files
c
c     mud2sp.f, mudcom.f
c
c ... purpose
c
c     mud24sp attempts to improve the estimate in phi, obtained by calling
c     mud2sp,  from second to fourth order accuracy.  see the file "mud2sp.d"
c     for a detailed discussion of the elliptic pde approximated and
c     arguments "work,phi" which are also part of the argument list for
c     mud2sp.
c
c ... assumptions
c
c     *  phi contains a second-order approximation from an earlier mud2sp call
c
c     *  arguments "work,phi" are the same used in calling mud2sp
c
c     *  "work,phi" have not changed since the last call to mud2sp
c
c     *  the finest grid level contains at least 6 points in each direction
c
c
c *** warning
c
c     if the first assumption is not true then a fourth order approximation
c     cannot be produced in phi.  the last assumption (adequate grid size)
c     is the only one checked. failure in any of the others can result in
c     in an undetectable error.
c
c ... language                                                                  
c
c     fortran90/fortran77
c
c ... portability                                                               
c
c     mudpack5.0 software has been compiled and tested with fortran77
c     and fortran90 on a variety of platforms.
c                                                                               
c ... methods
c
c     details of the methods employeed by the solvers in mudpack are given
c     in [1,9].  [1,2,9] contain performance measurements on a variety of
c     elliptic pdes (see "references" in the file "readme").  in summary:
c
c
c *** higher order solution (fourth-order solvers) (see [9,19,21])
c
c     if the multigrid cycling results in a second-order estimate (i.e.,
c     discretization level error is reached) then this can be improved to a
c     fourth-order estimate using the technique of "deferred corrections"
c     the values in the solution array are used to generate a fourth-order
c     approximation to the truncation error.  second-order finite difference
c     formula are used to approximate third and fourth partial derivatives
c     of the solution function [3].  the truncation error estimate is
c     transferred down to all grid levels using weighted averaging where
c     it serves as a new right hand side.  the default multigrid options
c     are used to compute the fourth-order correction term which is added
c     to the original solution array.
c
c
c ... references (partial)
c
c     for a complete list see "references" in the mudpack information and
c     directory file "readme" or visit the mudpack web site
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c     [1] J. Adams, "MUDPACK: Multigrid Fortran Software for the Efficient
c     Solution of Linear Elliptic Partial Differential Equations,"
c     Applied Math. and Comput. vol.34, Nov 1989, pp.113-146.
c
c     [2] J. Adams,"FMG Results with the Multigrid Software Package MUDPACK,"
c     proceedings of the fourth Copper Mountain Conference on Multigrid, SIAM,
c     1989, pp.1-12.
c     .
c     .
c     .
c     [7] J. Adams, R. Garcia, B. Gross, J. Hack, D. Haidvogel, and V. Pizzo,
c     "Applications of Multigrid Software in the Atmospheric Sciences,"
c      Mon. Wea. Rev.,vol. 120 # 7, July 1992, pp. 1447-1458.
c     .
c     .
c     .
c     [9] J. Adams, "Recent Enhancements in MUDPACK, a Multigrid Software
c     package for Elliptic Partial Differential Equations," Applied Math.
c     and Comp., 1991, vol. 43, May 1991, pp. 79-94.
c
c     [10]J. Adams, "MUDPACK-2: Multigrid Software for Approximating
c     Elliptic Partial Differential Equations on Uniform Grids with
c     any Resolution," Applied Math. and Comp., 1993, vol. 53, February
c     1993, pp. 235-249
c     .
c     .
c     .
c
c ... error argument
c
c     = 0 if no error is detected
c
c     = 30 if min0(nx,ny) < 6 where nx,ny are the fine grid sizes
c          in the x,y directions.
c                                                                               
c
c  ***********************************************************************
c  ***********************************************************************
c
c     end of mud24sp documentation
c
c  ***********************************************************************
c  ***********************************************************************
c
