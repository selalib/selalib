!! Copyright (C) 2009,2010,2011,2012  Marco Restelli
!!
!! This file is part of:
!!   FEMilaro -- Finite Element Method toolkit
!!
!! FEMilaro is free software; you can redistribute it and/or modify it
!! under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 3 of the License, or
!! (at your option) any later version.
!!
!! FEMilaro is distributed in the hope that it will be useful, but
!! WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with FEMilaro; If not, see <http://www.gnu.org/licenses/>.
!!
!! author: Marco Restelli                   <marco.restelli@gmail.com>

module mod_sparse

  !General comments: provides minimal support for working with matrices
  !in compact column format, copact row format and triplet format.
  !
  !Notice:
  !1) Although it is technically possible, it is strongly discouraged
  ! modifying the integer fields of the sparse types: n, m and nz. It is
  ! also discouraged changing the allocation status of the vector
  ! fields. For these operations, users should rely on the subroutines
  ! provided in this module. On the other hand, it is fine to make
  ! assignements at the elements of the allocatable components, since
  ! this can not destroy the internal consistency of the object.
  !
  !2) When refering to the present module, all indexes start from 0. In
  ! particular, indexes in t_col and t_tri start from 0.
  !
  !The column storage is such that the matrix
  !       
  !      [2  3  0 0 0]
  !      [3  0  4 0 6]
  !  a = [0 -1 -3 2 0]
  !      [0  0  1 0 0]
  !      [0  4  2 0 1]
  !
  !is stored as
  !
  !type t_col
  !  integer :: ierr = 0
  !  integer :: n = 5
  !  integer :: m = 5
  !  integer :: nz = 12
  !  integer, allocatable  :: ap(:)={0, 2, 5, 9, 10, 12}
  !  integer, allocatable  :: ai(:)={0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4}
  !  sll_real64, allocatable :: ax(:)={2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1}
  !end type t_col
  !
  !See /usr/share/doc/umfpack-5.2.0/UserGuide.pdf for more details.
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !here are adaptations for selalib
#include "sll_working_precision.h"
  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

  ! Module interface

  public :: &
       mod_sparse_constructor, &
       mod_sparse_destructor,  &
       mod_sparse_initialized, &
                                ! sparse types
       t_intar,     &
       t_col,       &
       t_tri,       &
       t_rp,        &
       t_pm_sk,     &
                                ! construction of new objects
       new_col,     &
       new_tri,     &
                                ! convertions
       col2tri,     &
       tri2col,     &
       tri2col_skeleton, &
       tri2col_skeleton_part, &
                                ! overloaded operators
       operator(+), &
       operator(-), &
       operator(*), &
       sum,         &
       transpose,   &
       matmul,      &
                                ! error codes
       wrong_n,     &
       wrong_m,     &
       wrong_nz,    &
       wrong_dim,   &
                                ! other functions
       nnz_col,     &
       nz_col,      &
       nz_col_i,    &
       get,         &
       set,         &
       diag,        &
       spdiag,      &
                                ! deallocate
       clear

  private

  !-----------------------------------------------------------------------

  ! Module types and parameters

  ! public members
  integer, parameter :: &
       wrong_n   = 1,   &
       wrong_m   = 2,   &
       wrong_nz  = 3,   &
       wrong_rel = 4,   &
       wrong_dim = 5

  !> column oriented storage
  type t_col
     integer :: ierr = 0
     integer :: n !< rows
     integer :: m !< columns
     integer :: nz=-1 !< non zero elements
     integer, allocatable  :: ap(:) !< (m+1)
     integer, allocatable  :: ai(:) !< (nz)
     sll_real64, allocatable :: ax(:) !< (nz)
  contains
     procedure, pass(a) :: check => check_col
  end type t_col

  !> triple <tt>(i,j,x)</tt> oriented storage
  type t_tri
     integer :: ierr = 0
     integer :: n !< rows
     integer :: m !< columns
     integer :: nz !< non zero elements
     integer, allocatable  :: ti(:) !< (nz)
     integer, allocatable  :: tj(:) !< (nz)
     sll_real64, allocatable :: tx(:) !< (nz)
  end type t_tri

  !> Pointer to real, used in \c t_pm_sk.
  type t_rp
     sll_real64, pointer :: p => null()
  end type t_rp

  !> Partitioned matrix skeleton
  !!
  !! A problem which is often encountered in finite element
  !! implementations is the following:
  !! <ul>
  !!  <li> build a matrix \f$M\f$ in compact column format, given a
  !!  collection of triplets \f$\{i,j,x\}\f$ (with repetitions of the
  !!  \f$i,j\f$ indexes implying summation);
  !!  <li> partition \f$M\f$ as
  !!  \f{displaymath}{
  !!  M = \left[\begin{array}{ccc}
  !!    M(idi_1,idj_1) & \ldots & M(idi_1,idj_{n_j}) \\\
  !!     \vdots & \ddots & \vdots \\\
  !!    M(idi_{n_i},idj_1) & \ldots & M(idi_{n_i},idj_{n_j})
  !!  \end{array}\right],
  !!  \f}
  !!  where \f$idi_i\f$ and \f$idj_j\f$ are index arrays.
  !! </ul>
  !! In practice, a significant gain in the execution time can be
  !! obtained by precomputing the correspondence between the given
  !! triplets and the representation of the submatrices. This is
  !! especially true if the same matrix partitioning is required for
  !! different values of the matrix entries \f$x\f$, as it is typically
  !! the case for nonlinear problems.
  !! 
  !! This type represents the skeleton of a partitioned matrix. More in
  !! details, it shows how the \c n_in input values shall be collocated
  !! into a collection of matrices, each of which is represented as a
  !! \c t_col object. The matrices can be assembled by
  !! \code
  !!  type(t_pm_sk), target :: pm
  !!
  !!  ! ...
  !!
  !!  ! initialization to zero
  !!  do j=1,shape(pm%m,2)
  !!    do i=1,shape(pm%m,1)
  !!      pm%m(i,j)%ax = 0.0d0
  !!    enddo
  !!  enddo
  !!  ! matrix assembling
  !!  do k=1,pm%n_in
  !!    if(ssociated(pm%t2c(k)%p) then
  !!      pm%t2c(k)%p = pm%t2c(k)%p + x(k)
  !!    endif
  !!  enddo
  !! \endcode
  !!
  !! \warning The field \c t2c is a collection of pointers to the
  !! elements of the field \c. This implies that a variable of type \c
  !! t_pm_sk must must always have the \c target attribute.
  !<
  type t_pm_sk
     integer :: n_in !< number of input values
     !> map from the \c t_tri representation to the \c t_col
     !! representations of the submatrices
     type(t_rp), allocatable :: t2c(:)
     !> matrix partition as a collection of submatrices
     type(t_col), allocatable :: m(:,:)
  end type t_pm_sk

  !> Integer array (useful to pack the index arrays required by
  !! tri2col_skeleton_part).
  !<
  type t_intar
     integer, allocatable :: i(:)
  end type t_intar

  ! private members

  ! Module variables

  ! public members
  logical, protected ::               &
       mod_sparse_initialized = .false.

  character(len=*), parameter :: &
       this_mod_name = 'mod_sparse'

  interface new_col
     module procedure new_col, new_col_data, new_col_const
  end interface new_col

  interface new_tri
     module procedure new_tri, new_tri_data, new_tri_const
  end interface new_tri

  interface operator(+)
     module procedure plus_tri, plus_col, &
          plus_col_tri, plus_tri_col
  end interface operator(+)

  interface operator(-)
     module procedure minus_tri, minus_col, &
          minus_col_tri, minus_tri_col
  end interface operator(-)

  interface operator(*)
     module procedure extract_column_col, extract_row_col, &
          scal_mult_col
  end interface operator(*)

  interface sum
     module procedure sum_tri, sum_tri_dim, &
          sum_col, sum_col_dim
  end interface sum

  interface transpose
     module procedure transpose_tri, transpose_col
  end interface transpose

  interface matmul
     module procedure matmul_col, matmul_col_transp, &
          matmul_mat_mat
  end interface matmul

  interface nnz_col
     module procedure nnz_col_col
  end interface nnz_col

  interface nz_col
     module procedure nz_col_col
  end interface nz_col

  interface nz_col_i
     module procedure nz_col_col_i
  end interface nz_col_i

  interface diag
     module procedure diag_col, diag_tri, diag_col_main, diag_tri_main
  end interface diag

  interface spdiag
     module procedure spdiag_col
  end interface spdiag

  interface clear
     module procedure clear_col, clear_tri, clear_pm_sk
  end interface clear

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  subroutine mod_sparse_constructor()

    character(len=*), parameter :: &
         this_sub_name = 'constructor'

!!$    !Consistency checks ---------------------------
!!$    if( (mod_messages_initialized.eqv..false.) .or. &
!!$                                !$  ( (detailed_timing_omp.eqv..true.).and. &
!!$                                !$    (mod_omp_utils_initialized.eqv..false.) ) .or. &
!!$         (mod_kinds_initialized.eqv..false.) ) then
!!$       call error(this_sub_name,this_mod_name, &
!!$            'Not all the required modules are initialized.')
!!$    endif
!!$    if(mod_sparse_initialized.eqv..true.) then
!!$       call warning(this_sub_name,this_mod_name, &
!!$            'Module is already initialized.')
!!$    endif
!!$    !----------------------------------------------

    mod_sparse_initialized = .true.
  end subroutine mod_sparse_constructor

  !-----------------------------------------------------------------------

  subroutine mod_sparse_destructor()
    character(len=*), parameter :: &
         this_sub_name = 'destructor'

!!$    !Consistency checks ---------------------------
!!$    if(mod_sparse_initialized.eqv..false.) then
!!$       call error(this_sub_name,this_mod_name, &
!!$            'This module is not initialized.')
!!$    endif
!!$    !----------------------------------------------

    mod_sparse_initialized = .false.
  end subroutine mod_sparse_destructor

  !-----------------------------------------------------------------------

  pure function new_col(n,m,nz) result(a)
    ! construct a new matrix
    type(t_col) :: a
    integer, intent(in) :: n,m,nz

    a%n = n
    a%m = m
    a%nz = nz
    allocate( a%ap(m+1),&
         a%ai(nz), &
         a%ax(nz)  )

  end function new_col

  !-----------------------------------------------------------------------

  pure function new_col_data(n,m,ap,ai,ax) result(a)
    ! construct a new matrix and initialize it
    type(t_col) :: a
    integer, intent(in) :: n,m,ap(:),ai(:)
    sll_real64, intent(in) :: ax(:)

    if(size(ap).ne.m+1) then
       a = new_col(0,0,0)
       a%ierr = wrong_m
    elseif(size(ai).ne.size(ax)) then
       a = new_col(0,0,0)
       a%ierr = wrong_rel
    else
       a = new_col(n,m,size(ax))
       a%ap = ap
       a%ai = ai
       a%ax = ax
    endif

  end function new_col_data

  !-----------------------------------------------------------------------

  pure function new_col_const(n,m,ap,ai,x) result(a)
    ! construct a new matrix and initialize it
    type(t_col) :: a
    integer, intent(in) :: n,m,ap(:),ai(:)
    sll_real64, intent(in) :: x

    if(size(ap).ne.m+1) then
       a = new_col(0,0,0)
       a%ierr = wrong_m
    else
       a = new_col(n,m,size(ai))
       a%ap = ap
       a%ai = ai
       a%ax = x
    endif

  end function new_col_const

  !-----------------------------------------------------------------------

  pure function new_tri(n,m,nz) result(t)
    ! construct a new matrix
    type(t_tri) :: t
    integer, intent(in) :: n,m,nz

    t%n = n
    t%m = m
    t%nz = nz
    allocate( t%ti(nz), &
         t%tj(nz), &
         t%tx(nz)  )

  end function new_tri

  !-----------------------------------------------------------------------

  pure function new_tri_data(n,m,ti,tj,tx) result(t)
    ! construct a new matrix and initialize it
    type(t_tri) :: t
    integer, intent(in) :: n,m,ti(:),tj(:)
    sll_real64, intent(in) :: tx(:)

    if( (size(ti).ne.size(tj)) .or. &
         (size(ti).ne.size(tx)) ) then
       t = new_tri(0,0,0)
       t%ierr = wrong_rel
    else
       t = new_tri(n,m,size(tx))
       t%ti = ti
       t%tj = tj
       t%tx = tx
    endif

  end function new_tri_data

  !-----------------------------------------------------------------------

  pure function new_tri_const(n,m,ti,tj,x) result(t)
    ! construct a new matrix and initialize it
    type(t_tri) :: t
    integer, intent(in) :: n,m,ti(:),tj(:)
    sll_real64, intent(in) :: x

    if( (size(ti).ne.size(tj)) ) then
       t = new_tri(0,0,0)
       t%ierr = wrong_rel
    else
       t = new_tri(n,m,size(ti))
       t%ti = ti
       t%tj = tj
       t%tx = x
    endif

  end function new_tri_const

  !-----------------------------------------------------------------------

  pure function col2tri(a) result(t)
    ! col to tri conversion
    type(t_tri) :: t
    type(t_col), intent(in) :: a

    integer :: j, p, inz
    character(len=*), parameter :: &
         this_sub_name = 'col2tri'

    t = new_tri(a%n,a%m,a%nz)
    t%ierr = a%ierr
    inz = 1
    coldo: do j=1,a%m ! column loop
       rowdo: do p=a%ap(j)+1,a%ap(j+1)
          t%ti(inz) = a%ai(p)
          t%tj(inz) = j-1
          t%tx(inz) = a%ax(p)
          inz = inz+1
       enddo rowdo
    enddo coldo

  end function col2tri

  !-----------------------------------------------------------------------

  pure function tri2col(t) result(a)
    ! tri to col conversion
    type(t_col) :: a
    type(t_tri), intent(in) :: t

    logical :: inserted
    integer :: nz, inz, i, j, ip, pp
    sll_real64 :: x
    ! matrix a as row of columns
    ! this type describes one entry in a column
    type t_col_el
       type(t_col_el), pointer :: next=>null()
       integer :: i
       sll_real64 :: x
    end type t_col_el
    type(t_col_el), allocatable, target :: matrix(:)
    type(t_col_el), pointer :: el1, el2, new_el

    character(len=*), parameter :: &
         this_sub_name = 'tri2col'

    allocate( matrix(t%m) )

    ! notice that we have to consider possible repetitions
    nz = 0
    do inz=1,t%nz
       i = t%ti(inz)
       j = t%tj(inz)
       x = t%tx(inz)
       el1 => matrix(j+1) ! point to the j-th column
       inserted = .false.
       do
          if(inserted) exit
          if(.not.associated(el1%next)) then ! add new element
             allocate(el1%next)
             el1 => el1%next
             el1%i = i
             el1%x = x
             inserted = .true.
             nz = nz+1
          else
             el2 => el1%next
             if(el2%i.gt.i) then ! insert new element
                allocate(new_el)
                new_el%next => el2
                new_el%i = i
                new_el%x = x
                el1%next => new_el
                inserted = .true.
                nz = nz+1
             elseif(el2%i.eq.i) then ! duplicate entry
                el2%x = el2%x + x
                inserted = .true.
             else
                el1 => el2
             endif
          endif
       enddo
    enddo

    ! fill a and clear matrix
    a = new_col(t%n,t%m,nz)
    a%ierr = t%ierr
    a%ap(1) = 0
    ip = 0
    inz = 0
    do j=1,a%m
       ip = ip + 1
       pp = 0
       el1 => matrix(j)%next
       do
          if(.not.associated(el1)) exit
          pp = pp + 1
          inz = inz + 1
          a%ai(inz) = el1%i
          a%ax(inz) = el1%x
          el2 => el1
          el1 => el2%next
          deallocate(el2)
       enddo
       a%ap(ip+1) = a%ap(ip) + pp
    enddo

    deallocate( matrix )

  end function tri2col

  !-----------------------------------------------------------------------

  pure subroutine tri2col_skeleton(a,t2c,t)
    ! This function is similar to tri2col, but it only works on the
    ! matrix sparsity pattern, not on the actual values.
    ! In practice, the matrix a has fields a%n, a%m, a%nz, a%ap, a%ai
    ! corresponding to t, and a%ax = 0.0d0. The index vector t2c is such
    ! that
    !   t%tx(i) -> a%ax(t2c(i))
    ! or, in other words,
    !
    !   call tri2col_skeleton(a,t2c,t)
    !   do i=1,t%nz
    !     a%ax(t2c(i)) = a%ax(t2c(i)) + t%tx(i)
    !   enddo
    !
    ! is equivalent to
    !
    !   a = tri2col(t)
    !
    ! This is useful to build many matrices with the same same pattern,
    ! as in the case of time dependent problems.
    type(t_col), intent(out) :: a
    integer, allocatable, intent(out) :: t2c(:)
    type(t_tri), intent(in) :: t

    logical :: inserted
    integer :: nz, inz, i, j, ip, pp
    type t_intlist ! integer list
       type(t_intlist), pointer :: next=>null()
       integer :: i
    end type  t_intlist
    type(t_intlist), pointer :: i_intlist
    ! matrix a as row of columns
    ! this type describes one entry in a column
    type t_col_el_sk
       type(t_col_el_sk), pointer :: next=>null()
       integer :: i
       ! similar to t_col_el, except for the following field
       type(t_intlist), pointer :: it=>null()
    end type t_col_el_sk
    type(t_col_el_sk), allocatable, target :: matrix(:)
    type(t_col_el_sk), pointer :: el1, el2, new_el

    character(len=*), parameter :: &
         this_sub_name = 'tri2col_skeleton'

    allocate( matrix(t%m) )

    ! notice that we have to consider possible repetitions
    nz = 0
    do inz=1,t%nz
       i = t%ti(inz)
       j = t%tj(inz)
       el1 => matrix(j+1) ! point to the j-th column
       inserted = .false.
       do
          if(inserted) exit
          if(.not.associated(el1%next)) then ! add new element
             allocate(el1%next)
             el1 => el1%next
             el1%i = i
             allocate(el1%it); el1%it%i = inz
             inserted = .true.
             nz = nz+1
          else
             el2 => el1%next
             if(el2%i.gt.i) then ! insert new element
                allocate(new_el)
                new_el%next => el2
                new_el%i = i
                allocate(new_el%it); new_el%it%i = inz
                el1%next => new_el
                inserted = .true.
                nz = nz+1
             elseif(el2%i.eq.i) then ! duplicate entry
                i_intlist => el2%it
                ! Here is the main difference compared to tri2col: there we
                ! can sum up the contributions, but here we need to keep
                ! the various indexes separate, so the single field x used
                ! in tri2col is substituted here with a list.
                do
                   if(.not.associated(i_intlist%next)) exit
                   i_intlist => i_intlist%next
                enddo
                allocate(i_intlist%next); i_intlist%next%i = inz
                inserted = .true.
             else
                el1 => el2
             endif
          endif
       enddo
    enddo

    ! fill a and clear matrix
    a = new_col(t%n,t%m,nz)
    allocate(t2c(t%nz))
    a%ierr = t%ierr
    a%ap(1) = 0
    ip = 0
    inz = 0
    do j=1,a%m
       ip = ip + 1
       pp = 0
       el1 => matrix(j)%next
       do
          if(.not.associated(el1)) exit
          pp = pp + 1
          inz = inz + 1
          a%ai(inz) = el1%i
          i_intlist => el1%it
          do ! collect all the indexes in t2c
             if(.not.associated(i_intlist)) exit
             t2c(i_intlist%i) = inz
             el1%it => i_intlist ! recycle this to clean up the list
             i_intlist => i_intlist%next
             deallocate(el1%it)
          enddo
          el2 => el1
          el1 => el2%next
          deallocate(el2)
       enddo
       a%ap(ip+1) = a%ap(ip) + pp
    enddo

    a%ax = 0.0d0

    deallocate( matrix )

  end subroutine tri2col_skeleton

  !-----------------------------------------------------------------------

  !> Build a partitioned matrix skeleton.
  !!
  !! This subroutine is similar to \c tri2col_skeleton, but it also
  !! takes into account the fact that the matrix must be partitioned.
  !! The partition indexes can not have repetitions, but some
  !! columns/rows of the original matrix can be dropped.
  !!
  !! \note As always, when working with sparse matrices, row and column
  !! indexes in \c idi and \c idj start from 0.
  !<
  pure subroutine tri2col_skeleton_part(pm,idi,idj,t)
    type(t_intar), intent(in) :: idi(:), idj(:)
    type(t_tri), intent(in) :: t
    type(t_pm_sk), target, intent(out) :: pm

    integer :: im, jm, i, j, k, mij, nij, nzij
    integer, allocatable :: i_trans(:), j_trans(:), used(:), ti(:), tj(:), &
         t2cij(:)
    type(t_tri) :: tij
    character(len=*), parameter :: &
         this_sub_name = 'tri2col_skeleton_part'

    ! allocations
    pm%n_in = t%nz
    allocate( pm%t2c(pm%n_in) )
    allocate( pm%m(size(idi),size(idj)) )
    ! work arrays
    allocate( i_trans(t%n) , j_trans(t%m) , used(t%nz) , ti(t%nz) , tj(t%nz) )

    ! generate the skeletons
    jm_do: do jm=1,size(idj)

       ! build the j translation array
       j_trans = -1
       mij = size(idj(jm)%i) ! submatrix # columns
       do k=1,mij
          j_trans(idj(jm)%i(k)+1) = k-1
       enddo

       im_do: do im=1,size(idi)

          ! build the i translation array
          i_trans = -1
          nij = size(idi(im)%i) ! submatrix # rows
          do k=1,nij
             i_trans(idi(im)%i(k)+1) = k-1
          enddo

          ! loop over the triplets
          nzij = 0
          do k=1,t%nz

             ! translate the indexes
             i = i_trans(t%ti(k)+1)
             j = j_trans(t%tj(k)+1)

             ! check whether the triplet is used
             if( (i.ne.-1).and.(j.ne.-1) ) then
                if(associated(pm%t2c(k)%p)) then
                   ! error: each triplet can be used at most once
                   pm%m(im,jm)%ierr = wrong_rel
                   exit im_do
                endif
                nzij = nzij+1
                used(nzij) = k
                ti(nzij) = i
                tj(nzij) = j
             endif

          enddo

          tij = new_tri(nij,mij,ti(1:nzij),tj(1:nzij),0.0d0)
          call tri2col_skeleton(pm%m(im,jm),t2cij,tij)

          ! fill t2c
          do k=1,nzij
             pm%t2c(used(k))%p => pm%m(im,jm)%ax(t2cij(k))
          enddo

       enddo im_do
    enddo jm_do

    ! deallocate working arrays
    deallocate( i_trans , j_trans , used , ti , tj , t2cij )
    call clear(tij)

  end subroutine tri2col_skeleton_part

  !-----------------------------------------------------------------------

  pure function plus_tri(x,y) result(z)
    ! +
    type(t_tri) :: z
    type(t_tri), intent(in) :: x,y

    if( (x%n.ne.y%n).or.(x%m.ne.y%m) ) then
       z = new_tri(0,0,0)
       z%ierr = wrong_dim
    else
       z = new_tri(x%n,x%m,     &
            (/ x%ti , y%ti /), &
            (/ x%tj , y%tj /), &
            (/ x%tx , y%tx /))
    endif

  end function plus_tri

  !-----------------------------------------------------------------------

  pure function plus_col(x,y) result(z)
    ! +
    type(t_col) :: z
    type(t_col), intent(in) :: x,y

    z = tri2col(col2tri(x)+col2tri(y))

  end function plus_col

  !----------------------------------------------------------------------


  pure function minus_tri(x,y) result(z)
    ! -
    type(t_tri) :: z
    type(t_tri), intent(in) :: x,y

    if( (x%n.ne.y%n).or.(x%m.ne.y%m) ) then
       z = new_tri(0,0,0)
       z%ierr = wrong_dim
    else
       z = new_tri(x%n,x%m,     &
            (/ x%ti , y%ti /), &
            (/ x%tj , y%tj /), &
            (/ x%tx , -y%tx /))
    endif

  end function minus_tri

  !-----------------------------------------------------------------------

  pure function minus_col(x,y) result(z)
    ! -
    type(t_col) :: z
    type(t_col), intent(in) :: x,y

    z = tri2col(col2tri(x)-col2tri(y))

  end function minus_col

  !-----------------------------------------------------------------------

  pure function scal_mult_col(x,y) result(z)
    ! *
    type(t_col) :: z
    sll_real64, intent(in) :: x
    type(t_col), intent(in) :: y

    z = new_col(y%n,y%m,y%ap,y%ai,x*y%ax)

  end function scal_mult_col

  !-----------------------------------------------------------------------

  pure function plus_col_tri(x,y) result(z)
    ! +
    type(t_col) :: z
    type(t_col), intent(in) :: x
    type(t_tri), intent(in) :: y

    z = tri2col(col2tri(x)+y)

  end function plus_col_tri

  !-----------------------------------------------------------------------

  pure function plus_tri_col(x,y) result(z)
    ! +
    type(t_col) :: z
    type(t_tri), intent(in) :: x
    type(t_col), intent(in) :: y

    z = y + x

  end function plus_tri_col

  !-----------------------------------------------------------------------

  pure function minus_col_tri(x,y) result(z)
    ! -
    type(t_col) :: z
    type(t_col), intent(in) :: x
    type(t_tri), intent(in) :: y

    z = tri2col(col2tri(x)-y)

  end function minus_col_tri

  !-----------------------------------------------------------------------

  pure function minus_tri_col(x,y) result(z)
    ! -
    type(t_col) :: z
    type(t_tri), intent(in) :: x
    type(t_col), intent(in) :: y

    z = y - x

  end function minus_tri_col

  !-----------------------------------------------------------------------

  pure function sum_tri(t) result(s)
    ! sum of nonzero entries
    sll_real64 :: s
    type(t_tri), intent(in) :: t

    s = sum(t%tx)

  end function sum_tri

  !-----------------------------------------------------------------------

  pure function sum_tri_dim(t,dim) result(s)
    ! sum of nonzero entries along dimension dim
    type(t_tri) :: s
    integer, intent(in) :: dim
    type(t_tri), intent(in) :: t

    select case (dim)
    case(1)
       s = new_tri(1,t%m,  &
            0*t%ti,t%tj,t%tx)
    case(2)
       s = new_tri(t%n,1,  &
            t%ti,0*t%tj,t%tx)
    case default
       s = t
    end select

  end function sum_tri_dim

  !-----------------------------------------------------------------------

  pure function sum_col(a) result(s)
    ! sum of nonzero entries
    sll_real64 :: s
    type(t_col), intent(in) :: a

    s = sum(col2tri(a))

  end function sum_col

  !-----------------------------------------------------------------------

  pure function sum_col_dim(a,dim) result(s)
    ! sum of nonzero entries along dimension dim
    type(t_col) :: s
    integer, intent(in) :: dim
    type(t_col), intent(in) :: a

    s = tri2col(sum(col2tri(a),dim)) !t2)

  end function sum_col_dim

  !-----------------------------------------------------------------------

  pure function transpose_tri(t) result(tt)
    type(t_tri), intent(in) :: t
    type(t_tri) tt

    tt = new_tri(t%m,t%n,t%tj,t%ti,t%tx)
  end function transpose_tri

  !-----------------------------------------------------------------------

  pure function transpose_col(a) result(at)
    type(t_col), intent(in) :: a
    type(t_col) at

    at = tri2col(transpose(col2tri(a)))
  end function transpose_col

  !-----------------------------------------------------------------------

  pure function matmul_col(x,a) result(b)
    ! Matrix vector multiplication. The most natural way when the matrix
    ! is in column major order is  b = x^T * A
    sll_real64, intent(in) :: x(:)
    type(t_col), intent(in) :: a
    sll_real64 :: b(a%m)

    integer :: j

    if(size(x).ne.a%n) then
       b = huge(0.0d0)
    else
!!$     !$ if(detailed_timing_omp) then
!!$     !$   call omput_push_key("MatmulCol")
!!$     !$   call omput_start_timer()
!!$     !$ endif
!!$     !$omp parallel do schedule(static) &
!!$     !$omp    default(none) private(j) shared(a,x,b)
       do j=1,a%m
          b(j) = dot_product(x(nz_col_i(a,j-1)+1),nz_col(a,j-1))
       enddo
!!$     !$omp end parallel do
!!$     !$ if(detailed_timing_omp) then
!!$     !$   call omput_write_time()
!!$     !$   call omput_close_timer()
!!$     !$   call omput_pop_key()
!!$     !$ endif
    endif

  end function matmul_col

  !-----------------------------------------------------------------------

  pure function matmul_col_transp(aa,x) result(b)
    type(t_col), intent(in) :: aa
    sll_real64, intent(in) :: x(:)
    sll_real64 :: b(aa%n)

    b = matmul(x,transpose(aa))
  end function matmul_col_transp

  !-----------------------------------------------------------------------

  pure function matmul_mat_mat(a,b) result(c)
    ! The product is computed one column at a time.
    type(t_col), intent(in) :: a, b
    type(t_col) :: c

    integer :: nnzbj, nnz, nnzn, i, j
    sll_real64 :: cij
    integer, allocatable :: bj_i(:)
    sll_real64, allocatable :: ai(:), bj(:)
    integer, allocatable :: tio(:), tjo(:), ti(:), tj(:), tin(:), tjn(:)
    sll_real64, allocatable :: txo(:), tx(:), txn(:)
    type(t_col) :: at

    if(a%m.ne.b%n) then
       c = tri2col(new_tri(a%n,b%m,(/1/),(/1/),huge(0.0d0)))
    else
       allocate( ai(a%m) )
       allocate(tin(a%n),tjn(a%n),txn(a%n))
       at = transpose(a)

       nnz = 0
       allocate(ti(nnz),tj(nnz),tx(nnz))
       do j=1,b%m ! column loop
          ! get the column of b
          nnzbj = nnz_col(b,j-1)
          allocate( bj(nnzbj) , bj_i(nnzbj) )
          bj   = nz_col  (b,j-1)
          bj_i = nz_col_i(b,j-1)+1
          nnzn = 0
          do i=1,a%n ! row loop
             ai = 0.0d0 ! i-th row of a
             ai( at%ai(at%ap(i)+1:at%ap(i+1))+1 ) = &
                  at%ax(at%ap(i)+1:at%ap(i+1))
             cij = dot_product( ai(bj_i) , bj )
             if(cij.ne.0.0d0) then
                nnzn = nnzn+1
                tin(nnzn) = i-1 ! index starts from 0
                tjn(nnzn) = j-1 ! index starts from 0
                txn(nnzn) = cij
             endif
          enddo
          deallocate(bj,bj_i)
          allocate(tio(nnz),tjo(nnz),txo(nnz))
          tio = ti; tjo = tj; txo = tx
          deallocate(ti,tj,tx)
          allocate(ti(nnz+nnzn),tj(nnz+nnzn),tx(nnz+nnzn))
          ti(1:nnz) = tio; ti(nnz+1:nnz+nnzn) = tin(1:nnzn)
          tj(1:nnz) = tjo; tj(nnz+1:nnz+nnzn) = tjn(1:nnzn)
          tx(1:nnz) = txo; tx(nnz+1:nnz+nnzn) = txn(1:nnzn)
          nnz = nnz + nnzn
          deallocate(tio,tjo,txo)
       enddo
       c = tri2col(new_tri(a%n,b%m,ti,tj,tx))
       deallocate(ti ,tj ,tx )
       deallocate(tin,tjn,txn)
       deallocate( ai )
    endif

  end function matmul_mat_mat

  !-----------------------------------------------------------------------

  pure function nnz_col_col(a,j) result(n)
    ! number of nonzero entries in column j
    integer :: n
    integer, intent(in) :: j
    type(t_col), intent(in) :: a

    n = a%ap(j+2) - a%ap(j+1)

  end function nnz_col_col

  !-----------------------------------------------------------------------

  pure function nz_col_col(a,j) result(col)
    ! nonzero entries in column j
    sll_real64, allocatable :: col(:)
    integer, intent(in) :: j
    type(t_col), intent(in) :: a

    allocate( col(nnz_col_col(a,j)) )
    col = a%ax(a%ap(j+1)+1:a%ap(j+2))

  end function nz_col_col

  !-----------------------------------------------------------------------

  pure function nz_col_col_i(a,j) result(ind)
    ! indexes of nonzero entries in column j
    integer, allocatable :: ind(:)
    integer, intent(in) :: j
    type(t_col), intent(in) :: a

    allocate( ind(nnz_col_col(a,j)) )
    ind = a%ai(a%ap(j+1)+1:a%ap(j+2))

  end function nz_col_col_i

  !-----------------------------------------------------------------------

  pure function get(a,i,j) result(x)
    ! extract element i,j
    sll_real64 :: x
    integer, intent(in) :: i,j 
    type(t_col), intent(in) :: a

    integer :: pos

    x = 0.0d0
    pos = search_sorted(i,nz_col_i(a,j))
    if(pos.gt.0) x = a%ax(a%ap(j+1)+pos)

  end function get

  !-----------------------------------------------------------------------

  pure subroutine set(a,i,j,x)
    ! set element i,j
    integer, intent(in) :: i,j 
    sll_real64, intent(in) :: x
    type(t_col), intent(inout) :: a

    integer :: pos

    pos = search_sorted(i,nz_col_i(a,j))
    if(pos.gt.0) then ! the element was already present
       a%ax(a%ap(j+1)+pos) = x
    else ! new nonzero entry: we need to reallocate a
       a = a + new_tri(a%n,a%m,(/i/),(/j/),x)
    endif

  end subroutine set

  !-----------------------------------------------------------------------

  pure function search_sorted(i,v) result(pos)
    ! search i in vector v assuming v is increasing
    integer :: pos
    integer, intent(in) :: i, v(:)

    pos = count(v.le.i)
    if(pos.ne.0) then
       if(v(pos).ne.i) pos = 0
    endif

  end function search_sorted

  !-----------------------------------------------------------------------

  pure function extract_column_col(a,ind) result(a_ind)
    ! Given a vector of indexes ind, extract the corresponding columns.
    ! The index vector can contain repetitions and can have arbitrary
    ! length. The order of the operands reflect the fact that a column
    ! can be extracted by right multiplication with an "identity" matrix
    ! where part of the diagonal coefficients are set to zero.
    type(t_col), intent(in) :: a
    integer, intent(in) :: ind(:)
    type(t_col) :: a_ind

    integer :: m, nz, j
    integer, allocatable :: ap(:), ai(:)
    sll_real64, allocatable :: ax(:)

    if(maxval(ind).ge.a%m) then
       a_ind = new_col(0,0,0)
       a_ind%ierr = wrong_m
    else

       m = size(ind)

       ! count the final nonzero entries
       nz = 0
       colloop1: do j=1,m
          nz = nz + nnz_col(a,ind(j))
       enddo colloop1

       ! contruct ap ai ax
       allocate(ap(m+1),ai(nz),ax(nz))
       ap(1) = 0
       colloop2: do j=1,m
          ap(j+1) = ap(j) + nnz_col(a,ind(j))
          ai(ap(j)+1:ap(j+1)) = nz_col_i(a,ind(j))
          ax(ap(j)+1:ap(j+1)) = nz_col(a,ind(j))
       enddo colloop2

       a_ind = new_col(a%n,m,ap,ai,ax)

       deallocate(ap,ai,ax)

    endif
  end function extract_column_col

  !-----------------------------------------------------------------------

  pure function extract_row_col(ind,a) result(a_ind)
    ! analogous to extract_column_col, but rows are etracted
    integer, intent(in) :: ind(:)
    type(t_col), intent(in) :: a
    type(t_col) :: a_ind

    a_ind = transpose( transpose(a) * ind )
  end function extract_row_col

  !-----------------------------------------------------------------------

  pure function diag_col(a,diags) result(d)
    ! Extract the diagonals of a indicated in diags. Use 0 for the main
    ! diagonal, negative values for the lower diagonals and positive
    ! values for the upper diagonals.
    ! When diags(id).ne.0, some of the last elements of d are outside the
    ! bounds of the matrix and are left uninitialized. For instance, for 
    !  diags = (/ -1 , 0 , 2 /), a%n = a%m = 6
    ! we have
    !      [ a(1,0) a(0,0) a(0,2) ]
    !      [ a(2,1) a(1,1) a(1,3) ]
    !      [ a(3,2) a(2,2) a(2,4) ]
    !  d = [ a(4,3) a(3,3) a(3,5) ]
    !      [ a(5,4) a(4,4)   **   ]
    !      [   **   a(5,5)   **   ]
    type(t_col), intent(in) :: a
    integer, intent(in) :: diags(:)
    sll_real64 :: d(min(a%n,a%m),size(diags))

    integer :: id, i, i_start, i_end, shift

    do id=1,size(diags)
       ! i is the row index (counting from 1), and must be ajusted for
       ! the secondary diagonals
       i_start = max(1-diags(id),1)
       i_end   = min(a%n,a%m-diags(id))
       shift = max(-diags(id),0)
       do i=i_start,i_end
          d(i-shift,id) = get(a,i-1,i+diags(id)-1)
       enddo
    enddo

  end function diag_col

  !-----------------------------------------------------------------------

  pure function diag_col_main(a) result(d)
    ! Analogous to diag_col, but only extracts the main diagonal
    type(t_col), intent(in) :: a
    sll_real64 :: d(min(a%n,a%m))

    d = reshape(diag(a,(/0/)),(/size(d)/))

  end function diag_col_main

  !-----------------------------------------------------------------------

  pure function diag_tri(t,diags) result(d)
    type(t_tri), intent(in) :: t
    integer, intent(in) :: diags(:)
    sll_real64 :: d(min(t%n,t%m),size(diags))

    d = diag(tri2col(t),diags)

  end function diag_tri

  !-----------------------------------------------------------------------

  pure function diag_tri_main(t) result(d)
    ! Analogous to diag_tri, but only extracts the main diagonal
    type(t_tri), intent(in) :: t
    sll_real64 :: d(min(t%n,t%m))

    d = reshape(diag(t,(/0/)),(/size(d)/))

  end function diag_tri_main

  !-----------------------------------------------------------------------

  pure function spdiag_col(a,diags) result(d)
    ! Analogous to diag, but the diagonals are written into a sparse
    ! matrix. An important feature is that the indexes in diags can
    ! exceed the size of a, in which case they are ignored.
    type(t_col), intent(in) :: a
    integer, intent(in) :: diags(:)
    type(t_col) :: d

    integer :: p, j, id, i, pos, ti(a%nz), tj(a%nz)
    sll_real64 :: tx(a%nz)
    integer, allocatable :: indi(:)
    sll_real64, allocatable :: xi(:)

    p = 0
    do j=0,a%m-1
       ! here we could use the reallocation on assignment on indi, xi
       if(allocated(indi)) deallocate(indi)
       allocate(indi(nnz_col(a,j)))
       indi = nz_col_i(a,j)
       if(allocated(xi)) deallocate(xi)
       allocate(xi(nnz_col(a,j)))
       xi = nz_col(a,j)
       do id=1,size(diags)
          i = j-diags(id)
          ! check whether a(i,j) is nonzero
          pos = search_sorted(i,indi)
          if(pos.ne.0) then
             p = p+1
             ti(p) = i
             tj(p) = j
             tx(p) = xi(pos)
          endif
       enddo
    enddo

    d = tri2col(new_tri(a%m,a%n,ti(1:p),tj(1:p),tx(1:p)))

  end function spdiag_col

  !-----------------------------------------------------------------------

  pure subroutine clear_col(a)
    ! deallocate a matrix
    type(t_col), intent(out) :: a

    a%ierr = 0
    a%n = -1
    a%m = -1
    a%nz = -1
    ! allocatable fields implicitly deallocated

  end subroutine clear_col

  !-----------------------------------------------------------------------

  pure subroutine clear_tri(t)
    ! deallocate a matrix
    type(t_tri), intent(out) :: t

    t%ierr = 0
    t%n = -1
    t%m = -1
    t%nz = -1
    ! allocatable fields implicitly deallocated

  end subroutine clear_tri

  !-----------------------------------------------------------------------

  pure subroutine clear_pm_sk(pm)
    ! deallocate a matrix
    type(t_pm_sk), intent(out) :: pm

    pm%n_in = -1
    ! allocatable fields implicitly deallocated

  end subroutine clear_pm_sk

  !-----------------------------------------------------------------------

  pure function check_col(a)
    class(t_col), intent(in) :: a
    logical :: check_col
      check_col = (a%nz.gt.0).and.(a%ierr.eq.0)
  end function check_col

  !-----------------------------------------------------------------------

end module mod_sparse

