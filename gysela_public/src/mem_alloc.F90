!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!---------------------------------------------
! file : mem_alloc.f90
! date : 31/05/2005
!  array allocation
!---------------------------------------------
module mem_alloc_module
  use prec_const
  use globals, only : pglobal_id, memory_test
  implicit none
      
  ! used for memory size calculation
  integer*8, public :: max_allocate, nb_allocate
       
  !*** surdefinition for temporary allocation ***
  interface temp_allocate
     module procedure temp_allocate1d_i, &
       temp_allocate1d_d, temp_allocate1d_c, &
       temp_allocate2d_i, temp_allocate2d_d, temp_allocate2d_c, &
       temp_allocate3d_d, temp_allocate3d_c, &
       temp_allocate4d_d, temp_allocate5d_d   
  end interface
      
  !*** surdefinition for temporary deallocation ***
  interface temp_deallocate
     module procedure temp_deallocate1d_i, temp_deallocate1d_d, &
       temp_deallocate1d_c, temp_deallocate2d_i, &
       temp_deallocate2d_d, temp_deallocate2d_c, &
       temp_deallocate3d_d, temp_deallocate3d_c, &
       temp_deallocate4d_d, temp_deallocate5d_d
  end interface
      
  !*** surdefinition for global allocation ***
  interface glob_allocate
     module procedure glob_allocate1d_i, glob_allocate1d_d, &
       glob_allocate1d_c, glob_allocate2d_i, &
       glob_allocate2d_d, glob_allocate2d_c, &
       glob_allocate3d_d, glob_allocate3d_c, &
       glob_allocate4d_d, glob_allocate4d_c, glob_allocate5d_d   
  end interface
      
  !*** surdefinition for global deallocation ***
  interface glob_deallocate
     module procedure glob_deallocate1d_i, &
       glob_deallocate1d_d, glob_deallocate1d_c, &
       glob_deallocate2d_i, glob_deallocate2d_d, &
       glob_deallocate2d_c, glob_deallocate3d_d, &
       glob_deallocate3d_c, glob_deallocate4d_d, &
       glob_deallocate4d_c, glob_deallocate5d_d
  end interface
      
  character(LEN=100), parameter, private :: &
    memory_filename = "gysela_mem.out"
  integer           , save     , private :: uout_mem = 30
      
  !******************************
  contains
  !******************************
  
  !---------------------------------------- 
  ! constructor
  !----------------------------------------
  subroutine new_allocate()
    if (pglobal_id.eq.0) then
      open(uout_mem, file = memory_filename, status = 'UNKNOWN', &
        form = 'FORMATTED')
      write(uout_mem,*) ' '
      write(uout_mem,*) &
        '*********************************************************'
      write(uout_mem,*) &
        '***           GYSELA memory required in Bytes         ***'
      write(uout_mem,*) &
        '*********************************************************'
      write(uout_mem,*) ' '
      close(uout_mem)
    end if
  end subroutine new_allocate
      
  !-------------------------------------------
  ! Write memory in the file memory_filename
  !-------------------------------------------
  subroutine Write_memory(size_array,type_name,var_name)
    integer*8    , intent(in)           :: size_array
    character*(*), intent(in)           :: type_name
    character*(*), intent(in), optional :: var_name
      
    if (pglobal_id.eq.0) then
      open(uout_mem, file = memory_filename, status = 'OLD', &
        position = 'APPEND', form = 'FORMATTED')
      write(uout_mem,'(A,I20,A,A25,A35)') &
        '-> ', &
        size_array,' Bytes for ',type_name,var_name
      close(uout_mem)    
    end if
  end subroutine Write_memory
      
  !---------------------------------------- 
  ! memory allocation for a 1D array
  !----------------------------------------
  subroutine temp_allocate1d_i(array1d,begin_dim1,end_dim1,var_name)
    integer, dimension(:)    , pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * KIND(0)
    call Write_memory(size_array,'int array 1D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array1d)
      allocate(array1d(begin_dim1:end_dim1),stat=err)
      if (err.eq.0) then
        do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array1d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate1d_i
      
  !---------------------------------------- 
  subroutine temp_allocate1d_d(array1d,begin_dim1,end_dim1,var_name)
    use globals, only : plocal_id
    real(RKIND), dimension(:), pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), optional  , intent(in) :: var_name
    
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
    
    size_array = int8(end_dim1-begin_dim1+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 1D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array1d)
      allocate(array1d(begin_dim1:end_dim1),stat=err)
      if (err.eq.0) then
        do i1 = begin_dim1,end_dim1
          array1d(i1) = 0._RKIND
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array1d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate1d_d
  
  !---------------------------------------- 
  subroutine temp_allocate1d_c(array1d,begin_dim1,end_dim1,var_name)
    complex(CKIND), dimension(:), pointer    :: array1d
    integer                     , intent(in) :: begin_dim1
    integer                     , intent(in) :: end_dim1
    character*(*)    , optional , intent(in) :: var_name
      
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 1D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array1d)
      allocate(array1d(begin_dim1:end_dim1),stat=err)
      if (err.eq.0) then
        do i1 = begin_dim1,end_dim1
          array1d(i1) = 0._CKIND
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array1d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate1d_c
      
    
  !---------------------------------------- 
  ! memory allocation for a 2D array
  !----------------------------------------
  subroutine temp_allocate2d_i(array2d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,var_name)
    integer, dimension(:,:), pointer    :: array2d
    integer                , intent(in) :: begin_dim1
    integer                , intent(in) :: end_dim1
    integer                , intent(in) :: begin_dim2
    integer                , intent(in) :: end_dim2
    character*(*), optional, intent(in) :: var_name
      
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * KIND(0) 
    call Write_memory(size_array,'integer array 2D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array2d)
      allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
        stat=err)
      if (err.eq.0) then
        do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
            array2d(i1,i2) = 0
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array2d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate2d_i
      
  !----------------------------------------
  subroutine temp_allocate2d_d(array2d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,var_name)
    real(RKIND), dimension(:,:), pointer    :: array2d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1  
    integer                    , intent(in) :: begin_dim2
    integer                    , intent(in) :: end_dim2
    character(LEN=30), optional, intent(in) :: var_name
      
    integer   :: i1, i2, err
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 2D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array2d)
      allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
        stat=err)
      if (err.eq.0) then
        do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
            array2d(i1,i2) = 0._RKIND
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array2d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate2d_d
  
  !----------------------------------------  
  subroutine temp_allocate2d_c(array2d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,var_name)
    complex(CKIND), dimension(:,:), pointer    :: array2d
    integer                       , intent(in) :: begin_dim1
    integer                       , intent(in) :: end_dim1  
    integer                       , intent(in) :: begin_dim2
    integer                       , intent(in) :: end_dim2
    character*(*)       , optional, intent(in) :: var_name
      
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 2D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array2d)
      allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
        stat=err)
      if (err.eq.0) then
        do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
            array2d(i1,i2) = 0._CKIND
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array2d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate2d_c
      
  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine temp_allocate3d_d(array3d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    real(RKIND), dimension(:,:,:), pointer    :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1  
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3    
    character*(*)      , optional, intent(in) :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1)* &
      int8(end_dim3-begin_dim3+1)* KIND(0.0D0)
    call Write_memory(size_array,'double array 3D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array3d)
      allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3),stat=err)
      if (err.eq.0) then
        do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
            do i1 = begin_dim1,end_dim1
              array3d(i1,i2,i3) = 0._RKIND
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array3d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate3d_d
      
  !----------------------------------------
  subroutine temp_allocate3d_c(array3d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    complex(CKIND), dimension(:,:,:), pointer    :: array3d
    integer                         , intent(in) :: begin_dim1
    integer                         , intent(in) :: end_dim1  
    integer                         , intent(in) :: begin_dim2
    integer                         , intent(in) :: end_dim2
    integer                         , intent(in) :: begin_dim3
    integer                         , intent(in) :: end_dim3
    character*(*)         , optional, intent(in) :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1)* &
      int8(end_dim3-begin_dim3+1)* 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 3D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array3d)
      allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3),stat=err)
      if (err.eq.0) then
        do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
            do i1 = begin_dim1,end_dim1
              array3d(i1,i2,i3) = 0._CKIND
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array3d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate3d_c
      
  !---------------------------------------- 
  ! memory allocation for a 4D array
  !----------------------------------------
  subroutine temp_allocate4d_d(array4d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2, begin_dim3,end_dim3, &
    begin_dim4,end_dim4,var_name)
    real(RKIND), dimension(:,:,:,:), pointer    :: array4d
    integer                        , intent(in) :: begin_dim1
    integer                        , intent(in) :: end_dim1
    integer                        , intent(in) :: begin_dim2
    integer                        , intent(in) :: end_dim2
    integer                        , intent(in) :: begin_dim3
    integer                        , intent(in) :: end_dim3
    integer                        , intent(in) :: begin_dim4
    integer                        , intent(in) :: end_dim4
    character*(*)        , optional, intent(in) :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * &
      int8(end_dim4-begin_dim4+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 4D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array4d)
      allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
      if (err.eq.0) then
        do i4 = begin_dim4,end_dim4
          do i3 = begin_dim3,end_dim3
            do i2 = begin_dim2,end_dim2
              do i1 = begin_dim1,end_dim1
                array4d(i1,i2,i3,i4) = 0._RKIND
              end do
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array4d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate4d_d
      
  !---------------------------------------- 
  ! memory allocation for a 5D array
  !----------------------------------------
  subroutine temp_allocate5d_d(array5d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3,begin_dim4,end_dim4, &
    begin_dim5,end_dim5,var_name)
    real(RKIND), dimension(:,:,:,:,:), pointer    :: array5d
    integer                          , intent(in) :: begin_dim1
    integer                          , intent(in) :: end_dim1
    integer                          , intent(in) :: begin_dim2
    integer                          , intent(in) :: end_dim2
    integer                          , intent(in) :: begin_dim3
    integer                          , intent(in) :: end_dim3
    integer                          , intent(in) :: begin_dim4
    integer                          , intent(in) :: end_dim4
    integer                          , intent(in) :: begin_dim5
    integer                          , intent(in) :: end_dim5
    character*(*)          , optional, intent(in) :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3, i4, i5
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * &
      int8(end_dim4-begin_dim4+1) * &
      int8(end_dim5-begin_dim5+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 5D',var_name)
    if (.not.memory_test) then 
      NULLIFY(array5d)
      allocate(array5d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3,begin_dim4:end_dim4, &
        begin_dim5:end_dim5),stat=err)
      if (err.eq.0) then
        do i5 = begin_dim5,end_dim5
          do i4 = begin_dim4,end_dim4
            do i3 = begin_dim3,end_dim3
              do i2 = begin_dim2,end_dim2
                do i1 = begin_dim1,end_dim1
                  array5d(i1,i2,i3,i4,i5) = 0._RKIND
                end do
              end do
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array5d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine temp_allocate5d_d
    
      
  !---------------------------------------- 
  ! memory deallocation of array 1D
  !----------------------------------------
  subroutine temp_deallocate1d_i(array1d)
    integer, dimension(:), pointer :: array1d
      
    nb_allocate = nb_allocate - size(array1d)
    if (.not.memory_test) &
      deallocate(array1d)
  end subroutine temp_deallocate1d_i
      
  !---------------------------------------- 
  subroutine temp_deallocate1d_d(array1d)
    real(RKIND), dimension(:), pointer :: array1d
      
    nb_allocate = nb_allocate - size(array1d)
    if (.not.memory_test) &
      deallocate(array1d)
  end subroutine temp_deallocate1d_d
      
  !---------------------------------------- 
  subroutine temp_deallocate1d_c(array1d)
    complex(CKIND), dimension(:), pointer :: array1d
      
    nb_allocate = nb_allocate - size(array1d)
    if (.not.memory_test) &
      deallocate(array1d)
  end subroutine temp_deallocate1d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 2D
  !----------------------------------------
  subroutine temp_deallocate2d_i(array2d)
    integer, dimension(:,:), pointer :: array2d
      
    nb_allocate = nb_allocate - size(array2d)
    if (.not.memory_test) &
      deallocate(array2d)
  end subroutine temp_deallocate2d_i
      
  !----------------------------------------
  subroutine temp_deallocate2d_d(array2d)
    real(RKIND), dimension(:,:), pointer :: array2d
      
    nb_allocate = nb_allocate - size(array2d)
    if (.not.memory_test) &
      deallocate(array2d)
  end subroutine temp_deallocate2d_d
      
  !----------------------------------------
  subroutine temp_deallocate2d_c(array2d)
    complex(CKIND), dimension(:,:), pointer :: array2d
      
    nb_allocate = nb_allocate - size(array2d)
    if (.not.memory_test) &
      deallocate(array2d)
  end subroutine temp_deallocate2d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 3D
  !----------------------------------------
  subroutine temp_deallocate3d_d(array3d)
    real(RKIND), dimension(:,:,:), pointer :: array3d
      
    nb_allocate = nb_allocate - size(array3d)
    if (.not.memory_test) &
      deallocate(array3d)
  end subroutine temp_deallocate3d_d
      
  !----------------------------------------
  subroutine temp_deallocate3d_c(array3d)
    complex(CKIND), dimension(:,:,:), pointer :: array3d
      
    nb_allocate = nb_allocate - size(array3d)
    if (.not.memory_test) &
      deallocate(array3d)
  end subroutine temp_deallocate3d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 4D
  !----------------------------------------
  subroutine temp_deallocate4d_d(array4d)
    real(RKIND), dimension(:,:,:,:), pointer :: array4d
      
    nb_allocate = nb_allocate - size(array4d)
    if (.not.memory_test) &
      deallocate(array4d)
  end subroutine temp_deallocate4d_d
      
  !---------------------------------------- 
  ! memory deallocation of array 5D
  !----------------------------------------
  subroutine temp_deallocate5d_d(array5d)
    real(RKIND), dimension(:,:,:,:,:), pointer :: array5d
      
    nb_allocate = nb_allocate - size(array5d)
    if (.not.memory_test) &
      deallocate(array5d)
  end subroutine temp_deallocate5d_d
      
  !---------------------------------------- 
  ! memory allocation for a 1D array
  !----------------------------------------
  subroutine glob_allocate1d_i(array1d,begin_dim1,end_dim1,var_name)
    integer, dimension(:)    , pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * KIND(0)
    call Write_memory(size_array,'integer array 1D',var_name)
    if (.not.memory_test) then 
      allocate(array1d(begin_dim1:end_dim1),stat=err)
      if (err.eq.0) then
        do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array1d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate1d_i
      
  !---------------------------------------- 
  subroutine glob_allocate1d_d(array1d,begin_dim1,end_dim1,var_name)
    real(RKIND), dimension(:), pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 1D',var_name)
    if (.not.memory_test) then 
      allocate(array1d(begin_dim1:end_dim1),stat=err)
      if (err.eq.0) then
        do i1 = begin_dim1,end_dim1
          array1d(i1) = 0._RKIND
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array1d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate1d_d
      
  !---------------------------------------- 
  subroutine glob_allocate1d_c(array1d,begin_dim1,end_dim1,var_name)
    complex(CKIND), dimension(:), pointer    :: array1d
    integer                     , intent(in) :: begin_dim1
    integer                     , intent(in) :: end_dim1
    character*(*)   , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 1D',var_name)
    if (.not.memory_test) then 
      allocate(array1d(begin_dim1:end_dim1),stat=err)
      if (err.eq.0) then
        do i1 = begin_dim1,end_dim1
          array1d(i1) = 0._CKIND
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array1d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate1d_c
      
  !---------------------------------------- 
  ! memory allocation for a 2D array
  !----------------------------------------
  subroutine glob_allocate2d_i(array2d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,var_name)
    integer    , dimension(:,:), pointer    :: array2d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1
    integer                    , intent(in) :: begin_dim2
    integer                    , intent(in) :: end_dim2
    character*(*)  , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * KIND(0) 
    call Write_memory(size_array,'integer array 2D',var_name)
    if (.not.memory_test) then 
      allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
        stat=err)
      if (err.eq.0) then
        do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
            array2d(i1,i2) = 0
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array2d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate2d_i
      
  !----------------------------------------
  subroutine glob_allocate2d_d(array2d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,var_name)
    real(RKIND)  , dimension(:,:), pointer    :: array2d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1  
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    character*(*)    , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 2D',var_name)
    if (.not.memory_test) then 
      allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
        stat=err)
      if (err.eq.0) then
        do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
            array2d(i1,i2) = 0._RKIND
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array2d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate2d_d
      
  !----------------------------------------  
  subroutine glob_allocate2d_c(array2d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,var_name)
    complex(CKIND), dimension(:,:), pointer    :: array2d
    integer                       , intent(in) :: begin_dim1
    integer                       , intent(in) :: end_dim1
    integer                       , intent(in) :: begin_dim2
    integer                       , intent(in) :: end_dim2
    character*(*)     , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 2D',var_name)
    if (.not.memory_test) then 
      allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
        stat=err)
      if (err.eq.0) then
        do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
            array2d(i1,i2) = 0._CKIND
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array2d'
          end if
          print*,'-> required memory (in Bytes) = ',size(array2d)
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)    
  end subroutine glob_allocate2d_c
      
  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine glob_allocate3d_d(array3d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    real(RKIND), dimension(:,:,:), pointer    :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 3D',var_name)
    if (.not.memory_test) then 
      allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
      begin_dim3:end_dim3),stat=err)
      if (err.eq.0) then
        do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
            do i1 = begin_dim1,end_dim1
              array3d(i1,i2,i3) = 0._RKIND
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array3d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate3d_d
      
  !----------------------------------------
  subroutine glob_allocate3d_c(array3d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    complex(CKIND), dimension(:,:,:), pointer    :: array3d
    integer                         , intent(in) :: begin_dim1
    integer                         , intent(in) :: end_dim1
    integer                         , intent(in) :: begin_dim2
    integer                         , intent(in) :: end_dim2
    integer                         , intent(in) :: begin_dim3
    integer                         , intent(in) :: end_dim3
    character*(*)       , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 3D',var_name)
    if (.not.memory_test) then 
      allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3),stat=err)
      if (err.eq.0) then
        do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
            do i1 = begin_dim1,end_dim1
              array3d(i1,i2,i3) = 0._CKIND
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array3d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate3d_c
      
  !---------------------------------------- 
  ! memory allocation for a 4D array
  !----------------------------------------
  subroutine glob_allocate4d_d(array4d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3, &
    begin_dim4,end_dim4,var_name)
    real(RKIND), dimension(:,:,:,:), pointer    :: array4d
    integer                        , intent(in) :: begin_dim1
    integer                        , intent(in) :: end_dim1
    integer                        , intent(in) :: begin_dim2
    integer                        , intent(in) :: end_dim2
    integer                        , intent(in) :: begin_dim3
    integer                        , intent(in) :: end_dim3
    integer                        , intent(in) :: begin_dim4
    integer                        , intent(in) :: end_dim4
    character*(*)      , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * &
      int8(end_dim4-begin_dim4+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 4D',var_name)
    if (.not.memory_test) then 
      allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
      if (err.eq.0) then
        do i4 = begin_dim4,end_dim4
          do i3 = begin_dim3,end_dim3
            do i2 = begin_dim2,end_dim2
              do i1 = begin_dim1,end_dim1
                array4d(i1,i2,i3,i4) = 0._RKIND
              end do
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array4d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate4d_d
      
  !----------------------------------------
  subroutine glob_allocate4d_c(array4d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3, &
    begin_dim4,end_dim4,var_name)
    complex(CKIND), dimension(:,:,:,:), pointer    :: array4d
    integer                           , intent(in) :: begin_dim1
    integer                           , intent(in) :: end_dim1
    integer                           , intent(in) :: begin_dim2
    integer                           , intent(in) :: end_dim2
    integer                           , intent(in) :: begin_dim3
    integer                           , intent(in) :: end_dim3
    integer                           , intent(in) :: begin_dim4
    integer                           , intent(in) :: end_dim4
    character*(*)         , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * &
      int8(end_dim4-begin_dim4+1) * 2*KIND(0.0D0)
    call Write_memory(size_array,'complex array 4D',var_name)
    if (.not.memory_test) then 
      allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
      if (err.eq.0) then
        do i4 = begin_dim4,end_dim4
          do i3 = begin_dim3,end_dim3
            do i2 = begin_dim2,end_dim2
              do i1 = begin_dim1,end_dim1
                array4d(i1,i2,i3,i4) = 0._CKIND
              end do
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array4d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate4d_c
      
  !---------------------------------------- 
  ! memory allocation for a 5D array
  !----------------------------------------
  subroutine glob_allocate5d_d(array5d,begin_dim1,end_dim1, &
    begin_dim2,end_dim2,begin_dim3,end_dim3, &
    begin_dim4,end_dim4,begin_dim5,end_dim5,var_name)
    real(RKIND), dimension(:,:,:,:,:), pointer    :: array5d
    integer                          , intent(in) :: begin_dim1
    integer                          , intent(in) :: end_dim1
    integer                          , intent(in) :: begin_dim2
    integer                          , intent(in) :: end_dim2
    integer                          , intent(in) :: begin_dim3
    integer                          , intent(in) :: end_dim3
    integer                          , intent(in) :: begin_dim4
    integer                          , intent(in) :: end_dim4
    integer                          , intent(in) :: begin_dim5
    integer                          , intent(in) :: end_dim5
    character*(*)        , intent(in), optional   :: var_name
      
    integer   :: err
    integer   :: i1, i2, i3, i4, i5
    integer*8 :: size_array 
      
    size_array = int8(end_dim1-begin_dim1+1) * &
      int8(end_dim2-begin_dim2+1) * &
      int8(end_dim3-begin_dim3+1) * &
      int8(end_dim4-begin_dim4+1) * &
      int8(end_dim5-begin_dim5+1) * KIND(0.0D0)
    call Write_memory(size_array,'double array 5D',var_name)
    if (.not.memory_test) then 
      allocate(array5d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
        begin_dim3:end_dim3,begin_dim4:end_dim4, &
        begin_dim5:end_dim5),stat=err)
      if (err.eq.0) then
        do i5 = begin_dim5,end_dim5
          do i4 = begin_dim4,end_dim4
            do i3 = begin_dim3,end_dim3
              do i2 = begin_dim2,end_dim2
                do i1 = begin_dim1,end_dim1
                  array5d(i1,i2,i3,i4,i5) = 0._RKIND
                end do
              end do
            end do
          end do
        end do
      else
        if (pglobal_id.eq.0) then
          if (present(var_name)) then
            print*,'problem for allocated ',var_name
          else
            print*,'problem for allocated array5d'
          end if
          print*,'-> required memory (in Bytes) = ',size_array
        end if
        stop
      end if
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine glob_allocate5d_d
      
  !---------------------------------------- 
  ! memory deallocation of array 1D
  !----------------------------------------
  subroutine glob_deallocate1d_i(array1d)
    integer, dimension(:), pointer :: array1d
      
    if (associated(array1d)) then
      nb_allocate = nb_allocate - size(array1d)
      if (.not.memory_test) &
        deallocate(array1d)
    end if
  end subroutine glob_deallocate1d_i
      
  !---------------------------------------- 
  subroutine glob_deallocate1d_d(array1d)
    real(RKIND), dimension(:), pointer :: array1d
      
    if (associated(array1d)) then
       nb_allocate = nb_allocate - size(array1d)
       if (.not.memory_test) &
         deallocate(array1d)
    end if
  end subroutine glob_deallocate1d_d
      
  !---------------------------------------- 
  subroutine glob_deallocate1d_c(array1d)
    complex(CKIND), dimension(:), pointer :: array1d
      
    if (associated(array1d)) then
       nb_allocate = nb_allocate - size(array1d)
       if (.not.memory_test) &
         deallocate(array1d)
    end if
  end subroutine glob_deallocate1d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 2D
  !----------------------------------------
  subroutine glob_deallocate2d_i(array2d)
    integer, dimension(:,:), pointer :: array2d
      
    if (associated(array2d)) then
       nb_allocate = nb_allocate - size(array2d)
       if (.not.memory_test) &
         deallocate(array2d)
    end if
  end subroutine glob_deallocate2d_i
      
  !----------------------------------------
  subroutine glob_deallocate2d_d(array2d)
    real(RKIND), dimension(:,:), pointer :: array2d
      
    if (associated(array2d)) then
       nb_allocate = nb_allocate - size(array2d)
       if (.not.memory_test) &
         deallocate(array2d)
    end if
  end subroutine glob_deallocate2d_d
      
  !----------------------------------------
  subroutine glob_deallocate2d_c(array2d)
    complex(CKIND), dimension(:,:), pointer :: array2d
      
    if (associated(array2d)) then
       nb_allocate = nb_allocate - size(array2d)
       if (.not.memory_test) &
         deallocate(array2d)
    end if
  end subroutine glob_deallocate2d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 3D
  !----------------------------------------
  subroutine glob_deallocate3d_d(array3d)
    real(RKIND), dimension(:,:,:), pointer :: array3d
      
    if (associated(array3d)) then
       nb_allocate = nb_allocate - size(array3d)
       if (.not.memory_test) &
         deallocate(array3d)
    end if
  end subroutine glob_deallocate3d_d
      
  !----------------------------------------
  subroutine glob_deallocate3d_c(array3d)
    complex(CKIND), dimension(:,:,:), pointer :: array3d
      
    if (associated(array3d)) then
       nb_allocate = nb_allocate - size(array3d)
       if (.not.memory_test) &
         deallocate(array3d)
    end if
  end subroutine glob_deallocate3d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 4D
  !----------------------------------------
  subroutine glob_deallocate4d_d(array4d)
    real(RKIND), dimension(:,:,:,:), pointer :: array4d
      
    if (associated(array4d)) then
       nb_allocate = nb_allocate - size(array4d)
       if (.not.memory_test) &
         deallocate(array4d)
    end if
  end subroutine glob_deallocate4d_d
      
  !----------------------------------------
  subroutine glob_deallocate4d_c(array4d)
    complex(CKIND), dimension(:,:,:,:), pointer :: array4d
      
    if (associated(array4d)) then
       nb_allocate = nb_allocate - size(array4d)
       if (.not.memory_test) &
         deallocate(array4d)
    end if
  end subroutine glob_deallocate4d_c
      
  !---------------------------------------- 
  ! memory deallocation of array 5D
  !----------------------------------------
  subroutine glob_deallocate5d_d(array5d)
    real(RKIND), dimension(:,:,:,:,:), pointer :: array5d
      
    if (associated(array5d)) then
       nb_allocate = nb_allocate - size(array5d)
       if (.not.memory_test) &
         deallocate(array5d)
    end if
  end subroutine glob_deallocate5d_d
      
  !***********************************************
  !  function for program analysis
  !***********************************************
  subroutine print_memory_size(max_memory_size)
    use globals, only : uout_res, outputproc, file_name_res
    implicit none
    real(RKIND), intent(in) :: max_memory_size 
      
    integer :: iprint
    integer :: uout
      
    if (pglobal_id.eq.outputproc) then
       do iprint = 1,2
         if (iprint.eq.1) then
           uout = uout_res
           open(uout_res, file = file_name_res, status = 'OLD', &
             position = 'APPEND', form = 'FORMATTED')
         else
           uout = uout_mem
           open(uout_mem, file = memory_filename, status = 'OLD', &
             position = 'APPEND', form = 'FORMATTED')
         end if
         print*,' '
         write(uout,'(A)') &
           '---------------------------------------------' // &
           '--------------------------------------------------'
         write(uout,*) '--> MEMORY SIZE COMPUTATION'
         if (max_memory_size.gt.1.e12_RKIND) then
           write(uout,'(A,1f10.2,A)'), &
             'max memory size required per MPI process = ', &
             max_memory_size/1.e12_RKIND, ' TBytes'
         else
           if (max_memory_size.gt.1.e9_RKIND) then
           write(uout,'(A,1f10.2,A)'), &
             'max memory size required per MPI process = ', &
             max_memory_size/1.e9_RKIND, ' GBytes'
           else if (max_memory_size.gt.1.e6_RKIND) then
             write(uout,'(A,1f10.2,A)'), &
               'max memory size required per MPI process = ', &
               max_memory_size/1.e6_RKIND, ' MBytes'
           else 
             write(uout,'(A,1f10.2,A)'), &
               'max memory size required per MPI process = ', &
               max_memory_size/1.e3_RKIND, ' Bytes'
           end if
         end if
         write(uout,'(A)') &
           '---------------------------------------------' // &
           '--------------------------------------------------'
         if (iprint.eq.2) close(uout_mem)
       end do
    end if
  end subroutine print_memory_size
end module mem_alloc_module
