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
! file : HDF5_io.f90
! date : 25/01/2006
!  array saving and reading in HDF5 format
!---------------------------------------------
module HDF5_io_module
  use prec_const
  use mem_alloc_module
      
  implicit none
      
  !******************************
  contains
  !******************************
  
  !---------------------------------------- 
  ! create HDF5 file 
  !----------------------------------------
  subroutine HDF5_create(filename,file_id,ierr)
    use HDF5
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr
      
    integer :: ierr_HDF5
      
    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)
      
    !*** Create a new file using default properties ***
    call H5Fcreate_f(trim(filename)//char(0), &
      H5F_ACC_TRUNC_F,file_id,ierr_HDF5)
    print *, pglobal_id, "create file", trim(filename)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_create
      
  !---------------------------------------- 
  ! open HDF5 file 
  !----------------------------------------
  subroutine HDF5_open(filename,file_id,ierr)
    use HDF5
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr
      
    integer :: ierr_HDF5
      
    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)
      
    !*** open the HDF5 file ***
    print *, pglobal_id, " open file ", trim(filename)
    call H5Fopen_f(trim(filename)//char(0), &
      H5F_ACC_RDONLY_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_open
      
  !---------------------------------------- 
  ! close HDF5 file 
  !----------------------------------------
  subroutine HDF5_close(file_id)
    use HDF5
    integer(HID_T), intent(in) :: file_id   ! file identifier
      
    integer :: error   ! error flag
      
    call H5Fclose_f(file_id,error)
  end subroutine HDF5_close
      
  !*************************************************
  !  HDF5 WRITING
  !*************************************************
  !---------------------------------------- 
  ! HDF5 saving for an integer 
  !----------------------------------------
  subroutine HDF5_integer_saving(file_id,int,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(in) :: int
    character(LEN=*), intent(in) :: dsetname  ! dataset name
      
    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
      
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,trim(dsetname), &
      H5T_NATIVE_INTEGER,dataspace,dataset,error)
      
    !*** Write the integer data to the dataset ***
    !***  using default transfer properties    ***
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)
      
    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_saving
      
  !---------------------------------------- 
  ! HDF5 saving for a real double
  !----------------------------------------
  subroutine HDF5_real_saving(file_id,rd,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    REAL(RKIND)     , intent(in) :: rd
    character(LEN=*), intent(in) :: dsetname  ! dataset name
      
    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
      
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,dsetname, &
      H5T_NATIVE_DOUBLE,dataspace,dataset,error)
      
    !*** Write the integer data to the dataset ***
    !***  using default transfer properties    ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)
      
    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_saving
      
  !---------------------------------------- 
  ! HDF5 saving for a 1D array of integer
  !----------------------------------------
  subroutine HDF5_array1D_saving_int(file_id,array1D,dim1,dsetname)
    use HDF5
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    integer, dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      dataspace,dataset,error)
      
    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
      
    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_int
      
  !---------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*4
  !----------------------------------------
  subroutine HDF5_array1D_saving_r4(file_id,array1D,dim1,dsetname)
    use HDF5
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    real(4), dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
      
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_REAL,&
      dataspace,dataset,error,property)
      
    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_REAL,array1D,dim,error)
      
    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_r4
      
  !-------------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*8
  !--------------------------------------------
  subroutine HDF5_array1D_saving(file_id,array1D,dim1,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)   , intent(in) :: file_id   ! file identifier
    real(RKIND), &
      dimension(:)   , intent(in) :: array1D
    integer          , intent(in) :: dim1
    character(LEN=*) , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Creates a new property for gzip dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE,&
      dataspace,dataset,error,property)
      
    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
      
    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving
      
  !---------------------------------------- 
  ! gzip HDF5 saving for a 2D array
  !----------------------------------------
  subroutine HDF5_array2D_saving(file_id,array2D,dim1,dim2,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)   , intent(in) :: file_id   ! file identifier
    real(RKIND), &
       dimension(:,:), intent(in) :: array2D
    integer          , intent(in) :: dim1, dim2
    character(LEN=*) , intent(in) :: dsetname  ! dataset name
      
    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer              :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(2)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: property   ! Property list identifier 
      
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)
      
    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)
      
    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_saving
      
  !---------------------------------------- 
  ! gzip HDF5 saving for a 3D array
  !----------------------------------------
  subroutine HDF5_array3D_saving(file_id,array3D, &
    dim1,dim2,dim3,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)    , intent(in) :: file_id   ! file identifier
    real(RKIND), &
      dimension(:,:,:), intent(in) :: array3D
    integer           , intent(in) :: dim1, dim2, dim3
    character(LEN=*)  , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
      
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    rank   = 3
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)
      
    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)
      
    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_saving
      
  !---------------------------------------- 
  ! gzip HDF5 saving for a 4D array
  !----------------------------------------
  subroutine HDF5_array4D_saving(file_id,array4d, &
    dim1,dim2,dim3,dim4,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)      , intent(in) :: file_id   ! file identifier
    real(RKIND), &
      dimension(:,:,:,:), intent(in) :: array4d
    integer             , intent(in) :: dim1, dim2, dim3, dim4
    character(LEN=*)    , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier
 
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    rank   = 4
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)
      
    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error)
      
    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array4D_saving
      
  !---------------------------------------- 
  ! gzip HDF5 saving for a 5D array
  !----------------------------------------
  subroutine HDF5_array5D_saving(file_id,array5d, &
    dim1,dim2,dim3,dim4,dim5,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)        , intent(in) :: file_id  ! file identifier
    real(RKIND), &
      dimension(:,:,:,:,:), intent(in) :: array5d
    integer               , intent(in) :: dim1, dim2
    integer               , intent(in) :: dim3, dim4, dim5
    character(LEN=*)      , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(5)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
      
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    dim(5) = dim5
    rank   = 5
    call H5Screate_simple_f(rank,dim,dataspace,error)
      
    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)
      
    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)
      
    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties         ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)
      
    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_saving
      
  !************************************************
  !  HDF5 READING
  !************************************************
  !---------------------------------------- 
  ! HDF5 reading for an integer 
  !----------------------------------------
  subroutine HDF5_integer_reading(file_id,int,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    integer         , intent(out) :: int
    character(LEN=*), intent(in)  :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type
      
    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   
      
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)
      
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_reading
      
  !---------------------------------------- 
  ! HDF5 reading for a real double
  !----------------------------------------
  subroutine HDF5_real_reading(file_id,rd,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    REAL(RKIND)     , intent(out) :: rd
    character(LEN=*), intent(in)  :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type
      
    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   
      
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)
      
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_reading
      
  !---------------------------------------- 
  ! HDF5 reading for an array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading(file_id,array1D,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    real(RKIND)     , &
        dimension(:), pointer    :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
      
    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
      
    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
      
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading
      
  !---------------------------------------- 
  ! HDF5 reading for an array 2D
  !----------------------------------------
  subroutine HDF5_array2D_reading(file_id,array2D,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)    , intent(in) :: file_id   ! file identifier
    real(RKIND)     , &
      dimension(:,:)  , pointer    :: array2D
    character(LEN=*)  , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
      
    !*** file opennning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array2D,1)
    dim(2) = size(array2D,2)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)
      
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading
      
  !---------------------------------------- 
  ! HDF5 reading for an array 3D
  !----------------------------------------
  subroutine HDF5_array3D_reading(file_id,array3D,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)    , intent(in) :: file_id   ! file identifier
    real(RKIND),  &
      dimension(:,:,:), pointer    :: array3D
    character(LEN=*)  , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
      
    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array3D,1)
    dim(2) = size(array3D,2)
    dim(3) = size(array3D,3)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)
      
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_reading
      
  !---------------------------------------- 
  ! HDF5 reading for an array 4D
  !----------------------------------------
  subroutine HDF5_array4D_reading(file_id,array4D,dsetname,ierr)
    use HDF5
    use prec_const
    integer(HID_T)      , intent(in)  :: file_id   ! file identifier
    real(RKIND)     , &
      dimension(:,:,:,:), pointer     :: array4D
    character(LEN=*)    , intent(in)  :: dsetname  ! dataset name
    integer   , optional, intent(out) :: ierr
      
    integer             :: ierr_HDF5      ! error flag
    integer             :: rank           ! dataset rank
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
      
    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,ierr_HDF5)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array4D,1)
    dim(2) = size(array4D,2)
    dim(3) = size(array4D,3)
    dim(4) = size(array4D,4)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,ierr_HDF5)
      
    !*** Closing ***
    call H5Dclose_f(dataset,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_array4D_reading
      
  !---------------------------------------- 
  ! HDF5 reading for an array 5D
  !----------------------------------------
  subroutine HDF5_array5D_reading(file_id,array5D,dsetname)
    use HDF5
    use prec_const
    integer(HID_T)        , intent(in) :: file_id   
    real(RKIND)     , &
      dimension(:,:,:,:,:), pointer    :: array5D
    character(LEN=*)      , intent(in) :: dsetname  ! dataset name
      
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
           dimension(5) :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
      
    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
      
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array5D,1)
    dim(2) = size(array5D,2)
    dim(3) = size(array5D,3)
    dim(4) = size(array5D,4)
    dim(5) = size(array5D,5)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)
      
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_reading
      
  !******************************************************
  !  HDF5 UTILS
  !******************************************************
  !----------------------------------------------------- 
  ! create the name of the HDF5 file on each processor 
  !   name+"_d<num_diag>"+"_p<proc_diag>".h5
  ! ( used for 3D storage)
  !-----------------------------------------------------
  function create_filename3D(name,diag_num,proc_num)
    use globals, only : numfmt
    character(LEN=*), intent(in)           :: name
    integer         , intent(in)           :: diag_num
    integer         , intent(in), optional :: proc_num
      
    character(LEN=50) :: create_filename3D
    character(LEN=30) :: diag_name, proc_name    
    
    write(diag_name,'('//numfmt//')') diag_num
    if (.not.present(proc_num)) then
      create_filename3D = trim(name)//&
        trim(diag_name)//".h5"
    else
      proc_name  = "_p    "
      write(proc_name(3:6),'(i4.4)') proc_num
      create_filename3D = trim(name)//trim(diag_name)&
        //trim(proc_name)//".h5"
    end if
  end function create_filename3D
      
  !----------------------------------------------------- 
  ! create the name of the HDF5 file on each processor 
  !   used for 5D storage
  !-----------------------------------------------------
  function create_filename5D(name,diag_num,proc_num)
    use globals, only : istart, Nr, Nbproc_r, &
      jstart, Ntheta, Nbproc_theta, &
      mu_id
    character(LEN=*), intent(in) :: name
    integer          , intent(in) :: diag_num
    integer          , intent(in) :: proc_num
      
    character(LEN=50) :: create_filename5D
    character(LEN=30) :: diag_name, proc_name
    character(LEN=30) :: dim1_name, dim2_name
    character(LEN=30) :: dim3_name, dim4_name, dim5_name, num
							
    write(num,'(i9)') diag_num
    diag_name ="_"//adjustl(num)
      
    write(num,'(i9)') ((istart)/((Nr+1)/Nbproc_r))
    dim1_name ="_"//adjustl(num)
      
    write(num,'(i9)') (jstart)/(Ntheta/Nbproc_theta)
    dim2_name ="-"//adjustl(num)
      
    dim3_name  = trim("-0")
      
    dim4_name  = trim("-0")
      
    write(num,'(i9)') mu_id
    dim5_name ="-"//adjustl(num)
      
    create_filename5D = trim(name)//trim(diag_name)//&
      trim(dim1_name)//trim(dim2_name)//trim(dim3_name)//&
      trim(dim4_name)//trim(dim5_name)//".h5"
  end function create_filename5D
      
  !------------------------------------------------ 
  ! create the name of the variable corresponding 
  !  to the tree in the HDF5 master file
  !------------------------------------------------
  function create_variable_name(tree,var_name)
    character(len=*), intent(in) :: tree, var_name
      
    character(LEN=50) :: create_variable_name
    
    create_variable_name = trim(tree)//"/"
    create_variable_name = trim(create_variable_name)//&
      trim(var_name)
  end function create_variable_name
      
  !------------------------------------------------ 
  !  Write a test file for 1D array
  !------------------------------------------------
  subroutine Write_HDF5_test1D(idiag_num,var1D_name,array1D)
    use HDF5
    integer                       , intent(in) :: idiag_num
    character(LEN=*)              , intent(in) :: var1D_name
    real(RKIND)     , dimension(:), pointer    :: array1D
      
    integer(HID_T)    :: file_id   
    integer           :: lbound1, ubound1, dim1
    character(LEN=50) :: file1D_name, var1Dname_tmp
    integer           :: nbchar
    
    var1Dname_tmp = trim(var1D_name)
    nbchar        = len_trim(var1D_name)
    if (idiag_num.ne.-1) then
      file1D_name = trim(var1D_name)//"_d    _p    .h5"
      write(file1D_name(nbchar+3:nbchar+6),'(i4.4)') idiag_num
      write(file1D_name(nbchar+9:nbchar+12),'(i4.4)') pglobal_id
    else
      file1D_name = trim(var1D_name)//"_dinit_p    .h5"
      write(file1D_name(nbchar+9:nbchar+12),'(i4.4)') pglobal_id
    end if
      
    lbound1 = lbound(array1D,1)
    ubound1 = ubound(array1D,1)
    dim1    = abs(ubound1-lbound1)+1 
    call HDF5_create(trim(file1D_name),file_id)
    call HDF5_integer_saving(file_id,lbound1,'lbound1'//char(0))
    call HDF5_integer_saving(file_id,ubound1,'ubound1'//char(0))
    call HDF5_array1D_saving(file_id,array1D,dim1,var1Dname_tmp)
    call HDF5_close(file_id)
  end subroutine Write_HDF5_test1D
      
  !------------------------------------------------ 
  !  Write a test file for 2D array
  !------------------------------------------------
  subroutine Write_HDF5_test2D(idiag_num,var2D_name,array2D)
    use HDF5
    integer                         , intent(in) :: idiag_num
    character(LEN=*)                , intent(in) :: var2D_name
    real(RKIND)     , dimension(:,:), pointer    :: array2D
      
    integer(HID_T)    :: file_id   
    integer           :: lbound1, ubound1, dim1
    integer           :: lbound2, ubound2, dim2
    character(LEN=50) :: file2D_name, var2Dname_tmp
    integer           :: nbchar
    
    var2Dname_tmp = trim(var2D_name)
    nbchar        = len_trim(var2D_name)
    if (idiag_num.ne.-1) then
      file2D_name = trim(var2D_name)//"_d    _p    .h5"
      write(file2D_name(nbchar+3:nbchar+6),'(i4.4)') idiag_num
      write(file2D_name(nbchar+9:nbchar+12),'(i4.4)') pglobal_id
    else
      file2D_name = trim(var2D_name)//"_dinit_p    .h5"
      write(file2D_name(nbchar+9:nbchar+12),'(i4.4)') pglobal_id
    end if
    lbound1 = lbound(array2D,1)
    ubound1 = ubound(array2D,1)
    dim1    = abs(ubound1-lbound1)+1 
    lbound2 = lbound(array2D,2)
    ubound2 = ubound(array2D,2)
    dim2    = abs(ubound2-lbound2)+1 
    call HDF5_create(trim(file2D_name),file_id)
    call HDF5_integer_saving(file_id,lbound1,'lbound1'//char(0))
    call HDF5_integer_saving(file_id,ubound1,'ubound1'//char(0))
    call HDF5_integer_saving(file_id,lbound2,'lbound2'//char(0))
    call HDF5_integer_saving(file_id,ubound2,'ubound2'//char(0))
    call HDF5_array2D_saving(file_id,array2D, &
      dim1,dim2,var2Dname_tmp)
    call HDF5_close(file_id)
  end subroutine Write_HDF5_test2D
      
  !------------------------------------------------ 
  !  Write a test file for 3D array
  !------------------------------------------------
  subroutine Write_HDF5_test3D(idiag_num,var3D_name,array3D)
    use HDF5
    integer                           , intent(in) :: idiag_num
    character(LEN=*)                  , intent(in) :: var3D_name
    real(RKIND)     , dimension(:,:,:), pointer    :: array3D
      
    integer(HID_T)    :: file_id   
    integer           :: lbound1, ubound1, dim1
    integer           :: lbound2, ubound2, dim2
    integer           :: lbound3, ubound3, dim3
    character(LEN=50) :: file3D_name, var3Dname_tmp
    integer           :: nbchar
    
    var3Dname_tmp = trim(var3D_name)
    nbchar        = len_trim(var3D_name)
    if (idiag_num.ne.-1) then
      file3D_name = trim(var3D_name)//"_d    _p    .h5"
      write(file3D_name(nbchar+3:nbchar+6),'(i4.4)') idiag_num
      write(file3D_name(nbchar+9:nbchar+12),'(i4.4)') pglobal_id
    else
      file3D_name = trim(var3D_name)//"_dinit_p    .h5"
      write(file3D_name(nbchar+9:nbchar+12),'(i4.4)') pglobal_id
    end if
    lbound1 = lbound(array3D,1)
    ubound1 = ubound(array3D,1)
    dim1    = abs(ubound1-lbound1)+1 
    lbound2 = lbound(array3D,2)
    ubound2 = ubound(array3D,2)
    dim2    = abs(ubound2-lbound2)+1 
    lbound3 = lbound(array3D,3)
    ubound3 = ubound(array3D,3)
    dim3    = abs(ubound3-lbound3)+1 
    call HDF5_create(trim(file3D_name),file_id)
    call HDF5_integer_saving(file_id,lbound1,'lbound1'//char(0))
    call HDF5_integer_saving(file_id,ubound1,'ubound1'//char(0))
    call HDF5_integer_saving(file_id,lbound2,'lbound2'//char(0))
    call HDF5_integer_saving(file_id,ubound2,'ubound2'//char(0))
    call HDF5_integer_saving(file_id,lbound3,'lbound3'//char(0))
    call HDF5_integer_saving(file_id,ubound3,'ubound3'//char(0))
    call HDF5_array3D_saving(file_id,array3D, &
      dim1,dim2,dim3,var3Dname_tmp)
    call HDF5_close(file_id)
  end subroutine Write_HDF5_test3D
end module HDF5_io_module
