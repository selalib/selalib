program vlasov_1d
  use cl

  implicit none

  type(cl_platform_id), allocatable :: platforms(:)
  type(cl_device_id),   allocatable :: devices(:)
  integer                           :: num_platforms
  integer                           :: num_devices
  integer                           :: iplat
  integer                           :: idev
  integer(8)                        :: val
  character(len=200)                :: info

  type(cl_context)       :: context
  type(cl_command_queue) :: command_queue
  type(cl_program)       :: prog
  type(cl_kernel)        :: kernel

  integer    :: ierr, irec, system_size, num_systems
  integer(8) :: size_in_bytes, globalsize, localsize
  integer, parameter :: iunit = 10
  integer, parameter :: source_length = 3000
  character(len = source_length) :: source
  real, allocatable  :: a(:,:), x(:,:), b(:,:), c(:,:), d(:,:)
  type(cl_mem)       :: cl_a, cl_x, cl_b, cl_c, cl_d

  integer :: i,j

  character (len=12), parameter :: kernel_file = "./tridiag.cl"

  !=====================
  ! INITIALIZATION
  !=====================

  num_systems = 128
  system_size = 128
  size_in_bytes = int(system_size*num_systems, 8)*4_8
  allocate(a(num_systems,system_size))
  allocate(b(num_systems,system_size))
  allocate(c(num_systems,system_size))
  allocate(d(num_systems,system_size))
  allocate(x(num_systems,system_size))

  ! get the number of platforms
  call clGetPlatformIDs(num_platforms, ierr)
  write(*, '(a,i1)')   'Number of CL platforms      : ', num_platforms
  write(*, '(a)')      ''

  allocate(platforms(1:num_platforms))

  ! get an array of platforms
  call clGetPlatformIDs(platforms, num_platforms, ierr)

  ! iterate over platforms
  do iplat = 1, num_platforms

    ! print some info
    write(*, '(a,i1)') 'Platform number             : ', iplat
    call clGetPlatformInfo(platforms(iplat), CL_PLATFORM_VENDOR, info, ierr)
    write(*, '(2a)')   'Vendor                      : ', trim(info)
    call clGetPlatformInfo(platforms(iplat), CL_PLATFORM_NAME, info, ierr)
    write(*, '(2a)')   'Name                        : ', trim(info)
    call clGetPlatformInfo(platforms(iplat), CL_PLATFORM_VERSION, info, ierr)
    write(*, '(2a)')   'Version                     : ', trim(info)
    ! get the device ID
    call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_ALL, num_devices, ierr)
    write(*, '(a,i1)') 'Number of devices           : ', num_devices
    write(*, '(a)')      ''

    allocate(devices(1:num_devices))

    call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_ALL, devices, num_devices, ierr)
    
    do idev = num_devices,1,-1

      write(*, '(a)')      ''
      write(*, '(a,i1)') '    Device number           : ', idev

      call clGetDeviceInfo(devices(idev), CL_DEVICE_TYPE, val, ierr)
      select case(val)
      case(CL_DEVICE_TYPE_CPU)
        info = 'CPU'
      case(CL_DEVICE_TYPE_GPU)
        info = 'GPU'
      case(CL_DEVICE_TYPE_ACCELERATOR)
        info = 'Accelerator'
      end select
   
      write(*, '(2a)')   '    Device type             : ', trim(info)
      call clGetDeviceInfo(devices(idev), CL_DEVICE_VENDOR, info, ierr)
      write(*, '(2a)')   '    Device vendor           : ', trim(info)
      call clGetDeviceInfo(devices(idev), CL_DEVICE_NAME, info, ierr)
      write(*, '(2a)')   '    Device name             : ', trim(info)
      call clGetDeviceInfo(devices(idev), CL_DEVICE_GLOBAL_MEM_SIZE, val, ierr)
      write(*, '(a,i4)') '    Device memory           : ', val/1024**2
      write(*, '(a)')      ''

      ! create the context and the command queue
      context = clCreateContext(platforms(iplat), devices(idev), ierr)
      command_queue = clCreateCommandQueue(context, devices(idev), CL_QUEUE_PROFILING_ENABLE, ierr)

      !=====================!
      ! BUILD THE KERNEL    !
      !=====================!

      ! read the source file
      open(unit = iunit, file = kernel_file, access='direct', &
           status = 'old', action = 'read', iostat = ierr, recl = 1)
      if (ierr /= 0) stop 'Cannot open file '//kernel_file

      source = ''
      irec = 1
      do
        read(unit = iunit, rec = irec, iostat = ierr) source(irec:irec)
        if (ierr /= 0) exit
        if(irec == source_length) stop 'Error: CL source file is too big'
        irec = irec + 1
      end do
      close(unit = iunit)

      ! create the program
      prog = clCreateProgramWithSource(context, source, ierr)
      if(ierr /= CL_SUCCESS) stop 'Error: cannot create program from source.'

      ! build
      call clBuildProgram(prog, '-cl-mad-enable', ierr)

      ! get the compilation log
      call clGetProgramBuildInfo(prog, devices(idev), CL_PROGRAM_BUILD_LOG, source, irec)
      if(len(trim(source)) > 0) print*, trim(source)

      if(ierr /= CL_SUCCESS) stop 'Error: program build failed.'

      ! finally get the kernel and release the program
      kernel = clCreateKernel(prog, 'tridiag', ierr)
      call clReleaseProgram(prog, ierr)

      !=====================!
      ! RUN THE KERNEL      !
      !=====================!
  
      a = 1.0
      b = -2.0
      c = 1.0
      d = 1.0
      x = 0.0

      b(:,1) = 1.; c(:,1) = 0.; d(:,1) = 0.

      a(:,system_size) = 0.; b(:,system_size) = 1.; d(:,system_size) = 0.

      ! allocate device memory
      cl_a = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
      cl_b = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
      cl_c = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
      cl_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
      cl_x = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)

      ! copy data to device memory
      call clEnqueueWriteBuffer(command_queue, cl_a, cl_bool(.true.), 0_8, size_in_bytes, a(1,1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_b, cl_bool(.true.), 0_8, size_in_bytes, b(1,1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_c, cl_bool(.true.), 0_8, size_in_bytes, c(1,1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_d, cl_bool(.true.), 0_8, size_in_bytes, d(1,1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_x, cl_bool(.true.), 0_8, size_in_bytes, x(1,1), ierr)

      ! set the kernel arguments
      call clSetKernelArg(kernel, 0, num_systems, ierr)
      call clSetKernelArg(kernel, 1, system_size, ierr)
      call clSetKernelArg(kernel, 2, cl_a, ierr)
      call clSetKernelArg(kernel, 3, cl_b, ierr)
      call clSetKernelArg(kernel, 4, cl_c, ierr)
      call clSetKernelArg(kernel, 5, cl_d, ierr)
      call clSetKernelArg(kernel, 6, cl_x, ierr)

      ! get the localsize for the kernel (note that the sizes are integer(8) variable)
      call clGetKernelWorkGroupInfo(kernel, devices(idev), CL_KERNEL_WORK_GROUP_SIZE, localsize, ierr)
      globalsize = int(system_size, 8)
      if(mod(globalsize, localsize) /= 0) globalsize = globalsize + localsize - mod(globalsize, localsize) 

      ! execute the kernel
      call clEnqueueNDRangeKernel(command_queue, kernel, (/globalsize/), (/localsize/), ierr)
      call clFinish(command_queue, ierr)

      ! read the resulting vector from device memory
      call clEnqueueReadBuffer(command_queue, cl_x, cl_bool(.true.), 0_8, size_in_bytes, x(1,1), ierr)
      
      open(11,file='plot'//char(idev+48)//'.dat')
      do j = 1, system_size
         do i = 1, num_systems
            write(11,*) i, j, x(i,j)
         end do
         write(11,*)
      end do
      close(11)

      !=====================!
      ! RELEASE EVERYTHING  !
      !=====================!

      call clReleaseKernel(kernel, ierr)
      call clReleaseCommandQueue(command_queue, ierr)
      call clReleaseContext(context, ierr)

    end do

    deallocate(devices)

  end do

  deallocate(platforms)
  write(*, '(a)')      'PASSED'

  a = 1.0
  b = -2.0
  c = 1.0
  d = 1.0

  b(:,1) = 1.
  c(:,1) = 0.
  d(:,1) = 0.

  a(:,system_size) = 0.
  b(:,system_size) = 1.
  d(:,system_size) = 0.

  do i = 1, num_systems

     call trian3( a(i,:), b(i,:), c(i,:), 1, system_size ) 
     call resol3( a(i,:), b(i,:), c(i,:), d(i,:), 1, system_size )

     do j = 1, system_size
        write(13,*)i,j, d(i,j)
     end do
     write(13,*)

  end do

end program vlasov_1d


subroutine trian3( tabi, tabd, tabs, m, n ) 
!         *      the first swept of the          *
!         *         tridiagonal solver           *

integer, intent(in) :: m, n
real(4), intent(in) :: tabi( n )
real(4), intent(inout) :: tabd( n ),tabs ( n )
integer :: m1, n1, i

m1 = m + 1
n1 = n - 1

tabs( m ) = tabs( m ) / tabd( m )
do i = m1, n1
   tabd( i ) = tabd( i ) - tabi( i ) * tabs( i-1 )
   tabs( i ) = tabs( i ) / tabd( i )
end do

tabd( n ) = tabd( n ) - tabi( n ) * tabs( n1 )

return
end subroutine trian3

subroutine resol3( tabi, tabd, tabs, sm, m, n )
!         *      the 2nd   swept of the          *
!         *         tridiagonal solver           *

integer, intent(in) :: m, n
real(4), intent(in) ::  tabd(n), tabi(n), tabs(n)
real(4), intent(inout) :: sm(n)
integer :: m1, n1, i, k
n1 = n - 1
m1 = m + 1
sm( m ) = sm( m ) / tabd( m )
do i = m1, n
   sm( i ) = ( sm( i ) - tabi( i ) * sm( i-1 ) ) / tabd( i )
end do
do k = m, n1
   i = n1 - k + m
   sm( i ) = sm( i ) - tabs( i ) * sm( i+1 )
end do
return
end subroutine resol3

