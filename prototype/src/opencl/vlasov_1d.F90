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

  integer    :: ierr, irec
  integer(8) :: size_in_bytes, globalsize, localsize
  integer, parameter :: iunit = 10
  integer, parameter :: source_length = 3000
  character(len = source_length) :: source
  real, allocatable  :: a(:,:), x(:,:), b(:,:), c(:,:), d(:,:)
  type(cl_mem)       :: cl_a, cl_x, cl_b, cl_c, cl_d

  integer            :: i, j, k

  integer :: i1, i2
  integer :: i_step, n_steps

  integer, parameter :: nc_eta1=128, nc_eta2=128

  real(4) :: eta1(nc_eta1+1)
  real(4) :: eta2(nc_eta2+1)
  real(4) :: eta1_new(nc_eta1+1)
  real(4) :: eta2_new(nc_eta2+1)
  real(4) :: delta_eta1
  real(4) :: delta_eta2
  real(4) :: eta1_min
  real(4) :: eta1_max
  real(4) :: eta2_min
  real(4) :: eta2_max
  real(4) :: delta_t

  real(4), dimension(nc_eta1+1)  :: advfield_v
  real(4), dimension(nc_eta2+1)  :: advfield_x

  real(4), dimension(nc_eta1+1,nc_eta2+1) :: field

  eta1_min =  -8.0; eta1_max =  8.0
  eta2_min =  -8.0; eta2_max =  8.0 

  delta_eta1 = (eta1_max-eta1_min) / nc_eta1
  delta_eta2 = (eta2_max-eta2_min) / nc_eta2

  eta1(1) = eta1_min 
  do i1 = 1, nc_eta1
     eta1(i1+1) = eta1(i1) + delta_eta1
  end do

  eta2(1) = eta2_min 
  do i2 = 1, nc_eta2
     eta2(i2+1) = eta2(i2) + delta_eta2
  end do

  forall(i1=1:nc_eta1,i2=1:nc_eta2)
     field(i1,i2) = exp(-0.5*(eta1(i1)*eta1(i1)+eta2(i2)*eta2(i2)))
  end forall
  
  advfield_x = 1.
  advfield_v = 1.

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 10
  delta_t = 3.0/n_steps

  size_in_bytes = int(nc_eta1*nc_eta2, 8)*4_8
  allocate(a(nc_eta2,nc_eta1))
  allocate(b(nc_eta2,nc_eta1))
  allocate(c(nc_eta2,nc_eta1))
  allocate(d(nc_eta2,nc_eta1))
  allocate(x(nc_eta2,nc_eta1))

  a =   1.0
  b = - 2.0
  c =   1.0
  d =   1.0
  x =   0.0

  b(:,1) = 1.; c(:,1) = 0.; d(:,1) = 0.

  a(:,nc_eta1) = 0.; b(:,nc_eta1) = 1.; d(:,nc_eta1) = 0.

  ! get the number of platforms
  call clGetPlatformIDs(num_platforms, ierr)
  write(*, '(a,i1)')   'Number of CL platforms      : ', num_platforms
  write(*, '(a)')      ''

  allocate(platforms(1:num_platforms))

  ! get an array of platforms
  call clGetPlatformIDs(platforms, num_platforms, ierr)

  ! iterate over platforms
  do iplat = 1, num_platforms

     call platform_info()

     ! get the device ID
     call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_ALL, num_devices, ierr)
     write(*, '(a,i1)') 'Number of devices           : ', num_devices
     write(*, '(a)')      ''

     allocate(devices(1:num_devices))

     call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_ALL, devices, num_devices, ierr)
    
     do idev = 1, num_devices ! Loop over device

        call device_info()

        ! create the context and the command queue
        context = clCreateContext(platforms(iplat), devices(idev), ierr)
        command_queue = clCreateCommandQueue(context, devices(idev), CL_QUEUE_PROFILING_ENABLE, ierr)

        !=====================!
        ! BUILD THE KERNEL    !
        !=====================!

        call read_kernel_source_file('./vlasov_1d.cl')

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
        kernel = clCreateKernel(prog, 'vlasov_1d', ierr)
        call clReleaseProgram(prog, ierr)

        ! allocate device memory
        cl_a = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
        cl_b = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
        cl_c = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
        cl_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  size_in_bytes, ierr)
        cl_x = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)

        do i_step = 1, n_steps

           do i2 = 1, nc_eta2

           !call compute_spline(spline1, field(:,i2), nc_eta1)

              do i1 = 1, nc_eta1+1
                 eta1_new(i1) = eta1(i1) - delta_t * advfield_x(i2)
                 if (eta1_new(i1) < eta1_min) then
                    eta1_new(i1) = eta1_new(i1) + eta1_max - eta1_min
                 else if (eta1_new(i1) > eta1_max) then
                    eta1_new(i1) = eta1_new(i1) - eta1_max + eta1_min
                 end if
              end do

           !call interpolate_array_values( spline1, eta1_new, field(:,i2) )

           end do

           call copy_data_to_device()

           call execute_kernel()

           call read_data_from_device()


           do i1 = 1, nc_eta1

              !call compute_spline(spline2, field(i1,:), nc_eta1)

              do i2 = 1, nc_eta2+1
                 eta2_new(i2) = eta2(i2) - delta_t * advfield_v(i1)
                 if (eta2_new(i2) < eta2_min) then
                    eta2_new(i2) = eta2_new(i2) + eta2_max - eta2_min
                 else if (eta2_new(i2) > eta2_max) then
                    eta2_new(i2) = eta2_new(i2) - eta2_max + eta2_min
                 end if
              end do

              !call interpolate_array_values( spline2, eta2_new, field(i1,:) )

           end do


        end do ! next time step

        call plot_field()

        !=====================!
        ! RELEASE EVERYTHING  !
        !=====================!

        call clReleaseKernel(kernel, ierr)
        call clReleaseCommandQueue(command_queue, ierr)
        call clReleaseContext(context, ierr)

     end do ! next device

  end do ! next platform

  deallocate(devices)
  deallocate(platforms)

  print *, 'Successful, exiting program.'
  write(*, '(a)')      'PASSED'


contains

  subroutine platform_info()

     ! print some info
     write(*, '(a,i1)') 'Platform number             : ', iplat
     call clGetPlatformInfo(platforms(iplat), CL_PLATFORM_VENDOR, info, ierr)
     write(*, '(2a)')   'Vendor                      : ', trim(info)
     call clGetPlatformInfo(platforms(iplat), CL_PLATFORM_NAME, info, ierr)
     write(*, '(2a)')   'Name                        : ', trim(info)
     call clGetPlatformInfo(platforms(iplat), CL_PLATFORM_VERSION, info, ierr)
     write(*, '(2a)')   'Version                     : ', trim(info)

  end subroutine platform_info

  subroutine device_info()

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

  end subroutine device_info

  subroutine read_kernel_source_file(filename)

     character(len=*) :: filename

     ! read the kernel source file
     open(unit = iunit, file = filename, access='direct', &
          status = 'old', action = 'read', iostat = ierr, recl = 1)
     if (ierr /= 0) then
        print*, 'Cannot open file '//filename
        stop
     end if

     source = ''
     irec = 1
     do
        read(unit = iunit, rec = irec, iostat = ierr) source(irec:irec)
        if (ierr /= 0) exit
        if(irec == source_length) stop 'Error: CL source file is too big'
        irec = irec + 1
     end do
     close(unit = iunit)

  end subroutine read_kernel_source_file

  subroutine execute_kernel()

     ! set the kernel arguments
     call clSetKernelArg(kernel, 0, nc_eta2, ierr)
     call clSetKernelArg(kernel, 1, nc_eta1, ierr)
     call clSetKernelArg(kernel, 2, cl_a, ierr)
     call clSetKernelArg(kernel, 3, cl_b, ierr)
     call clSetKernelArg(kernel, 4, cl_c, ierr)
     call clSetKernelArg(kernel, 5, cl_d, ierr)
     call clSetKernelArg(kernel, 6, cl_x, ierr)

     ! get the localsize for the kernel (note that the sizes are integer(8) variable)
     call clGetKernelWorkGroupInfo(kernel, devices(idev), CL_KERNEL_WORK_GROUP_SIZE, localsize, ierr)
     globalsize = int(nc_eta1, 8)
     if(mod(globalsize, localsize) /= 0) globalsize = globalsize + localsize - mod(globalsize, localsize) 

     ! execute the kernel
     call clEnqueueNDRangeKernel(command_queue, kernel, (/globalsize/), (/localsize/), ierr)
     call clFinish(command_queue, ierr)

  end subroutine execute_kernel

  subroutine copy_data_to_device()

     call clEnqueueWriteBuffer(command_queue, cl_a, cl_bool(.true.), 0_8, size_in_bytes, a(1,1), ierr)
     call clEnqueueWriteBuffer(command_queue, cl_b, cl_bool(.true.), 0_8, size_in_bytes, b(1,1), ierr)
     call clEnqueueWriteBuffer(command_queue, cl_c, cl_bool(.true.), 0_8, size_in_bytes, c(1,1), ierr)
     call clEnqueueWriteBuffer(command_queue, cl_d, cl_bool(.true.), 0_8, size_in_bytes, d(1,1), ierr)
     call clEnqueueWriteBuffer(command_queue, cl_x, cl_bool(.true.), 0_8, size_in_bytes, x(1,1), ierr)

  end subroutine copy_data_to_device


  subroutine read_data_from_device()

     ! read the resulting vector from device memory
     call clEnqueueReadBuffer(command_queue, cl_x, cl_bool(.true.), 0_8, size_in_bytes, field(1,1), ierr)

  end subroutine read_data_from_device

  subroutine plot_field()

     integer :: kk0, kk1, kk2, kk3, kk4
     character(len=4) :: fin
     character(len=1) :: aa,bb,cc,dd

     kk0 = i_step
     kk1 = kk0/1000
     aa  = char(kk1 + 48)
     kk2 = (kk0 - kk1*1000)/100
     bb  = char(kk2 + 48)
     kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
     cc  = char(kk3 + 48)
     kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
     dd  = char(kk4 + 48)
     fin = aa//bb//cc//dd

     open( 80, file = fin//".dat" )
     do i1=1,nc_eta1+1
        do i2=1,nc_eta2+1
           !write(80,"(3e10.2)") eta1(i1), eta2(i2), field(i1,i2)
           write(80,*) i1, i2, field(i1,i2)
        end do
        write(80,*) 
     end do
     close(80)
   
     open( 90, file = 'plot.gnu', position="append" )
     if ( i_step == 1 ) then
        rewind(90)
        write(90,*)"set surf"
        write(90,*)"set term x11"
     end if
     write(90,*)"set title 'Time = ",i_step*delta_t,"'"
     write(90,*)"splot '"//fin//".dat' w l"
     close(90)

  end subroutine plot_field

  subroutine compute_error()
  real(4) :: error, point
  ! compute error when Gaussian arrives at center (t=1)
  error = 0.0
  eta1 = eta1_min
  do i1 = 1, nc_eta1
     eta2 = eta2_min
     do i2 = 1, nc_eta2
        point = field(i1, i2) 
        error = max(error,abs(point-exp(-0.5*(eta1(i1)*eta1(i1)+eta2(i2)*eta2(i2)))))
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do

  print*, (nc_eta1+1)*(nc_eta2+1),' nodes 10 time steps. Error= ', error

  end subroutine compute_error
end program vlasov_1d
