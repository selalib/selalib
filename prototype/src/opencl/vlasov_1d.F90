!! Copyright (C) 2011 X. Andrade
!!
!! FortranCL is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! FortranCL is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

program vlasov_1d
  use cl

  implicit none

  type(cl_platform_id), allocatable :: platforms(:)
  type(cl_device_id),   allocatable :: devices(:)
  integer :: num_platforms, num_devices, iplat, idev
  integer(8) :: val
  character(len=200) :: info

  type(cl_context)       :: context
  type(cl_command_queue) :: command_queue
  type(cl_program)       :: prog
  type(cl_kernel)        :: kernel

  integer    :: ierr, irec, size_system, num_systems
  integer(8) :: size_in_bytes, globalsize, localsize
  integer, parameter :: iunit = 10
  integer, parameter :: source_length = 1000
  character(len = source_length) :: source
  real, allocatable  :: a(:), x(:), b(:), c(:), d(:)
  type(cl_mem)       :: cl_a, cl_x, cl_b, cl_c, cl_d

  !=====================
  ! INITIALIZATION
  !=====================

  num_systems = 1
  size_system = 10
  size_in_bytes = int(size_system, 8)*4_8
  allocate(a(1:size_system))
  allocate(b(1:size_system))
  allocate(c(1:size_system))
  allocate(d(1:size_system))
  allocate(x(1:size_system))

  ! get the number of platforms
  call clGetPlatformIDs(num_platforms, ierr)

  allocate(platforms(1:num_platforms))

  write(*, '(a,i1)')   'Number of CL platforms      : ', num_platforms
  write(*, '(a)')      ''

  
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

    ! get the device ID
    call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_ALL, devices, num_devices, ierr)
    
    do idev = 1, num_devices
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

      !=====================
      ! BUILD THE KERNEL
      !=====================

      ! read the source file
      open(unit = iunit, file = './vlasov_1d.cl', access='direct', &
           status = 'old', action = 'read', iostat = ierr, recl = 1)
      if (ierr /= 0) stop 'Cannot open file vlasov_1d.cl'

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
      kernel = clCreateKernel(prog, 'vlasov_1d', ierr)
      call clReleaseProgram(prog, ierr)

      !=====================
      ! RUN THE KERNEL
      !=====================
  
      a = 1.0
      b = 2.0
      c = 3.0
      d = 4.0
      x = 0.0

      ! allocate device memory
      cl_a = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
      cl_b = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
      cl_c = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
      cl_d = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
      cl_x = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)

      ! copy data to device memory
      call clEnqueueWriteBuffer(command_queue, cl_a, cl_bool(.true.), 0_8, size_in_bytes, a(1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_b, cl_bool(.true.), 0_8, size_in_bytes, b(1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_c, cl_bool(.true.), 0_8, size_in_bytes, c(1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_d, cl_bool(.true.), 0_8, size_in_bytes, d(1), ierr)
      call clEnqueueWriteBuffer(command_queue, cl_x, cl_bool(.true.), 0_8, size_in_bytes, x(1), ierr)

      ! set the kernel arguments
      call clSetKernelArg(kernel, 0, size_system, ierr)
      call clSetKernelArg(kernel, 1, num_systems, ierr)
      call clSetKernelArg(kernel, 2, cl_a, ierr)
      call clSetKernelArg(kernel, 3, cl_b, ierr)
      call clSetKernelArg(kernel, 4, cl_c, ierr)
      call clSetKernelArg(kernel, 5, cl_d, ierr)
      call clSetKernelArg(kernel, 6, cl_x, ierr)

      ! get the localsize for the kernel (note that the sizes are integer(8) variable)
      call clGetKernelWorkGroupInfo(kernel, devices(idev), CL_KERNEL_WORK_GROUP_SIZE, localsize, ierr)
      globalsize = int(size_system, 8)
      if(mod(globalsize, localsize) /= 0) globalsize = globalsize + localsize - mod(globalsize, localsize) 

      ! execute the kernel
      call clEnqueueNDRangeKernel(command_queue, kernel, (/globalsize/), (/localsize/), ierr)
      call clFinish(command_queue, ierr)

      ! read the resulting vector from device memory
      call clEnqueueReadBuffer(command_queue, cl_x, cl_bool(.true.), 0_8, size_in_bytes, x(1), ierr)

      write(*,*) x

      !=====================
      ! RELEASE EVERYTHING
      !=====================

      call clReleaseKernel(kernel, ierr)
      call clReleaseCommandQueue(command_queue, ierr)
      call clReleaseContext(context, ierr)

    end do

    deallocate(devices)

  end do


  deallocate(platforms)
  write(*, '(a)')      'PASSED'

end program vlasov_1d
