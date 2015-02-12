program interpolation

  use spline_periodic_1d

  implicit none
  integer :: i1, i2
  integer :: i_step, n_steps

  integer, parameter :: nc_eta1 = 100
  integer, parameter :: nc_eta2 = 100

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

  real(4) :: error, val

  type(spline_t) :: spline1, spline2


  integer :: kk0, kk1, kk2, kk3, kk4
  character(len=4) :: fin
  character(len=1) :: aa,bb,cc,dd


  eta1_min =  -8.0; eta1_max =  8.0
  eta2_min =  -8.0; eta2_max =  8.0 

  delta_eta1 = (eta1_max-eta1_min) / nc_eta1
  delta_eta2 = (eta2_max-eta2_min) / nc_eta2

  eta1(1) = eta1_min 
  do i1 = 1, nc_eta1
     eta1(i1+1) = eta1(i1) + delta_eta1
  end do

  call new_spline( spline1, eta1, nc_eta1)

  eta2(1) = eta2_min 
  do i2 = 1, nc_eta2
     eta2(i2+1) = eta2(i2) + delta_eta2
  end do

  call new_spline( spline2, eta2, nc_eta2)

  forall(i1=1:nc_eta1,i2=1:nc_eta2)
     field(i1,i2) = exp(-0.5*(eta1(i1)*eta1(i1)+eta2(i2)*eta2(i2)))
  end forall
  
  advfield_x = 1.
  advfield_v = 1.

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 100
  delta_t = 3.0/n_steps

  do i_step = 1, n_steps

     do i2 = 1, nc_eta2

        call compute_spline(spline1, field(:,i2), nc_eta1)

        do i1 = 1, nc_eta1+1
           eta1_new(i1) = eta1(i1) - delta_t * advfield_x(i2)
           if (eta1_new(i1) < eta1_min) then
              eta1_new(i1) = eta1_new(i1) + eta1_max - eta1_min
           else if (eta1_new(i1) > eta1_max) then
              eta1_new(i1) = eta1_new(i1) - eta1_max + eta1_min
           end if
        end do

        call interpolate_array_values( spline1, eta1_new, field(:,i2) )

     end do

     do i1 = 1, nc_eta1

        call compute_spline(spline2, field(i1,:), nc_eta1)

        do i2 = 1, nc_eta2+1
           eta2_new(i2) = eta2(i2) - delta_t * advfield_v(i1)
           if (eta2_new(i2) < eta2_min) then
              eta2_new(i2) = eta2_new(i2) + eta2_max - eta2_min
           else if (eta2_new(i2) > eta2_max) then
              eta2_new(i2) = eta2_new(i2) - eta2_max + eta2_min
           end if
        end do

        call interpolate_array_values( spline2, eta2_new, field(i1,:) )

     end do

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
           write(80,"(3e10.2)") eta1(i1), eta2(i2), field(i1,i2)
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

  end do
    
  ! compute error when Gaussian arrives at center (t=1)

  error = 0.0
  eta1 = eta1_min
  do i1 = 1, nc_eta1
     eta2 = eta2_min
     do i2 = 1, nc_eta2
        val = field(i1, i2) 
        error = max(error,abs(val-exp(-0.5*(eta1(i1)*eta1(i1)+eta2(i2)*eta2(i2)))))
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do

  print*, ' 100 nodes, 10 time steps. Error= ', error

  print *, 'Successful, exiting program.'

end program 
