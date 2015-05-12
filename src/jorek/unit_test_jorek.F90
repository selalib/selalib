!! test case 1 :
!!   - mesh: square, collela or square-periodic
!!   - r0=z0=1 (lenght square)
!!   - cylindrical or cartesian coordindate
!!   - dirichet or periodic boundary condition
!!   - f(x,y)= 8pi**2 * sin(2pi*r)*sin(2pi*z)
!! test case 2 : validate this case
!! test case 3 : validate this case
!! test case 4 :
!!   - mesh: square-periodic
!!   - r0=z0=1 (lenght square)
!!   - cylindrical or cartesian coordindate
!!   - periodic boundary condition
!!   - f(x,y)= 8pi**2 * cos(2pi*r)*cos(2pi*z)
program test_jorek 

  use sll_jorek
  type(sll_jorek_solver) :: jorek_solver

          
  call sll_create(jorek_solver)
  call sll_solve(jorek_solver)

  call plot_field(jorek_solver, "jorek", "field")
  call plot_jorek_field_2d_with_plotmtv(jorek_solver, "jorek")

  call sll_delete(jorek_solver)


end program test_jorek
