module diagnostics
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gp_plot2d(rang, nproc, f, nx, ny, hx, hy, sx, sy )

integer, intent(in) :: rang, nx, ny, nproc
real(8), intent(in), dimension(:,:) :: f
integer :: sx, sy
real(8) :: hx, hy
integer :: i, j

print*, sx, sy
!write domains
open( 80, file = "p"//char(rang+48)//".dat" )
   do i= 1, nx
      do j= 1, ny
         write(80,"(3e15.5)") (sx+i-1)/hx, (sy+j-1)/hy, f(i,j)
      end do
      write(80,*) 
   end do
close(80)
   
!write master file
if (rang == 0) then
   open( 90, file = 'p.gnu', position="append" )
   write(90,"(a)",advance='no')"splot 'p0.dat' w lines"
   do j = 1, nproc - 1
      write(90,"(a)",advance='no') ",'p"//char(j+48)//".dat' w lines" 
   end do
   write(90,*)
   close(90)
end if

end subroutine gp_plot2d

end module diagnostics
