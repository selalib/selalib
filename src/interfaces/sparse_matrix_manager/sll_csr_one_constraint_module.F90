! We add one constraint for a matrix in CSR format
! form a matrix A
! we get a matrix
!   A   b
!   b^T 0
! where b in the vector of constraint (constraint_vec)
! so that all the component of b are nonzero
! b^T stands for the transposition of b
! we suppose that A is square

!contact: Adnane Hamiaz (hamiaz@math.unistra.fr)
!         Michel Mehrenberger (mehrenbe@math.unistra.fr)
module csr_one_constraint_module

contains 
  subroutine csr_add_one_constraint( &
    ia_in, &
    ja_in, &
    a_in, &
    num_rows_in, &
    num_nz_in, &
    constraint_vec, &
    ia_out, &
    ja_out, &
    a_out)
    integer, dimension(:), intent(in) :: ia_in  
    integer, dimension(:), intent(in) :: ja_in  
    real(8), dimension(:), intent(in) :: a_in
    integer, intent(in) :: num_rows_in
    integer, intent(in) :: num_nz_in
    real(8), dimension(:), intent(in) :: constraint_vec
    integer, dimension(:), intent(out) :: ia_out
    integer, dimension(:), intent(out) :: ja_out
    real(8), dimension(:), intent(out) :: a_out
    integer :: num_rows_out
    integer :: num_nz_out
    integer :: i
    integer :: s
    
    
    num_rows_out = num_rows_in+1
    num_nz_out = num_nz_in+2*num_rows_in
    
    if(size(ia_in)<num_rows_in+1) then
      print *, '#problem of size of ia_in', size(ia_in),num_rows_in+1
      stop
    endif
    if(size(ja_in)<num_nz_in) then
      print *, '#problem of size of ja_in', size(ja_in),num_nz_in
      stop
    endif
    if(size(a_in)<num_nz_in) then
      print *, '#problem of size of a_in', size(a_in),num_nz_in
      stop
    endif
    if(size(ia_out)<num_rows_out+1) then
      print *, '#problem of size of ia_out', size(ia_out),num_rows_out+1
      stop
    endif
    if(size(ja_out)<num_nz_out) then
      print *, '#problem of size of ja_out', size(ja_out),num_nz_out
      stop
    endif
    if(size(a_out)<num_nz_out) then
      print *, '#problem of size of a_out', size(a_out),num_nz_out
      stop
    endif
    if(ia_in(num_rows_in+1).ne.num_nz_in+1)then
      print *,'#bad value of ia_in(num_rows_in+1)', ia_in(num_rows_in+1),num_nz_in+1
      stop
    endif
    
    s = 1
    do i=1,num_rows_in
      ia_out(i) = s
      do k = ia_in(i), ia_in(i+1)-1
        a_out(s) = a_in(k)
        ja_out(s) = ja_in(k)
        s = s+1
      enddo
      a_out(s) = constraint_vec(i)
      ja_out(s) = num_rows_out
      s = s+1
    enddo
    ia_out(num_rows_in+1) = s
    do i=1,num_rows_in
      a_out(s) = constraint_vec(i)
      ja_out(s) = i
      s = s+1      
    enddo
    ia_out(num_rows_in+2) = s
     
    if(ia_out(num_rows_out+1).ne.num_nz_out+1)then
      print *,'#bad value of ia_out(num_rows_out+1)',ia_out(num_rows_out+1),num_nz_out+1
      stop
    endif
    
  end subroutine csr_add_one_constraint
  
  
  
  
end module csr_one_constraint_module
