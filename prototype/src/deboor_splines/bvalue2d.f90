! Just for the record: what is ar_x, ar_y, ai_nx, etc., etc.? What is this
! function supposed to do? This is a deeply frustrating file.
subroutine bvalue2d(&
     ar_x,&
     ar_y,&
     ai_nx,&
     ai_kx,&
     ai_ny,&
     ai_ky,&
     apr_Bcoef,&
     apr_tx,&
     sz_tx,&
     apr_ty,&
     sz_ty,&
     val )
  implicit none
  ! INPUT
  real(8) :: ar_x, ar_y
  real(8) :: bvalue,val
  integer  :: ai_nx, ai_kx, ai_ny, ai_ky,sz_tx,sz_ty
  real(8), dimension(:), pointer :: apr_tx !  ai_nx + ai_kx 
  real(8), dimension(:), pointer :: apr_ty !  ai_ny + ai_ky	
  real(8), dimension(:,:),pointer :: apr_Bcoef!( ai_nx,ai_ny)			
  ! LOCAL VARIABLES
  integer  :: li_i, li_j, li_mflag, li_lefty
  real(8), dimension(:),pointer :: lpr_coef ! ai_ky
  real(8), dimension(:),pointer :: tmp_tab
  real(8), dimension(:),pointer :: tmp_ty
  integer :: ierr
  
  call interv ( apr_ty, sz_ty, ar_y, li_lefty, li_mflag )
  
  if ( li_mflag .NE. 0 ) then
     val = 0.0_8
     return 
  end if

  ALLOCATE(lpr_coef(ai_ky))
  ALLOCATE(tmp_tab(ai_nx))
  ALLOCATE(tmp_ty( 2*ai_ky ))
  do li_j = 1, ai_ky
     
     
     tmp_tab = apr_bcoef ( 1:ai_nx , li_lefty - ai_ky + li_j )
     
     lpr_coef ( li_j ) = bvalue(&
          apr_tx,&
          tmp_tab,&
          ai_nx,&
          ai_kx,&
          ar_x,&
          0 )
     
     
  end do
  
  tmp_ty =  apr_ty ( li_lefty - ai_ky + 1 : li_lefty + ai_ky)
  val = bvalue(&
       tmp_ty,&
       lpr_coef,&
       ai_ky,&
       ai_ky,&
       ar_y,&
       0 )
  
  DEALLOCATE(lpr_coef)
  DEALLOCATE(tmp_tab)
  DEALLOCATE(tmp_ty)
end subroutine bvalue2d
