  read(file_id, cos_onegaussian)
  
  allocate( params%kx(params%dims(1),params%n_cos) )
  allocate( params%alpha(params%n_cos) )   
  allocate( params%v_thermal(params%dims(2),1) )
  allocate( params%v_mean(params%dims(2),1) )
  allocate( params%normal(1) )
  allocate( params%delta(1) )

  params%n_gaussians = 1
  if ( params%n_cos == 1 ) then
     params%kx(:,1) = kx
  else
     params%kx = 0.0_f64
     do j=1, params%n_cos
        params%kx(j,j) = kx(j)
     end do
  end if
  params%alpha = alpha
  params%v_thermal(:,1) = v_thermal
  params%v_mean(:,1) = v_mean
  
  params%normal = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
       product(params%v_thermal(:,1)))
  
  params%delta(1) = 1.0_f64
  
end subroutine
