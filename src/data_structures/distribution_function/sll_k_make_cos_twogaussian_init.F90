    read(file_id, cos_twogaussian)
    allocate( params%kx(params%dims(1),params%n_cos) )
    allocate( params%modnum(params%dims(1),params%n_cos) )
    allocate( params%alpha(params%n_cos) );
    allocate( params%phase_shift(params%n_cos) )
    allocate( params%v_thermal(params%dims(2),2) )
    allocate( params%v_mean(params%dims(2),2) )
    allocate( params%normal(2) )
    allocate( params%delta(2) )
    
    params%n_gaussians = 2
    if ( params%n_cos == 1 ) then
       params%kx(:,1) = kx
       params%modnum(:,1) = modnum
    else
       params%modnum = 0.0_f64
       do j=1, params%n_cos
          params%kx(j,j) = kx(j)
          params%modnum(j,j) = modnum(j)
       end do
    end if
    params%alpha = alpha
    params%phase_shift=phase_shift
    params%v_thermal(:,1) = v_thermal_1
    params%v_mean(:,1) = v_mean_1
    params%v_thermal(:,2) = v_thermal_2
    params%v_mean(:,2) = v_mean_2
    params%delta(1) = delta
    params%delta(2) = 1.0_f64 - delta

    do j=1,2
       params%normal(j) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
            product(params%v_thermal(:,j)))

    end do
    
  end subroutine
