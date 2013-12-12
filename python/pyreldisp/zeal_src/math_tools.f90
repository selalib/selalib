module math_tools_module
  use Function_Input_Module
  
  implicit none

  !**********************
  contains
  !**********************

  subroutine check_overflow(omega,fomega)
    !----------------------------------------------------------
    ! checks whether call to fdf at omega generates an overflow
    !----------------------------------------------------------
    complex(kind=dp), intent(in)  :: omega
    complex(kind=dp), intent(out) :: fomega
    complex(kind=dp)  :: f, df

    call fdf(omega,f,df)
    if (abs(f)>abs(df)) then
      fomega = f
    else
      fomega = df
    end if
    return
  end subroutine check_overflow
end module math_tools_module
