
module finite_volume2
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use cubic_non_uniform_splines
  !use utils
  implicit none
contains  


subroutine Compute_flux2(a1,a2,f,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
  use sll_constants
  implicit none
  
  
  sll_int,intent(in)::N_x1,N_x2,order
  !sll_real64,dimension(:,:), pointer :: chi,sigma
   sll_real64,dimension(N_x1+1,N_x2+1):: a1,a2,f,jac_array, Flux_x1,Flux_x2,flux
  sll_real64,dimension(N_x1+1+2+1,N_x2+1+2+1) ::chi,sigma
  sll_real64,dimension(N_x1+1)::abar_x1,abar_x1m
  sll_real   ,dimension(N_x1+1)::abar_x2,abar_x2m!E
  sll_real64,intent(in)::x1_min,x2_min,delta_x1,delta_x2!L
  sll_int::i1,i2,i1m2,i1m1, i2m2,i2m1,i1p1,i1p2,i2p1,i2p2,i1m3,i2m3
  sll_real64::x1,x2,Flux_p05,Flux_m05,alpha,beta,df,coef1,tmp!,eold,enew,dx2,tmp
  sll_int   :: test
    !compute chi(uij) and sigma(uij)
     alpha=(1._f64/delta_x1)
     beta =(1._f64/delta_x2)
     !alpha=real(N_x1,f64)
     !beta =real(N_x2,f64)
   do i1=1,N_x1+1
        do i2=1,N_x2+1   
           chi(i1,i2)   =f(i1,i2)/jac_array(i1,i2)
           sigma(i1,i2) =f(i1,i2)/jac_array(i1,i2)
       enddo
   enddo

if(order==3) then
   
    chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2) 
    chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)

  ! For fixed  i2 
 do i2=1,N_x2+1
     do i1=1,N_x1
       abar_x1(i1) = a1(i1,i2)
     enddo
     
     
     abar_x1(N_x1+1)=abar_x1(1)


     

     abar_x1m(1)=abar_x1(N_x1)
      do i1=2,N_x1+1
       abar_x1m(i1)=abar_x1(i1-1)
      enddo

  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m2=i1-2
       i1m1=i1-1

       i1p1=modulo(i1+1-1,N_x1)+1
       i1p2=modulo(i1+2-1,N_x1)+1
       i1m2=modulo(i1-2-1,N_x1)+1
       i1m1=modulo(i1-1-1,N_x1)+1

     
       if (i1m1 <=0) then 
         i1m1=i1m1+N_x1
       endif
       if (i1m2 <=0) then 
         i1m2=i1m2+N_x1     
       endif
        
       !if (i1p1 >=N_x1+1) then 
        ! i1p1=i1p1+N_x1+1
       !endif
       !if (i1p2 >=N_x1+1) then 
        ! i1p2=i1p2+N_x1+1
       !endif
    
       if(abar_x1(i1)>=0) then   
        Flux_p05=-(1._f64/6._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1p1,i2)! stencil i1-1,i1,i1+1  
        !Flux_m05=-(1._f64/6._f64)*chi(i1m2,i2)+(5._f64/6._f64)*chi(i1m1,i2)+(1._f64/3._f64)*chi(i1,i2)! stencil i1-1,i1,i1+1   
       else 
        Flux_p05=(1._f64/3._f64)*chi(i1,i2)+(5._f64/6._f64)*chi(i1p1,i2)-(1._f64/6._f64)*chi(i1p2,i2)! stencil i1,i1+1,i1+2
        !Flux_m05=(1._f64/3._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)-(1._f64/6._f64)*chi(i1p1,i2)! stencil i1,i1+1,i1+2
      endif  

       if(abar_x1m(i1)>=0) then   
        !Flux_p05=-(1._f64/6._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1p1,i2)! stencil i1-1,i1,i1+1  
        Flux_m05=-(1._f64/6._f64)*chi(i1m2,i2)+(5._f64/6._f64)*chi(i1m1,i2)+(1._f64/3._f64)*chi(i1,i2)! stencil i1-1,i1,i1+1   
       else 
        !Flux_p05=(1._f64/3._f64)*chi(i1,i2)+(5._f64/6._f64)*chi(i1p1,i2)-(1._f64/6._f64)*chi(i1p2,i2)! stencil i1,i1+1,i1+2
        Flux_m05=(1._f64/3._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)-(1._f64/6._f64)*chi(i1p1,i2)! stencil i1,i1+1,i1+2
      endif  
       
       Flux_x1(i1,i2)=Flux_p05!-Flux_m05
     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     
     do i2=1,N_x2
       abar_x2(i2) = a2(i1,i2)
     enddo

     
     
     abar_x2(N_x2+1)=abar_x2(1)!special uniform

      abar_x2m(1)=abar_x2(N_x2)
      do i2=2,N_x2+1
       abar_x2m(i2)=abar_x2(i2-1)
      enddo
     
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1

       i2p1=modulo(i2+1-1,N_x2)+1
       i2p2=modulo(i2+2-1,N_x2)+1
       i2m2=modulo(i2-2-1,N_x2)+1
       i2m1=modulo(i2-1-1,N_x2)+1

     
     
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2      
     endif
        
      ! if (i2p1 >=N_x2+1) then 
        ! i2p1=i2p1+N_x2+1
       !endif
       !if (i2p2 >=N_x2+1) then 
       !  i2p2=i2p2+N_x2+1
    !   endif
   
    
     if(abar_x2(i2)>=0) then   
        Flux_p05=-(1._f64/6._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  +(1._f64/3._f64)*sigma(i1,i2+1)! stencil i2-1,i2,i2+1  
        !Flux_m05=-(1._f64/6._f64)*sigma(i1,i2m2)+(5._f64/6._f64)*sigma(i1,i2m1)+(1._f64/3._f64)*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
     else 
        Flux_p05=(1._f64/3._f64)*sigma(i1,i2  )+(5._f64/6._f64)*sigma(i1,i2p1)-(1._f64/6._f64)*sigma(i1,i2p2)! stencil i2,i2+1,i2+2
        !Flux_m05=(1._f64/3._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  -(1._f64/6._f64)*sigma(i1,i2p1)! stencil i2,i2+1,i2+2
     endif
     
      if(abar_x2m(i2)>=0) then   
        !Flux_p05=-(1._f64/6._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  +(1._f64/3._f64)*sigma(i1,i2+1)! stencil i2-1,i2,i2+1  
        Flux_m05=-(1._f64/6._f64)*sigma(i1,i2m2)+(5._f64/6._f64)*sigma(i1,i2m1)+(1._f64/3._f64)*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
     else 
        !Flux_p05=(1._f64/3._f64)*sigma(i1,i2  )+(5._f64/6._f64)*sigma(i1,i2p1)-(1._f64/6._f64)*sigma(i1,i2p2)! stencil i2,i2+1,i2+2
        Flux_m05=(1._f64/3._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  -(1._f64/6._f64)*sigma(i1,i2p1)! stencil i2,i2+1,i2+2
     endif

      Flux_x2(i1,i2)=Flux_p05!-Flux_m05

  enddo
 enddo
endif!order 3








if(order==4) then
   
    chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2) 
    chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)

  ! For fixed  i2 
 do i2=1,N_x2+1
      x2 = x2_min+real(i2-1,f64)*delta_x2


  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m2=i1-2
       i1m1=i1-1

       i1p1=modulo(i1+1-1,N_x1)+1
       i1p2=modulo(i1+2-1,N_x1)+1
       i1m2=modulo(i1-2-1,N_x1)+1
       i1m1=modulo(i1-1-1,N_x1)+1

     
       if (i1m1 <=0) then 
         i1m1=i1m1+N_x1
       endif
       if (i1m2 <=0) then 
         i1m2=i1m2+N_x1     
       endif
        
       !if (i1p1 >=N_x1+1) then 
        ! i1p1=i1p1+N_x1+1
       !endif
       !if (i1p2 >=N_x1+1) then 
        ! i1p2=i1p2+N_x1+1
       !endif
    
        Flux_p05=-(1._f64/6._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1p1,i2)! stencil i1-1,i1,i1+1  
        tmp=(1._f64/3._f64)*chi(i1,i2)+(5._f64/6._f64)*chi(i1p1,i2)-(1._f64/6._f64)*chi(i1p2,i2)! stencil i1,i1+1,i1+2
        Flux_p05 = 0.5_f64*(Flux_p05+tmp)

        Flux_m05=-(1._f64/6._f64)*chi(i1m2,i2)+(5._f64/6._f64)*chi(i1m1,i2)+(1._f64/3._f64)*chi(i1,i2)! stencil i1-1,i1,i1+1   
        tmp=(1._f64/3._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)-(1._f64/6._f64)*chi(i1p1,i2)! stencil i1,i1+1,i1+2
        Flux_m05 = 0.5_f64*(Flux_m05+tmp)
        
        Flux_x1(i1,i2)=Flux_p05!-Flux_m05
     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1

       i2p1=modulo(i2+1-1,N_x2)+1
       i2p2=modulo(i2+2-1,N_x2)+1
       i2m2=modulo(i2-2-1,N_x2)+1
       i2m1=modulo(i2-1-1,N_x2)+1
     
     
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2      
     endif
        
      ! if (i2p1 >=N_x2+1) then 
        ! i2p1=i2p1+N_x2+1
       !endif
       !if (i2p2 >=N_x2+1) then 
       !  i2p2=i2p2+N_x2+1
    !   endif
   
    
        Flux_p05=-(1._f64/6._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  +(1._f64/3._f64)*sigma(i1,i2+1)! stencil i2-1,i2,i2+1  
        tmp=(1._f64/3._f64)*sigma(i1,i2  )+(5._f64/6._f64)*sigma(i1,i2p1)-(1._f64/6._f64)*sigma(i1,i2p2)! stencil i2,i2+1,i2+2
        Flux_p05 = 0.5_f64*(Flux_p05+tmp)
     
        Flux_m05=-(1._f64/6._f64)*sigma(i1,i2m2)+(5._f64/6._f64)*sigma(i1,i2m1)+(1._f64/3._f64)*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
        tmp=(1._f64/3._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  -(1._f64/6._f64)*sigma(i1,i2p1)! stencil i2,i2+1,i2+2
        Flux_m05 = 0.5_f64*(Flux_m05+tmp)

      Flux_x2(i1,i2)=Flux_p05!-Flux_m05

  enddo
 enddo
endif!order 4









 

if(order==5) then
     
    

    chi(N_x1+1,1:N_x2+4)=chi(1,1:N_x2+4)  
    chi(N_x1+2,1:N_x2+4)=chi(2,1:N_x2+4) 
    chi(N_x1+3,1:N_x2+4)=chi(3,1:N_x2+4)
    chi(N_x1+4,1:N_x2+4)=chi(4,1:N_x2+4) 

    !sigma(1:N_x1+4,N_x2+1)=sigma(1:N_x1+4,1)
    sigma(1:N_x1+4,N_x2+1)=0!sigma(1:N_x1+4,1)
    sigma(1:N_x1+4,N_x2+2)=0!sigma(1:N_x1+4,2)
    sigma(1:N_x1+4,N_x2+3)=0!sigma(1:N_x1+4,3)
    sigma(1:N_x1+4,N_x2+4)=0!sigma(1:N_x1+4,4)
    coef1=(9._f64/20._f64)
  ! For fixed  i2 
 do i2=1,N_x2+1
      x2 = x2_min+real(i2-1,f64)*delta_x2

     do i1=1,N_x1
       abar_x1(i1) = a1(i1,i2)
     enddo

     abar_x1(N_x1+1)=abar_x1(1)

      abar_x1m(1)=abar_x1(N_x1)
      do i1=2,N_x1+1
       abar_x1m(i1)=abar_x1(i1-1)
      enddo

  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m2=i1-2
       i1m1=i1-1
       i1m3=i1-3

       i1p1=modulo(i1+1-1,N_x1)+1
       i1p2=modulo(i1+2-1,N_x1)+1
       i1m2=modulo(i1-2-1,N_x1)+1
       i1m1=modulo(i1-1-1,N_x1)+1
       i1m3=modulo(i1-3-1,N_x1)+1


       if (i1m1 <=0) then 
         i1m1=i1m1+N_x1
       endif
       if (i1m2 <=0) then 
         i1m2=i1m2+N_x1     
       endif
        if (i1m3 <=0) then 
         i1m3=i1m3+N_x1     
       endif
     

       if(abar_x1(i1)>=0.) then   
        Flux_p05=(1._f64/30._f64)*chi(i1m2,i2)-(13._f64/60._f64)*chi(i1m1,i2)+(47._f64/60._f64)*chi(i1,i2)&
                +coef1*chi(i1+1,i2)-(1._f64/20._f64)*chi(i1+2,i2)! stencil i1-1,i1,i1+1 
        !Flux_m05=(1._f64/30._f64)*chi(i1m3,i2)-(13._f64/60._f64)*chi(i1m2,i2)+(47._f64/60._f64)*chi(i1m1,i2)&
            ! +coef1*chi(i1,i2)-(1._f64/20._f64)*chi(i1+1,i2)! stencil i1-1,i1,i1+1  
       else 
        Flux_p05=-(1._f64/20._f64)*chi(i1m1,i2)+ coef1*chi(i1,i2)+(47._f64/60._f64)*chi(i1+1,i2)&
                 -(13._f64/60._f64)*chi(i1+2,i2)+(1._f64/30._f64)*chi(i1+3,i2)! stencil i1,i1+1,i1+2
       ! Flux_m05=-(1._f64/20._f64)*chi(i1m2,i2)+ coef1*chi(i1m1,i2)+(47._f64/60._f64)*chi(i1,i2)-&
                 ! (13._f64/60._f64)*chi(i1+1,i2)+(1._f64/30._f64)*chi(i1+2,i2)!! stencil i1,i1+1,i1+2
       endif

       if(abar_x1m(i1)>=0.) then   
        !Flux_p05=(1._f64/30._f64)*chi(i1m2,i2)-(13._f64/60._f64)*chi(i1m1,i2)+(47._f64/60._f64)*chi(i1,i2)&
               ! +coef1*chi(i1+1,i2)-(1._f64/20._f64)*chi(i1+2,i2)! stencil i1-1,i1,i1+1 
        Flux_m05=(1._f64/30._f64)*chi(i1m3,i2)-(13._f64/60._f64)*chi(i1m2,i2)+(47._f64/60._f64)*chi(i1m1,i2)&
             +coef1*chi(i1,i2)-(1._f64/20._f64)*chi(i1+1,i2)! stencil i1-1,i1,i1+1  
       else 
        !Flux_p05=-(1._f64/20._f64)*chi(i1m1,i2)+ coef1*chi(i1,i2)+(47._f64/60._f64)*chi(i1+1,i2)&
                 !-(13._f64/60._f64)*chi(i1+2,i2)+(1._f64/30._f64)*chi(i1+3,i2)! stencil i1,i1+1,i1+2
        Flux_m05=-(1._f64/20._f64)*chi(i1m2,i2)+ coef1*chi(i1m1,i2)+(47._f64/60._f64)*chi(i1,i2)-&
                  (13._f64/60._f64)*chi(i1+1,i2)+(1._f64/30._f64)*chi(i1+2,i2)!! stencil i1,i1+1,i1+2
       endif

        Flux_x1(i1,i2)=Flux_p05!-Flux_m05
 

  enddo
 enddo 

  !For fixed i1
 do i1=1,N_x1+1

     do i2=1,N_x2
       abar_x2(i2) = a2(i1,i2)
     enddo


     abar_x2(N_x2+1)=abar_x2(1)!special uniform

      abar_x2m(1)=abar_x2(N_x2)
      do i2=2,N_x2+1
       abar_x2m(i2)=abar_x2(i2-1)
      enddo
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1
     i2m3=i2-3

       i2p1=modulo(i2+1-1,N_x2)+1
       i2p2=modulo(i2+2-1,N_x2)+1
       i2m2=modulo(i2-2-1,N_x2)+1
       i2m1=modulo(i2-1-1,N_x2)+1
       i2m3=modulo(i2-3-1,N_x2)+1

     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2      
     endif
        if (i2m3 <=0) then 
         i2m3=i2m3+N_x2     
       endif 
      ! if (i2p1 >=N_x2+1) then 
        ! i2p1=i2p1+N_x2+1
       !endif
       !if (i2p2 >=N_x2+1) then 
       !  i2p2=i2p2+N_x2+1
    !   endif
   
    
       if(abar_x2(i2)>=0) then   
        Flux_p05=(1._f64/30._f64)*sigma(i1,i2m2)-(13._f64/60._f64)*sigma(i1,i2m1)+(47._f64/60._f64)*sigma(i1,i2)+&
                coef1*sigma(i1,i2+1)-(1._f64/20._f64)*sigma(i1,i2+2)! stencil i1-1,i1,i1+1  
        !Flux_m05=(1._f64/30._f64)*sigma(i1,i2m3)-(13._f64/60._f64)*sigma(i1,i2m2)+(47._f64/60._f64)*sigma(i1,i2m1)+&
                ! coef1*sigma(i1,i2)-(1._f64/20._f64)*sigma(i1,i2+1)! stencil i1-1,i1,i1+1  
       else 
        Flux_p05=-(1._f64/20._f64)*sigma(i1,i2m1)+ coef1*sigma(i1,i2)+(47._f64/60._f64)*sigma(i1,i2+1)-&
                   (13._f64/60._f64)*sigma(i1,i2+2)+(1._f64/30._f64)*sigma(i1,i2+3)! stencil i1,i1+1,i1+2
        !Flux_m05=-(1._f64/20._f64)*sigma(i1,i2m2)+ coef1*sigma(i1,i2m1)+(47._f64/60._f64)*sigma(i1,i2)-&
                 !(13._f64/60._f64)*sigma(i1,i2+1)+(1._f64/30._f64)*sigma(i1,i2+2)!! stencil i1,i1+1,i1+2
      endif  

        if(abar_x2m(i2)>=0) then   
        !Flux_p05=(1._f64/30._f64)*sigma(i1,i2m2)-(13._f64/60._f64)*sigma(i1,i2m1)+(47._f64/60._f64)*sigma(i1,i2)+&
                !coef1*sigma(i1,i2+1)-(1._f64/20._f64)*sigma(i1,i2+2)! stencil i1-1,i1,i1+1  
        Flux_m05=(1._f64/30._f64)*sigma(i1,i2m3)-(13._f64/60._f64)*sigma(i1,i2m2)+(47._f64/60._f64)*sigma(i1,i2m1)+&
                 coef1*sigma(i1,i2)-(1._f64/20._f64)*sigma(i1,i2+1)! stencil i1-1,i1,i1+1  
       else 
        !Flux_p05=-(1._f64/20._f64)*sigma(i1,i2m1)+ coef1*sigma(i1,i2)+(47._f64/60._f64)*sigma(i1,i2+1)-&
                 !  (13._f64/60._f64)*sigma(i1,i2+2)+(1._f64/30._f64)*sigma(i1,i2+3)! stencil i1,i1+1,i1+2
        Flux_m05=-(1._f64/20._f64)*sigma(i1,i2m2)+ coef1*sigma(i1,i2m1)+(47._f64/60._f64)*sigma(i1,i2)-&
                 (13._f64/60._f64)*sigma(i1,i2+1)+(1._f64/30._f64)*sigma(i1,i2+2)!! stencil i1,i1+1,i1+2
      endif  
       Flux_x2(i1,i2)=Flux_p05!-Flux_m05
  enddo
 enddo
endif

!Order 6
 if(order==6) then
   
    chi(N_x1+2,1:N_x2+4)=chi(2,1:N_x2+4) 
    chi(N_x1+3,1:N_x2+4)=chi(3,1:N_x2+4)
    chi(N_x1+4,1:N_x2+4)=chi(4,1:N_x2+4)
    sigma(1:N_x1+4,N_x2+2)=sigma(1:N_x1+4,2)
    sigma(1:N_x1+4,N_x2+3)=sigma(1:N_x1+4,3)
    sigma(1:N_x1+4,N_x2+4)=sigma(1:N_x1+4,4)
    
    
    !chi(N_x1+1,1:N_x2+4)=chi(1,1:N_x2+4)  
    !chi(N_x1+2,1:N_x2+4)=chi(2,1:N_x2+4) 
    !chi(N_x1+3,1:N_x2+4)=chi(3,1:N_x2+4)
    !chi(N_x1+4,1:N_x2+4)=chi(4,1:N_x2+4) 

    !!sigma(1:N_x1+4,N_x2+1)=sigma(1:N_x1+4,1)
    !sigma(1:N_x1+4,N_x2+1)=0!sigma(1:N_x1+4,1)
    !sigma(1:N_x1+4,N_x2+2)=0!sigma(1:N_x1+4,2)
    !sigma(1:N_x1+4,N_x2+3)=0!sigma(1:N_x1+4,3)
    !sigma(1:N_x1+4,N_x2+4)=0!sigma(1:N_x1+4,4)

    
    
    
    
    coef1=(9._f64/20._f64)
  ! For fixed  i2 
 do i2=1,N_x2+1
      x2 = x2_min+real(i2-1,f64)*delta_x2


  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m2=i1-2
       i1m1=i1-1
       i1m3=i1-3
       
       i1p1=modulo(i1+1-1,N_x1)+1
       i1p2=modulo(i1+2-1,N_x1)+1
       i1m2=modulo(i1-2-1,N_x1)+1
       i1m1=modulo(i1-1-1,N_x1)+1
       i1m3=modulo(i1-3-1,N_x1)+1

       
       if (i1m1 <=0) then 
         i1m1=i1m1+N_x1
       endif
       if (i1m2 <=0) then 
         i1m2=i1m2+N_x1     
       endif
        if (i1m3 <=0) then 
         i1m3=i1m3+N_x1     
       endif
        
        !*******
         
        Flux_p05=(1._f64/30._f64)*chi(i1m2,i2)-(13._f64/60._f64)*chi(i1m1,i2)+(47._f64/60._f64)*chi(i1,i2)&
                +coef1*chi(i1+1,i2)-(1._f64/20._f64)*chi(i1+2,i2)! stencil i1-1,i1,i1+1 
        
       
        tmp=-(1._f64/20._f64)*chi(i1m1,i2)+ coef1*chi(i1,i2)+(47._f64/60._f64)*chi(i1+1,i2)&
                 -(13._f64/60._f64)*chi(i1+2,i2)+(1._f64/30._f64)*chi(i1+3,i2)! stencil i1,i1+1,i1+2      
        Flux_p05 = 0.5_f64*(Flux_p05+tmp)
        
         
        Flux_m05=(1._f64/30._f64)*chi(i1m3,i2)-(13._f64/60._f64)*chi(i1m2,i2)+(47._f64/60._f64)*chi(i1m1,i2)&
             +coef1*chi(i1,i2)-(1._f64/20._f64)*chi(i1+1,i2)! stencil i1-1,i1,i1+1  
       
        tmp=-(1._f64/20._f64)*chi(i1m2,i2)+ coef1*chi(i1m1,i2)+(47._f64/60._f64)*chi(i1,i2)-&
                  (13._f64/60._f64)*chi(i1+1,i2)+(1._f64/30._f64)*chi(i1+2,i2)!! stencil i1,i1+1,i1+2
      
        Flux_m05 = 0.5_f64*(Flux_m05+tmp)
        
        Flux_x1(i1,i2)=Flux_p05!-Flux_m05
 
        !******
    
      !  Flux_p05=-(1._f64/6._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1p1,i2)! stencil i1-1,i1,i1+1  
       ! tmp=(1._f64/3._f64)*chi(i1,i2)+(5._f64/6._f64)*chi(i1p1,i2)-(1._f64/6._f64)*chi(i1p2,i2)! stencil i1,i1+1,i1+2
       ! Flux_p05 = 0.5_f64*(Flux_p05+tmp)

       ! Flux_m05=-(1._f64/6._f64)*chi(i1m2,i2)+(5._f64/6._f64)*chi(i1m1,i2)+(1._f64/3._f64)*chi(i1,i2)! stencil i1-1,i1,i1+1   
       ! tmp=(1._f64/3._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)-(1._f64/6._f64)*chi(i1p1,i2)! stencil i1,i1+1,i1+2
       ! Flux_m05 = 0.5_f64*(Flux_m05+tmp)
        
       ! Flux_x1(i1,i2)=Flux_p05-Flux_m05
     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     
  do i2=1,N_x2+1
    
     !******
      i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1
     i2m3=i2-3
     

       i2p1=modulo(i2+1-1,N_x2)+1
       i2p2=modulo(i2+2-1,N_x2)+1
       i2m2=modulo(i2-2-1,N_x2)+1
       i2m1=modulo(i2-1-1,N_x2)+1
       i2m3=modulo(i2-3-1,N_x2)+1
     
     
     
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2      
     endif
        if (i2m3 <=0) then 
         i2m3=i2m3+N_x2     
       endif 
    
   
    
       
        Flux_p05=(1._f64/30._f64)*sigma(i1,i2m2)-(13._f64/60._f64)*sigma(i1,i2m1)+(47._f64/60._f64)*sigma(i1,i2)+&
                coef1*sigma(i1,i2+1)-(1._f64/20._f64)*sigma(i1,i2+2)! stencil i1-1,i1,i1+1  
        
        tmp=-(1._f64/20._f64)*sigma(i1,i2m1)+ coef1*sigma(i1,i2)+(47._f64/60._f64)*sigma(i1,i2+1)-&
                   (13._f64/60._f64)*sigma(i1,i2+2)+(1._f64/30._f64)*sigma(i1,i2+3)! stencil i1,i1+1,i1+2
        

         Flux_p05 = 0.5_f64*(Flux_p05+tmp)
         
        Flux_m05=(1._f64/30._f64)*sigma(i1,i2m3)-(13._f64/60._f64)*sigma(i1,i2m2)+(47._f64/60._f64)*sigma(i1,i2m1)+&
                 coef1*sigma(i1,i2)-(1._f64/20._f64)*sigma(i1,i2+1)! stencil i1-1,i1,i1+1  
       
        tmp=-(1._f64/20._f64)*sigma(i1,i2m2)+ coef1*sigma(i1,i2m1)+(47._f64/60._f64)*sigma(i1,i2)-&
                 (13._f64/60._f64)*sigma(i1,i2+1)+(1._f64/30._f64)*sigma(i1,i2+2)!! stencil i1,i1+1,i1+2
     
       Flux_m05 = 0.5_f64*(Flux_m05+tmp)
        
        Flux_x2(i1,i2)=Flux_p05!-Flux_m05
     !*****
   
    
       ! Flux_p05=-(1._f64/6._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  +(1._f64/3._f64)*sigma(i1,i2+1)! stencil i2-1,i2,i2+1  
        !tmp=(1._f64/3._f64)*sigma(i1,i2  )+(5._f64/6._f64)*sigma(i1,i2p1)-(1._f64/6._f64)*sigma(i1,i2p2)! stencil i2,i2+1,i2+2
       ! Flux_p05 = 0.5_f64*(Flux_p05+tmp)
     
        !Flux_m05=-(1._f64/6._f64)*sigma(i1,i2m2)+(5._f64/6._f64)*sigma(i1,i2m1)+(1._f64/3._f64)*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
        !tmp=(1._f64/3._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  -(1._f64/6._f64)*sigma(i1,i2p1)! stencil i2,i2+1,i2+2
        !Flux_m05 = 0.5_f64*(Flux_m05+tmp)

      !Flux_x2(i1,i2)=Flux_p05-Flux_m05

  enddo
 enddo
endif!order 4


  do i1=1,N_x1+1
     do i2=1,N_x2+1
       i1m1=modulo(i1-1-1+N_x1,N_x1)+1
       i2m1=modulo(i2-1-1+N_x2,N_x2)+1
       i1p1=modulo(i1+1-1+N_x1,N_x1)+1
       i2p1=modulo(i2+1-1+N_x2,N_x2)+1
       
       !Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
       Flux(i1,i2)=-(alpha*(a1(i1,i2)*Flux_x1(i1,i2)-a1(i1m1,i2)*Flux_x1(i1m1,i2))&
       +beta*(a2(i1,i2)*Flux_x2(i1,i2)-a2(i1,i2m1)*Flux_x2(i1,i2m1)))

       Flux(i1,i2)=Flux(i1,i2)&
       -(alpha*(1._f64/real(N_x2,f64))**2/48._f64*((a1(i1,i2p1)-a1(i1,i2m1))*(Flux_x1(i1,i2p1)-Flux_x1(i1,i2m1))&
       -(a1(i1m1,i2p1)-a1(i1m1,i2m1))*(Flux_x1(i1m1,i2p1)-Flux_x1(i1m1,i2m1)))&
       +beta*(1._f64/real(N_x1,f64))**2/48._f64*(&
       (a2(i1p1,i2)-a2(i1m1,i2))*(Flux_x2(i1p1,i2)-Flux_x2(i1m1,i2))&
       -(a2(i1p1,i2m1)-a2(i1m1,i2m1))*(Flux_x2(i1p1,i2m1)-Flux_x2(i1m1,i2))&
       ))


     enddo
   enddo






end subroutine compute_flux2




subroutine compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)
  use sll_constants
  implicit none

  sll_int,intent(in) :: nc_eta1,nc_eta2
  sll_real64,dimension(1:nc_eta1+1) :: rho
  sll_real64,dimension(1:nc_eta1+1) :: phi_poisson
  sll_real64,dimension(1:nc_eta1+1) :: E
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: x1n_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: x2n_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: jac_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: a1
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: a2
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: psi
  sll_real64,dimension(:), pointer :: buf_1d
  sll_real64 :: x1_min,x1_max,x2_min,x2_max
  sll_int :: i1,i2,ii,err,i1m1,i1p1,i1p2,i1m2,i1m3,i2m1,i2p1,i2p2,i2m2
  sll_real64 :: tmp,x1,x2,phi_val,xx,a(-10:10)
  sll_real64,intent(in) :: delta_eta1,delta_eta2
  sll_real64,intent(in) :: geom_x(2,2)
  sll_int,intent(in) :: div_case
  sll_int :: N
  
  N = max(nc_eta1,nc_eta2)
  
  SLL_ALLOCATE(buf_1d(N+1),err)
  
  !a[-2] = 0, a[-1] = -1/6, a[0] = 1, a[1] = -1/2, a[2] = -1/3
  x1_min = geom_x(1,1)
  x1_max = geom_x(2,1)
  x2_min = geom_x(1,2)
  x2_max = geom_x(2,2)



    E=rho-1._f64
    call poisson1dpertrap(E,x1_max-x1_min,nc_eta1)
    phi_poisson = E
    call poisson1dpertrap(phi_poisson,x1_max-x1_min,nc_eta1)
    tmp = phi_poisson(1)
    do i1=1,nc_eta1
      phi_poisson(i1) = -phi_poisson(i1) + tmp
    enddo
    phi_poisson(nc_eta1+1) = phi_poisson(1)
    
 do i1=2,nc_eta1+1
   i1m3 = modulo(i1-3-1+nc_eta1,nc_eta1)+1
   i1m2 = modulo(i1-2-1+nc_eta1,nc_eta1)+1
   i1m1 = modulo(i1-1-1+nc_eta1,nc_eta1)+1
   i1p1 = modulo(i1+1-1+nc_eta1,nc_eta1)+1
   i1p2 = modulo(i1+2-1+nc_eta1,nc_eta1)+1
   if(modulo(div_case,10)==0)then
     buf_1d(i1) = 0.5_f64*(phi_poisson(i1)+phi_poisson(i1-1))
   endif
   if(modulo(div_case,10)==1)then   
     buf_1d(i1)=(1._f64/3._f64)*phi_poisson(i1m1)+(5._f64/6._f64)*phi_poisson(i1)&
     -(1._f64/6._f64)*phi_poisson(i1p1)
   endif  
   if(modulo(div_case,10)==2)then   
     buf_1d(i1)=(1._f64/30._f64)*phi_poisson(i1m3)-(13._f64/60._f64)*phi_poisson(i1m2)&
     +(47._f64/60._f64)*phi_poisson(i1m1)+&
                 (9._f64/20._f64)*phi_poisson(i1)-(1._f64/20._f64)*phi_poisson(i1p1)
   endif  

 enddo
 buf_1d(1)=buf_1d(nc_eta1+1)
 phi_poisson(1:nc_eta1+1) = buf_1d(1:nc_eta1+1)

    
    do i1=1,nc_eta1+1
      do i2=1,nc_eta2+1
        x1 = x1n_array(i1,i2)
        x2 = x2n_array(i1,i2)
        phi_val = 0._f64
        xx = (x1-x1_min)/(x1_max-x1_min)
        if(xx<=0._f64)then
          xx = 0._f64
        endif
        if(xx>=1._f64)then
          xx = xx-1._f64!1._f64-1e-15_f64
        endif
        if(xx<=0._f64)then
          xx = 0._f64
        endif      
        xx = xx*real(nc_eta1,f64)
        ii = floor(xx)
        xx = xx-real(ii,f64)      
        !phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)     
        phi_val = (0.5_f64*xx**3-xx**2-0.5_f64*xx+1._f64)*phi_poisson(ii+1)+(-0.5_f64*xx**3+0.5_f64*xx**2+xx)*phi_poisson(ii+2)     
        phi_val = phi_val+(-1._f64/6._f64*xx**3+0.5_f64*xx**2-1._f64/3._f64*xx)*phi_poisson(modulo(ii-1+nc_eta1,nc_eta1)+1)&
        +((xx**3-xx)/6._f64)*phi_poisson(modulo(ii+3-1+nc_eta1,nc_eta1)+1)     
       psi( i1, i2 ) = ( 0.5_f64*x2**2+phi_val) 
      enddo  
    enddo
    
    !print *,psi(1,:)-psi(nc_eta1+1,:)
    !print *,psi(:,nc_eta2+1)-psi(:,1)
    psi(nc_eta1+1,:)=psi(1,:)
    psi(:,nc_eta2+1)=psi(:,1)
   
    !stop
    
   
       do i1=1,nc_eta1
         do i2=1,Nc_eta2

       i1m2=modulo(i1-2-1+nc_eta1,nc_eta1)+1
       i1m1=modulo(i1-1-1+nc_eta1,nc_eta1)+1
       i2m2=modulo(i2-2-1+nc_eta2,nc_eta2)+1
       i2m1=modulo(i2-1-1+nc_eta2,nc_eta2)+1
       i1p1=modulo(i1+1-1+nc_eta1,nc_eta1)+1
       i1p2=modulo(i1+2-1+nc_eta1,nc_eta1)+1
       i2p1=modulo(i2+1-1+nc_eta2,nc_eta2)+1
       i2p2=modulo(i2+2-1+nc_eta2,nc_eta2)+1

           
           
           !a1(i1,i2)=((psi(i1,i2+1)-psi(i1,modulo(i2-1+nc_eta2,nc_eta2)+1))&
           !/(delta_eta2))*(x1_max-x1_min)/(0.5_f64*(jac_array(i1,i2)+jac_array(i1+1,i2)))
           !a2(i1,i2)=-((psi(i1+1,i2)-psi(modulo(i1-1+nc_eta1,nc_eta1)+1,i2))&
           !/(delta_eta1))*(x2_max-x2_min)/(0.5_f64*(jac_array(i1,i2)+jac_array(i1,i2+1)))
           a1(i1,i2)=((psi(i1+1,i2+1)-psi(i1+1,i2))&
           /(delta_eta2))*(x1_max-x1_min)!/(0.5_f64*(jac_array(i1,i2)+jac_array(i1+1,i2)))
           a2(i1,i2)=-((psi(i1+1,i2+1)-psi(i1,i2+1))&
           /(delta_eta1))*(x2_max-x2_min)!/(0.5_f64*(jac_array(i1,i2)+jac_array(i1,i2+1)))


           if(div_case>=10)then
             a1(i1,i2)=1._f64/24._f64*(psi(i1+1,i2m1)-psi(i1+1,i2p2))
             a1(i1,i2)=a1(i1,i2)+9._f64/8._f64*(psi(i1+1,i2p1)-psi(i1+1,i2))
             a1(i1,i2)=a1(i1,i2)/delta_eta2*(x1_max-x1_min)
             a2(i1,i2)=1._f64/24._f64*(psi(i1m1,i2+1)-psi(i1p2,i2+1))
             a2(i1,i2)=a2(i1,i2)+9._f64/8._f64*(psi(i1p1,i2+1)-psi(i1,i2+1))           
             a2(i1,i2)=-a2(i1,i2)/delta_eta1*(x2_max-x2_min)
           endif  
         enddo
       enddo
    
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo

            
         
  SLL_DEALLOCATE(buf_1d,err)


  
end subroutine compute_psi2

end module finite_volume2




program bgk_fv2
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"

  use sll_constants
  use distribution_function
  use sll_diagnostics
  use sll_csl
  use sll_splines
  use contrib_rho_module
  use finite_volume2
  
  implicit none
  !external Compute_flux,compute_psi
  !external poisson1dpertrap
  sll_real64 :: x2_min,x2_max,x1_min,x1_max,x1,x2,delta_x1,delta_x2,tmp
  sll_real64 :: mu,xi,L,H
  sll_int    :: i,j,N_phi,err,N_x1,N_x2,i1,i2,N,nb_step
  LOGICAL :: ex
  sll_real64,dimension(:), pointer :: phi,node_positions_x1,node_positions_x2,phi_poisson
  sll_real64,dimension(:), pointer :: new_node_positions,buf_1d,rho,rho2,E,rho_exact,E_store
  sll_real64,dimension(:,:), pointer :: f,f_init,f_equil, f_store,f_tmp
  sll_real64 :: phi_val,delta_x1_phi,xx,dt,alpha,val,max_f
  sll_int :: ii,step,div_case
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2

  sll_int32 :: nc_eta1, nc_eta2
  sll_real64, dimension(:,:), pointer :: x1c_array, x2c_array, jac_array
  sll_real64, dimension(:,:), pointer :: x1n_array, x2n_array
  sll_real64, dimension(:,:), pointer :: a1,a2,flux,psi,K1,K2,K3,K4
   sll_real64 :: b1,b2,b3,b4, a21, a32,a42,a43,a41,a31
  sll_int  ::  rk,order
  sll_real64 :: eta1_min, eta1_max,  eta2_min, eta2_max, eta1, eta2, eta1c, eta2c
  sll_real64 :: delta_eta1, delta_eta2,alpha_mesh,coef1,tmp2
  sll_int  :: mesh_case,ierr,visu_step,test_case,rho_case
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: mesh
  type(sll_distribution_function_2D_t), pointer :: dist_func
  character(32), parameter  :: name = 'distribution_function'
  sll_real64,dimension(:,:,:), pointer :: integration_points
  sll_real64,dimension(:,:), pointer :: integration_points_val
  sll_real64 :: geom_eta(2,2),geom_x(2,2)
  character*80,str,str2

  mesh_case = 3
  alpha_mesh =1.e-2_f64 !0.1_f64!0._f64!
  visu_step = 100
  test_case = 4
  rho_case = 2
  div_case=12
  
   rk=4
   order=5
  N_x1 = 64!256
  N_x2 = 64!256
  dt = 0.02_f64!0.005_f64
  nb_step =3000!120000! 6000!6000
  
  N = max(N_x1,N_x2)
  
  print *,'#max(N_x1,N_x2)=',N
  
  
  x2_max = 10._f64
  x2_min = -x2_max
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_init(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_equil(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  SLL_ALLOCATE(new_node_positions(N+1),err)
  SLL_ALLOCATE(buf_1d(N+1),err)
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(rho2(N_x1+1),err)
  SLL_ALLOCATE(rho_exact(N_x1+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(phi_poisson(N_x1+1),err)
  SLL_ALLOCATE(integration_points(3,N_x1+1,N_x2+1),err)  
  SLL_ALLOCATE(integration_points_val(2,N_x2),err) 
  SLL_ALLOCATE(f_store(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(a1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(a2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(Flux(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(psi(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_tmp(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K3(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K4(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(E_store(N_x1+1),err) 
  
  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, PERIODIC_SPLINE)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, PERIODIC_SPLINE)
  
  f_equil=0._f64
  
  !physical parameters
  mu=0.92_f64
  xi=0.90_f64
  L=14.71_f64
  
  inquire(file='half_phi.dat', exist=ex) 
  if(.not.(ex))then
    print *,'file half_phi.dat does not exist'
    stop
  endif  
  open(unit=900,file='half_phi.dat')  
    read(900,*) N_phi,L
    N_phi = 2*N_phi
    SLL_ALLOCATE(phi(N_phi+1),err)
    do j=1,N_phi/2+1
      read(900,*) i,x1,x2
      phi(i)=x1
    enddo
    do j=N_phi/2+2,N_phi+1
      phi(j)=phi(N_phi+2-j)
    enddo
  close(900)
  
  !L = 4._f64*sll_pi
  
  x1_min = 0._f64
  x1_max = L
  
  
  if(test_case>=4)then
    L = 4._f64*sll_pi
    x1_min = 0._f64
    x1_max = L
    x2_min = -6._f64
    x2_max = -x2_min
  endif

!  x1_min = 0._f64
!  x1_max = 1._f64
!  x2_min = 0._f64
!  x2_max = 1._f64
  
  
  
  
  
  
  delta_x1_phi = (x1_max-x1_min)/real(N_phi,f64)
  
  !nb_step = 10
  !dt = L/real(nb_step,f64)
  

  open(unit=900,file='phi0.dat')  
    do i1=1,N_phi+1
      x1 = x1_min+real(i1-1,f64)*delta_x1_phi
      write(900,*) x1,phi(i1)
    enddo
  close(900)
  
  

  
  print *,'#N_phi=',N_phi
  
  delta_x1 = (x1_max-x1_min)/real(N_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(N_x2,f64)

  do i1=1,N_x1+1
    x1 = x1_min+real(i1-1,f64)*delta_x1
    node_positions_x1(i1) = x1
  enddo
  do i2=1,N_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    node_positions_x2(i2) = x2
  enddo
  
  
  !definition of mesh

  nc_eta1 = N_x1
  nc_eta2 = N_x2

  eta1_min =  0.0_f64
  eta1_max = 1.0_f64 ! 0.15_f64*x1_max! 1.0_f64
  eta2_min =  0.0_f64
  eta2_max =  1.0_f64

  !eta1_min =  -0.5_f64
  !eta1_max = 0.5_f64 ! 0.15_f64*x1_max! 1.0_f64
  !eta2_min =  -0.5_f64
  !eta2_max =  0.5_f64


  !eta1_min =  x1_min
  !eta1_max = x1_max
  !eta2_min = x2_min
  !eta2_max = x2_max



  delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
  delta_eta2 = (eta1_max-eta1_min)/real(nc_eta2,f64)
  SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x1c_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x2c_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), ierr)
  

  geom_x(1,1)=x1_min
  geom_x(2,1)=x1_max 
  geom_x(1,2)=x2_min 
  geom_x(2,2)=x2_max 

  
  if(mesh_case==1)then
    do i2=1,nc_eta2+1
      do i1=1,nc_eta1+1
        x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
        x2n_array(i1,i2) = x2_min+real(i2-1,f64)*delta_x2
        x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
        x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
        jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)
        !jac_array(i1,i2) = 1._f64!(x1_max-x1_min)*(x2_max-x2_min)
      enddo
    enddo
    geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    do i2=1,nc_eta2+1
      do i1=1,nc_eta1+1
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
      enddo
    enddo
    
    
    
  endif
  


  

  if(mesh_case==2)then
     eta2 = 0.0_f64 
     eta2c =  0.5_f64*delta_eta2
     do i2= 1, nc_eta2 + 1
        eta1 = 0.0_f64
        eta1c = 0.5_f64*delta_eta1
        do i1 = 1, nc_eta1 + 1
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x1c_array(i1,i2) = eta1c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
           x2c_array(i1,i2) = eta2c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
           !x1n_array(i1,i2) = (x1n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2n_array(i1,i2) = (x2n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x1c_array(i1,i2) = (x1c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2c_array(i1,i2) = (x2c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !jacobian defined on center
           jac_array(i1,i2) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)) * &
             (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c)) - &
             alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c) * &
             alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)
           !jacobian defined on nodes
           jac_array(i1,i2) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1) * sin (2*sll_pi*eta2)) * &
             (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1) * cos (2*sll_pi*eta2)) - &
             alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1) * cos (2*sll_pi*eta2) * &
             alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1) * sin (2*sll_pi*eta2)


           eta1 = eta1 + delta_eta1
           eta1c = eta1c + delta_eta1
           x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
           x1c_array(i1,i2) = x1_min+x1c_array(i1,i2)*(x1_max-x1_min)
           x2c_array(i1,i2) = x2_min+x2c_array(i1,i2)*(x2_max-x2_min)
           jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
        end do
        eta2 = eta2 + delta_eta2
        eta2c = eta2c + delta_eta2
     end do

    geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
    dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    val = 0._f64
    do i2=1,nc_eta2
      do i1=1,nc_eta1
       ! x1 = (real(i1,f64)-0.5_f64)/real(nc_eta1,f64)
         x1 = (real(i1,f64)-1._f64)/real(nc_eta1,f64)
        !eta2 = (real(i2,f64)-0.5_f64)/real(nc_eta2,f64)
        eta2 = (real(i2,f64)-1._f64)/real(nc_eta2,f64)
        tmp = alpha_mesh*sin(2._f64*sll_pi*eta2)
        do i=1,100
          val = val-(val+tmp*sin(2._f64*sll_pi*val)-x1)/&
          (1._f64+2._f64*sll_pi*tmp*cos(2._f64*sll_pi*val))
        enddo
        if(abs(val+tmp*sin(2._f64*sll_pi*val)-x1)>1.e-14)then
          print *,i1,i2,val+tmp*sin(2._f64*sll_pi*val)-x1,val
          print *,'Problem of convergence of Newton'
          stop
        endif
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = val
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+(x1-val+eta2)*(x2_max-x2_min)        
      enddo
    enddo
    !do i1=1,nc_eta1
    !  print *,i1,integration_points(i1,1), integration_points(i1,nc_eta2)
    !enddo
    !stop


  endif



 if(mesh_case==3)then
     eta2 = eta2_min 
     eta2c = eta2_min + 0.5_f64*delta_eta2
     do i2= 1, nc_eta2 + 1
        eta1 = eta1_min    
        eta1c = 0.5_f64*delta_eta1
        do i1 = 1, nc_eta1 + 1
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)**2
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x1c_array(i1,i2) = eta1c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)**2
           x2c_array(i1,i2) = eta2c + alpha_mesh * sin(2*sll_pi*eta1c)* sin(2*sll_pi*eta2c)
           !jacobian defined on center
           jac_array(i1,i2) = 1._f64+2._f64*sll_pi*alpha_mesh*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)&
           +2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1c)&
           -2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)**2&
           -4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)*cos(2._f64*sll_pi*eta1c)&
           +4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)**3*cos(2._f64*sll_pi*eta1c)

           !jacobian defined on nodes           
           !jac_array(i1,i2) = 1._f64+2._f64*sll_pi*alpha_mesh*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)&
           !+2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1)&
           !-2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)**2&
           !-4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)*cos(2._f64*sll_pi*eta1)&
           !+4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1)*cos(2._f64*sll_pi*eta2)**3*cos(2._f64*sll_pi*eta1)


           eta1 = eta1 + delta_eta1
           eta1c = eta1c + delta_eta1
           x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
           x1c_array(i1,i2) = x1_min+x1c_array(i1,i2)*(x1_max-x1_min)
           x2c_array(i1,i2) = x2_min+x2c_array(i1,i2)*(x2_max-x2_min)
           jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
        end do
        eta2 = eta2 + delta_eta2
        eta2c = eta2c + delta_eta2
     end do

    geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)




    dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    val = 0._f64
    do i2=1,nc_eta2
      do i1=1,nc_eta1
        x1 = x1_min+(real(i1,f64)-0.5_f64)/real(nc_eta1,f64)
        eta2  = (real(i2,f64)-0.5_f64)/real(nc_eta2,f64)
        !x1   = (real(i1,f64)-1._f64)/real(nc_eta1,f64)
        !eta2 = (real(i2,f64)-1._f64)/real(nc_eta2,f64)
        if(abs(eta2)<1e-14)eta2=1e-8;
        !x1   =x1_min+ (real(i1,f64))/real(nc_eta1,f64)
       ! eta2 =eta1_min+ (real(i2,f64))/real(nc_eta2,f64)
        tmp = alpha_mesh*sin(2._f64*sll_pi*eta2)**2
        do i=1,100
          val = val-(val+tmp*sin(2._f64*sll_pi*val)-x1)/&
          (1._f64+2._f64*sll_pi*tmp*cos(2._f64*sll_pi*val))
        enddo
        if(abs(val+tmp*sin(2._f64*sll_pi*val)-x1)>1.e-14)then
          print *,i1,i2,val+tmp*sin(2._f64*sll_pi*val)-x1,val
          print *,'Problem of convergence of Newton'
          stop
        endif
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = val
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+((x1-val)/sin(2._f64*sll_pi*eta2)+eta2)*(x2_max-x2_min)        
      enddo
    enddo
    
    !do i1=1,nc_eta1
       !do i2=1,nc_eta2
        !print*,i1,i2,integration_points(1,i1,i2), integration_points(2,i1,i2)
       !enddo
    !enddo
    !stop


  endif

  open(unit=900,file='intersect_points.dat')  
    do i1=1,N_x1
      x1 = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
      !x1 = x1_min+(real(i1,f64)-1._f64)*delta_x1
      do i2=1,N_x2
        !write(900,*) x1,integration_points(2,i1,i2),x1c_array(i1,i2),x2c_array(i1,i2)
        write(900,*) x1,integration_points(2,i1,i2),x1n_array(i1,i2),x2n_array(i1,i2)
      enddo  
    enddo
  close(900)
  
  

  
  !call sll_init_distribution_function_2D( dist_func, GAUSSIAN )
  do i1=1,nc_eta1+1
    do i2=1,nc_eta2+1
      x1 = x1c_array(i1,i2)
      x2 = x2c_array(i1,i2)
      phi_val = 0._f64
      xx = (x1-x1_min)/(x1_max-x1_min)
      if(xx<=0._f64)then
        xx = 0._f64
      endif
      if(xx>=1._f64)then
        xx = xx-1._f64!1._f64-1e-15_f64
      endif
      if(xx<=0._f64)then
        xx = 0._f64
      endif
      
      xx = xx*real(N_phi,f64)
      ii = floor(xx)
      xx = xx-real(ii,f64)      
      phi_val = (1._f64-xx)*phi(ii+1)+xx*phi(ii+2)
      !phi_val = 0._f64
      H = 0.5_f64*x2*x2 + phi_val
      if(test_case==1)then
        val = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
      endif
      if(test_case==2)then
        val = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
        val = val*(1._f64+0.1_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==3)then
        val = exp(-0.5_f64*40._f64*((x1-.5_f64)**2+(x2-.5_f64)**2))
      endif
      if(test_case==4)then
        !linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        f_equil(i1,i2) = val*jac_array(i1,i2)
        val = val*(1._f64+0.001_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==5)then
        !non linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.5_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==6.or.test_case==7)then
        !gaussian equilibrium
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        !f_equil(i1,i2) = val*jac_array(i1,i2)
      endif
      f_init(i1,i2) = val*jac_array(i1,i2)
      f(i1,i2) = val*jac_array(i1,i2)   
       !Flux( i1, i2 ) = ( 0.5_f64*x2**2+phi_val)   
      call sll_set_df_val(dist_func, i1, i2, f_init(i1,i2))      
    enddo
  enddo
   f_store=f
  do i1=1,N_x1+1
      xx = (real(i1,f64)-0.5_f64)/real(N_x1,f64)
      !xx = (real(i1,f64)-1._f64)/real(N_x1,f64)
      if(xx<=0._f64)then
        xx = 0._f64
      endif
      if(xx>=1._f64)then
        xx = xx-1._f64!1._f64-1e-15_f64
      endif
      if(xx<=0._f64)then
        xx = 0._f64
      endif      
      xx = xx*real(N_phi,f64)
      ii = floor(xx)
      xx = xx-real(ii,f64)      
      phi_val = (1._f64-xx)*phi(ii+1)+xx*phi(ii+2)
      rho_exact(i1) = mu*(3._f64-2._f64*xi+2*phi_val)/(3._f64-2._f64*xi)*exp(-phi_val)
  enddo
  
  
  call write_mesh_2D(mesh)
   
  call write_distribution_function ( dist_func )


  

  ! initialize CSL  
  !csl_work => new_csl_workspace( dist_func )
  !uniform_field => new_field_2D_vec1(mesh)
  !uniform_field_new => new_field_2D_vec1(mesh)
  !uniform_field_velocity => new_field_2D_vec1(mesh)

  geom_eta(1,1) = eta1_min
  geom_eta(2,1) = eta1_max
  geom_eta(1,2) = eta2_min
  geom_eta(2,2) = eta2_max

  call compute_rho_mapped_mesh2(rho,f_store,integration_points,rho_case,nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
  
  call compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
  geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)
  E_store=E
    
    
    
     
     !diagnostic
      do i1=1,nc_eta1+1
        x1 = x1_min+real(i1-0.5_f64)*delta_x1
        !write(*,*) x1,E(i1),rho(i1)!,rho_exact(i1)
      enddo

     !write(*,*) delta_eta1,delta_eta2,delta_x1, delta_x2,jac_array(1,2),(x1_max-x1_min)*(x2_max-x2_min)
    !stop
   f_equil=1._f64
     do i2=1,nc_eta2+1
      do i1 = 1,nc_eta1+2
      f_equil(i1,i2)=jac_array(i1,i2)
      enddo
     enddo
     
    
  do step = 1, nb_step

  if(rk==1) then
    f_tmp=f
    call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        f(i1,i2)=f_store(i1,i2)+dt*Flux(i1,i2)
     enddo
   enddo
  endif
 
    if(rk==2) then
   
        !Flux=(-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2))
        b1=0._f64
        b2=1._f64
        a21=1._f64/2._f64
        do i1=1,N_x1+1
          do i2=1,N_x2+1
        f_tmp(i1,i2)=f_store(i1,i2)!jac_array(i1,i2)
        enddo
          enddo
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K1=Flux
         
        do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a21
          enddo
         enddo

        call compute_rho_mapped_mesh2(rho,f_tmp,integration_points,rho_case,&
        nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
  
        call compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
        geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)

         
         
         
         
         
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K2=Flux
         
         
         
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        
        !K2(i1,i2)=0._f64/3._f64
        f(i1,i2)=(f_store(i1,i2)+dt*(b1*K1(i1,i2)+b2*K2(i1,i2)))!*jac_array(i1,i2)
      !-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
    enddo
    

  endif
    if(rk==4) then
   
        !Flux=(-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2))
        b1=1._f64/6._f64
        b2=1._f64/3._f64
        b3=1._f64/3._f64
        b4=1._f64/6._f64
       
        a21=1._f64/2._f64
        a32=1._f64/2._f64
        a31=0._f64
        a41=0._f64
        a42=0._f64
        a43=1._f64
        f_tmp=f_store
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K1=Flux
        do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a21
          enddo
         enddo
        
        !update field psi
        call compute_rho_mapped_mesh2(rho,f_tmp,integration_points,rho_case,nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
  
        call compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
        geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)

         
                
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K2=Flux

         do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a31+K2(i1,i2)*dt*a32
          enddo
         enddo
         
        !update field psi
        call compute_rho_mapped_mesh2(rho,f_tmp,integration_points,rho_case,&
        nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
  
        call compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
        geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)


         
         
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K3=Flux
         do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a41+K2(i1,i2)*dt*a42+K3(i1,i2)*dt*a43
          enddo
         enddo

        !update field psi
        call compute_rho_mapped_mesh2(rho,f_tmp,integration_points,rho_case,&
        nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
  
        call compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
        geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)

         
         
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K4=Flux
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        f(i1,i2)=f_store(i1,i2)+dt*(b1*K1(i1,i2)+b2*K2(i1,i2)+b3*K3(i1,i2)+b4*K4(i1,i2))
      !-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
    enddo
  endif
    f_store=f
 !diagnostic********************************
  
   if(rk==2) then
   
        !Flux=(-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2))
        b1=0._f64
        b2=1._f64
        a21=1._f64/2._f64
        do i1=1,N_x1+1
          do i2=1,N_x2+1
        f_tmp(i1,i2)=f_equil(i1,i2)!jac_array(i1,i2)
        enddo
          enddo
        call Compute_flux2(a1,a2,f_tmp,f_equil,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K1=Flux
        do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_equil(i1,i2)+K1(i1,i2)*dt*a21
          enddo
         enddo
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K2=Flux
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        
        !K2(i1,i2)=0._f64/3._f64
        f_init(i1,i2)=(f_equil(i1,i2)+dt*(b1*K1(i1,i2)+b2*K2(i1,i2)))!*jac_array(i1,i2)
      !-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
    enddo
  endif
    if(rk==4) then
   
        !Flux=(-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2))
        b1=1._f64/6._f64
        b2=1._f64/3._f64
        b3=1._f64/3._f64
        b4=1._f64/6._f64
       
        a21=1._f64/2._f64
        a32=1._f64/2._f64
        a31=0._f64
        a41=0._f64
        a42=0._f64
        a43=1._f64
        f_tmp=f_equil
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K1=Flux
        do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_equil(i1,i2)+K1(i1,i2)*dt*a21
          enddo
         enddo
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K2=Flux

         do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_equil(i1,i2)+K1(i1,i2)*dt*a31+K2(i1,i2)*dt*a32
          enddo
         enddo
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K3=Flux
         do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_equil(i1,i2)+K1(i1,i2)*dt*a41+K2(i1,i2)*dt*a42+K3(i1,i2)*dt*a43
          enddo
         enddo
        call Compute_flux2(a1,a2,f_tmp,jac_array,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K4=Flux
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        f_init(i1,i2)=f_equil(i1,i2)+dt*(b1*K1(i1,i2)+b2*K2(i1,i2)+b3*K3(i1,i2)+b4*K4(i1,i2))
      !-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
    enddo
  endif
   max_f = 0._f64 
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        if(f(i1,i2)>max_f)then
          max_f = f(i1,i2)
        endif
          if(f(i1,i2)<1.e-2)then
          f_init(i1,i2)=jac_array(i1,i2)
        endif
     enddo
   enddo
    f_equil=f_init
    
    
 !diagnostic********************************

        !update field psi
        call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,nc_eta1,&
        nc_eta2,geom_eta,jac_array,spl_per_x1)
  
        call compute_psi2(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
        geom_x,x1n_array,x2n_array,jac_array,delta_eta1,delta_eta2,div_case)


     
     
     E_store=E
      
    f_store = f
     !diagnostic
    !do i1=1,N_x1+1
    !  x1 = x1_min+real(i1-1,f64)*delta_x1
    !  write(900,*) x1,E(i1),rho(i1)
    !enddo
   

  
  
   if(test_case==7)then
     phi_poisson = 0._f64
   endif


    tmp=sum(rho(1:N_x1))*delta_x1
  
    val=0._f64
    do i1=1,N_x1
      val = val+E(i1)*E(i1)
    enddo
    val = val/real(N_x1,f64)
    
    tmp2=0._f64
    do i1=1,N_x1
      do i2=1,N_x2
        if(abs(f_init(i1,i2)/jac_array(i1,i2)-1._f64)>tmp2)then
          tmp2 = abs(f_init(i1,i2)/jac_array(i1,i2)-1._f64)
        endif
        !print *,i1,i2,f_init(i1,i2)/jac_array(i1,i2)
      enddo
    enddo
    
    
    
    
    print *,(real(step,f64)-1._f64)*dt,val,tmp/(x1_max-x1_min)-1._f64,tmp2,max_f
     
  enddo! time loop*********************
       stop
  !diagnostic
    do i2=1,nc_eta2+1
      do i1 = 1,nc_eta1+1 
      print*,i1,i2,f_init(i1,i2), jac_array(i1,i2),f(i1,i2)
      enddo
    enddo
  stop
    !do i1 = 1, nc_eta1 
     ! do i2 = 1, nc_eta2    
       ! call sll_set_df_val(dist_func, i1, i2,f(i1,i2))
      !  
    !enddo
    !if(mod(step,visu_step)==0)then
     ! call write_distribution_function ( dist_func )
    !endif
   ! if(mod(step,visu_step)==0)then
     ! call write_distribution_function ( dist_func )
    !endif
  !end do time loop
   
  !val = 0._f64
 ! do i1=1,nc_eta1+1
   ! do i2=1,nc_eta2+1
    !  val = max(abs(sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)&
    !  -f_init(i1,i2)/jac_array(i1,i2)),val)
   ! enddo
  !enddo  
  
  !print *,'#',nc_eta1,nc_eta2,dt,nb_step,val
  
  !open(unit=900,file='field_final.dat')  
    !do i1=1,N_x1+1
     ! x1 = x1_min+real(i1-1,f64)*delta_x1
     ! write(900,*) x1,E(i1),rho(i1)
   ! enddo
  !close(900)


  
end program




subroutine poisson1dpertrap(E,L,N)
  use sll_constants
  implicit none
  sll_int,intent(in)::N
  sll_real64,dimension(N+1),intent(inout)::E
  sll_real64,intent(in)::L
  sll_int::i
  sll_real64::eold,enew,dx2,tmp
  !ensures at first that Ein is of mean zero
  
  tmp=0._f64
  do i=1,N
    tmp=tmp+E(i)
  enddo
  tmp=-tmp/real(N,f64)
  do i=1,N
    E(i)=E(i)+tmp
  enddo
  
  dx2=0.5_f64*L/real(N,f64)
  eold=E(1)
  E(1)=0._f64
  tmp=0._f64
  do i=1,N-1
    enew=E(i+1)
    E(i+1)=E(i)+dx2*(eold+enew)
    tmp=tmp+E(i+1)
    eold=enew
  enddo
  tmp=-tmp/real(N,f64)
  do i=1,N
    E(i)=E(i)+tmp
  enddo
  E(N+1)=E(1)
end subroutine poisson1dpertrap

