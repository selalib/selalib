module conservative_finite_difference
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use cubic_non_uniform_splines
  !use utils
  implicit none

contains 
subroutine Compute_flux(a1,a2,f,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
  use sll_constants
  implicit none
  
  
  sll_int,intent(in)::N_x1,N_x2,order
  !sll_real64,dimension(:,:), pointer :: chi,sigma
   sll_real64,dimension(N_x1+1,N_x2+1):: a1,a2,f,f_store, Flux_x1,Flux_x2,flux
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
   do i1=1,N_x1+1
        do i2=1,N_x2+1   
           chi(i1,i2)   =a1(i1,i2)*f(i1,i2)
           sigma(i1,i2) =a2(i1,i2)*f(i1,i2)
       enddo
   enddo
   
   
if(order==2) then

    chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2) 
    chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)

  ! For fixed  i2 
 do i2=1,N_x2+1

     
     
 
  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m2=i1-2
       i1m1=i1-1
     
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
    
        Flux_p05=0.5_f64*chi(i1,i2)+0.5_f64*chi(i1p1,i2)! stencil i1-1,i1,i1+1  

        Flux_m05=0.5_f64*chi(i1m1,i2)+0.5_f64*chi(i1,i2)! stencil i1-1,i1,i1+1   

       
       Flux_x1(i1,i2)=Flux_p05-Flux_m05
     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     

     
     
     
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1
     
     
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2      
     endif
        
   
    
        Flux_p05=0.5_f64*sigma(i1,i2)+0.5_f64*sigma(i1,i2+1)! stencil i2-1,i2,i2+1  
 
     
        Flux_m05=0.5_f64*sigma(i1,i2m1)+0.5_f64*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
   
      Flux_x2(i1,i2)=Flux_p05-Flux_m05

  enddo
 enddo
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
     enddo
   enddo



endif   

if(order==3) then
   
    chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2) 
    chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)

  ! For fixed  i2 
 do i2=1,N_x2+1
      x2 = x2_min+real(i2-1,f64)*delta_x2
     !do i1=1,N_x1
     !abar_x1(i1)=(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))!a1(i1,i2)!special uniform     
     !enddo
     do i1=1,N_x1
      df=abs(f_store(i1+1,i2)-f_store(i1,i2))
     if(df>1.d-14) then
     abar_x1(i1)=(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))!a1(i1,i2)!special uniform
      else
     abar_x1(i1)=chi(i1+1,i2)-chi(i1,i2)
     endif
     enddo
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
       
       Flux_x1(i1,i2)=Flux_p05-Flux_m05
     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     do i2=1,N_x2
      df=abs(f_store(i1,i2+1)-f_store(i1,i2))
     if(df>1.d-14) then
      abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     else
      abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))
     endif
     !abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))!a1(i1,i2)!(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     enddo
     
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

      Flux_x2(i1,i2)=Flux_p05-Flux_m05

  enddo
 enddo
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
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
        
        Flux_x1(i1,i2)=Flux_p05-Flux_m05
     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1
     
     
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

      Flux_x2(i1,i2)=Flux_p05-Flux_m05

  enddo
 enddo
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
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
     !do i1=1,N_x1
     !abar_x1(i1)=(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))!a1(i1,i2)!special uniform     
     !enddo
     do i1=1,N_x1
      df=abs(f_store(i1+1,i2)-f_store(i1,i2))
     if(df>1.d-14) then
     abar_x1(i1)=(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))!a1(i1,i2)!special uniform
      else
     abar_x1(i1)=chi(i1+1,i2)-chi(i1,i2)
     endif
     enddo

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

        Flux_x1(i1,i2)=Flux_p05-Flux_m05
 

  enddo
 enddo 

  !For fixed i1
 do i1=1,N_x1+1
     do i2=1,N_x2
      df=abs(f_store(i1,i2+1)-f_store(i1,i2))
     if(df>1.d-14) then
      abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     else
      abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))
     endif
     !abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))!a1(i1,i2)!(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     enddo

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
       Flux_x2(i1,i2)=Flux_p05-Flux_m05
  enddo
 enddo
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
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
        
        Flux_x1(i1,i2)=Flux_p05-Flux_m05
 
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
        
        Flux_x2(i1,i2)=Flux_p05-Flux_m05
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
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
     enddo
   enddo
endif!order 4








end subroutine compute_flux



subroutine compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
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
  sll_real64 :: x1_min,x1_max,x2_min,x2_max
  sll_int :: i1,i2,ii
  sll_real64 :: tmp,x1,x2,phi_val,xx,a(-10:10)
  sll_real64,intent(in) :: delta_eta1,delta_eta2
  sll_real64,intent(in) :: geom_x(2,2)
  sll_int,intent(in) :: div_case
  
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
        phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)     
       psi( i1, i2 ) = ( 0.5_f64*x2**2+phi_val)!& utilisation de tableau abusive 
      enddo  
    enddo
   
     if(div_case==0)then
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
           a1(i1,i2)=((psi(i1,i2+1)-psi(i1,modulo(i2-1-1+nc_eta2,nc_eta2)+1))/(2._f64*delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           a2(i1,i2)=-((psi(i1+1,i2)-psi(modulo(i1-1-1+nc_eta1,nc_eta1)+1,i2))/(2._f64*delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif
     
     if(div_case==1)then
       a(-1) =-1._f64/3._f64
       a(0) =-1._f64/2._f64
       a(1) =1._f64
       a(2)=-1._f64/6._f64
       a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif
     if(div_case==2)then
       a(-2)= 1._f64/6._f64
       a(-1) =-1._f64
       a(0) =1._f64/2._f64
       a(1) =1._f64/3._f64
       a(2) = 0._f64
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
             a1(i1,i2)=a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1) + &
             a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
             a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
             a2(i1,i2)=-(a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif


     if(div_case==3)then
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
           a1(i1,i2)=((psi(i1,i2+1)-psi(i1,modulo(i2-1-1+nc_eta2,nc_eta2)+1))/(2._f64*delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           a2(i1,i2)=-((psi(i1+1,i2)-psi(modulo(i1-1-1+nc_eta1,nc_eta1)+1,i2))/(2._f64*delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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



       a(-2)= 1._f64/6._f64
       a(-1) =-1._f64
       a(0) =1._f64/2._f64
       a(1) =1._f64/3._f64
       a(2) = 0._f64
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            if(1+0*a1(i1,i2)>0._f64)then
            a1(i1,i2)=a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            endif
            if(1+0*a2(i1,i2)>0._f64)then  
            a2(i1,i2)=-(a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
            endif 
         enddo
       enddo


       a(-1) =-1._f64/3._f64
       a(0) =-1._f64/2._f64
       a(1) =1._f64
       a(2)=-1._f64/6._f64
       a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
          if(a1(i1,i2)<0._f64)then
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
             a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           endif     
          if(a2(i1,i2)<0._f64)then
            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif



     if(div_case==30)then



       a(-2)= 1._f64/6._f64
       a(-1) =-1._f64
       a(0) =1._f64/2._f64
       a(1) =1._f64/3._f64
       a(2) = 0._f64
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=-(a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          enddo
       enddo


       a(-1) =-1._f64/3._f64
       a(0) =-1._f64/2._f64
       a(1) =1._f64
       a(2)=-1._f64/6._f64
       a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            tmp=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
             tmp=(tmp/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a1(i1,i2)=0.5_f64*(a1(i1,i2)+tmp)
            tmp=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             tmp=(tmp/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
             a2(i1,i2)=0.5_f64*(a2(i1,i2)+tmp)
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

    
     endif



    if(div_case==4)then
       a(-3) = -1._f64/30._f64
       a(-2) = 1._f64/4._f64
       a(-1) = -1._f64
       a(0)  = 1._f64/3._f64
       a(1)  = 1._f64/2._f64
       a(2)  = -1._f64/20._f64
       a(3)  = 0._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)

            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif
      if(div_case==5)then
       a(-3) = 0._f64
       a(-2) = 1._f64/20._f64
       a(-1) = -1._f64/2._f64
       a(0)  = -1._f64/3._f64
       a(1)  = 1._f64
       a(2)  = -1._f64/4._f64
       a(3)  = 1._f64/30._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)

            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif



      if(div_case==50)then
       a(-3) = 0._f64
       a(-2) = 1._f64/20._f64
       a(-1) = -1._f64/2._f64
       a(0)  = -1._f64/3._f64
       a(1)  = 1._f64
       a(2)  = -1._f64/4._f64
       a(3)  = 1._f64/30._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)

            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo

       a(-3) = -1._f64/30._f64
       a(-2) = 1._f64/4._f64
       a(-1) = -1._f64
       a(0)  = 1._f64/3._f64
       a(1)  = 1._f64/2._f64
       a(2)  = -1._f64/20._f64
       a(3)  = 0._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            tmp=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)
            a1(i1,i2)=0.5_f64*(a1(i1,i2)+tmp)
            
            tmp=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             a2(i1,i2)=0.5_f64*(a2(i1,i2)+tmp)
         enddo
       enddo




        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
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
     endif

            
         



  
end subroutine compute_psi


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



end module conservative_finite_difference
