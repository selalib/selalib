#ifndef _sll_accumulators_h_
#define _sll_accumulators_h_

  use sll_accumulators

  ! Note that this macro depends on having the precision module active in its
  ! environment.
  
  ! tmp1 and tmp2 are just auxiliary double precision reals. The choice of the
  ! do-exit construct is to emulate C's do... while(0) idiom, thus a call like
  ! 
  ! if(condition)
  !   SLL_ACCUMULATE_PARTICLE_CHARGE(...)
  ! else
  !   ...
  ! endif
  !
  ! will be legitimate. Else syntax errors may occur. We can't expect the 
  ! laudable discipline of always using 'then' after the conditional 
  ! expression unfortunately.
#define SLL_ACCUMULATE_CHARGE( acc, dx, dy, icell, q, tmp1, tmp2) \
  do; \
    tmp1 = (1.0_f64 - dx); \
    tmp2 = (1.0_f64 - dy); \
    acc%q_acc(icell)%q_sw = acc%q_acc(icell)%q_sw + q*tmp1*tmp2; \
    acc%q_acc(icell)%q_se = acc%q_acc(icell)%q_se + q*  dx*tmp2; \
    acc%q_acc(icell)%q_nw = acc%q_acc(icell)%q_nw + q*tmp1*  dy; \
    acc%q_acc(icell)%q_ne = acc%q_acc(icell)%q_ne + q*  dx*  dy; \
    exit; \
 end do

#define SLL_ACCUMULATE_PARTICLE_CHARGE(acc,p,tmp1,tmp2) \
  SLL_ACCUMULATE_CHARGE(acc,p%dx,p%dy,p%ic,p%q,tmp1,tmp2)

#define SLL_ACCUMULATE_CHARGE_CS( acc, dx, dy, icell, q, tmp, temp) \
  do; \
    tmp(1,1) = dx**3;\
    tmp(1,2) = dy**3;\
    tmp(2,1) = (1.0_f64 - dx)**3; \
    tmp(2,2) = (1.0_f64 - dy)**3; \
    tmp(3,1) = (1.0_f64 + dx)**3 - 4.0_f64*tmp(1,1); \
    tmp(3,2) = (1.0_f64 + dy)**3 - 4.0_f64*tmp(1,2); \
    tmp(4,1) = (2.0_f64 - dx)**3 - 4.0_f64*tmp(2,1); \
    tmp(4,2) = (2.0_f64 - dy)**3 - 4.0_f64*tmp(2,2); \
    temp = q/36._f64; \
    acc%q_acc(icell)%q_im1j = acc%q_acc(icell)%q_im1j + temp*tmp(2,1)*tmp(4,2); \
    acc%q_acc(icell)%q_ij   = acc%q_acc(icell)%q_ij   + temp*tmp(4,1)*tmp(4,2); \
    acc%q_acc(icell)%q_ip1j = acc%q_acc(icell)%q_ip1j + temp*tmp(3,1)*tmp(4,2); \
    acc%q_acc(icell)%q_ip2j = acc%q_acc(icell)%q_ip2j + temp*tmp(1,1)*tmp(4,2); \
    acc%q_acc(icell)%q_im1jm1 = acc%q_acc(icell)%q_im1jm1 + temp*tmp(2,1)*tmp(2,2); \
    acc%q_acc(icell)%q_ijm1   = acc%q_acc(icell)%q_ijm1   + temp*tmp(4,1)*tmp(2,2); \
    acc%q_acc(icell)%q_ip1jm1 = acc%q_acc(icell)%q_ip1jm1 + temp*tmp(3,1)*tmp(2,2); \
    acc%q_acc(icell)%q_ip2jm1 = acc%q_acc(icell)%q_ip2jm1 + temp*tmp(1,1)*tmp(2,2); \
    acc%q_acc(icell)%q_im1jp1 = acc%q_acc(icell)%q_im1jp1 + temp*tmp(2,1)*tmp(3,2); \
    acc%q_acc(icell)%q_ijp1   = acc%q_acc(icell)%q_ijp1   + temp*tmp(4,1)*tmp(3,2); \
    acc%q_acc(icell)%q_ip1jp1 = acc%q_acc(icell)%q_ip1jp1 + temp*tmp(3,1)*tmp(3,2); \
    acc%q_acc(icell)%q_ip2jp1 = acc%q_acc(icell)%q_ip2jp1 + temp*tmp(1,1)*tmp(3,2); \
    acc%q_acc(icell)%q_im1jp2 = acc%q_acc(icell)%q_im1jp2 + temp*tmp(2,1)*tmp(1,2); \
    acc%q_acc(icell)%q_ijp2   = acc%q_acc(icell)%q_ijp2   + temp*tmp(4,1)*tmp(1,2); \
    acc%q_acc(icell)%q_ip1jp2 = acc%q_acc(icell)%q_ip1jp2 + temp*tmp(3,1)*tmp(1,2); \
    acc%q_acc(icell)%q_ip2jp2 = acc%q_acc(icell)%q_ip2jp2 + temp*tmp(1,1)*tmp(1,2); \
    exit; \
 end do

#define SLL_ACCUMULATE_PARTICLE_CHARGE_CS(acc,p,tmp,temp) \
 SLL_ACCUMULATE_CHARGE_CS(acc,p%dx,p%dy,p%ic,p%q,tmp,temp)

#define SLL_GET_CHARGE_ACC_VALUE(acc,i,j,slot) \
 acc%q_acc(i+(j-1)*acc%mesh%num_cells1)%slot

#define SLL_INTERPOLATE_FIELD(part_Ex, part_Ey, cell_E, p, tmp3, tmp4) \
  do; \
    tmp3 = (1.0_f64 - p%dx); \
    tmp4 = (1.0_f64 - p%dy); \
    part_Ex = p%dx*p%dy*cell_E(p%ic)%Ex_sw + tmp3*p%dy*cell_E(p%ic)%Ex_se + p%dx*tmp4*cell_E(p%ic)%Ex_nw + tmp3*tmp4*cell_E(p%ic)%Ex_ne; \
    part_Ey = p%dx*p%dy*cell_E(p%ic)%Ey_sw + tmp3*p%dy*cell_E(p%ic)%Ey_se + p%dx*tmp4*cell_E(p%ic)%Ey_nw + tmp3*tmp4*cell_E(p%ic)%Ey_ne; \
    exit; \
 end do

#define SLL_INTERPOLATE_FIELD_CS(part_Ex, part_Ey, cellE, p, tmp) \
 do; \
    tmp(1,1) = p%dx**3; \
    tmp(1,2) = p%dy**3; \
    tmp(2,1) = (1.0_f64 - p%dx)**3; \
    tmp(2,2) = (1.0_f64 - p%dy)**3; \
    tmp(3,1) = (1.0_f64 + p%dx)**3 - 4.0_f64*tmp(1,1); \
    tmp(3,2) = (1.0_f64 + p%dy)**3 - 4.0_f64*tmp(1,2); \
    tmp(4,1) = (2.0_f64 - p%dx)**3 - 4.0_f64*tmp(2,1); \
    tmp(4,2) = (2.0_f64 - p%dy)**3 - 4.0_f64*tmp(2,2); \
    part_Ex = tmp(2,1)*tmp(4,2)*cellE(p%ic)%Ex_im1j + tmp(4,1)*tmp(4,2)*cellE(p%ic)%Ex_ij   + tmp(3,1)*tmp(4,2)*cellE(p%ic)%Ex_ip1j   + tmp(1,1)*tmp(4,2)*cellE(p%ic)%Ex_ip2j +\
            tmp(2,1)*tmp(2,2)*cellE(p%ic)%Ex_im1jm1 + tmp(4,1)*tmp(2,2)*cellE(p%ic)%Ex_ijm1 + tmp(3,1)*tmp(2,2)*cellE(p%ic)%Ex_ip1jm1 + tmp(1,1)*tmp(2,2)*cellE(p%ic)%Ex_ip2jm1 +\
            tmp(2,1)*tmp(3,2)*cellE(p%ic)%Ex_im1jp1 + tmp(4,1)*tmp(3,2)*cellE(p%ic)%Ex_ijp1 + tmp(3,1)*tmp(3,2)*cellE(p%ic)%Ex_ip1jp1 + tmp(1,1)*tmp(3,2)*cellE(p%ic)%Ex_ip2jp1 +\
            tmp(2,1)*tmp(1,2)*cellE(p%ic)%Ex_im1jp2 + tmp(4,1)*tmp(1,2)*cellE(p%ic)%Ex_ijp2 + tmp(3,1)*tmp(1,2)*cellE(p%ic)%Ex_ip1jp2 + tmp(1,1)*tmp(1,2)*cellE(p%ic)%Ex_ip2jp2; \
    part_Ex = part_Ex/36.0_f64; \
    part_Ey = tmp(2,1)*tmp(4,2)*cellE(p%ic)%Ey_im1j + tmp(4,1)*tmp(4,2)*cellE(p%ic)%Ey_ij   + tmp(3,1)*tmp(4,2)*cellE(p%ic)%Ey_ip1j   + tmp(1,1)*tmp(4,2)*cellE(p%ic)%Ey_ip2j +\
            tmp(2,1)*tmp(2,2)*cellE(p%ic)%Ey_im1jm1 + tmp(4,1)*tmp(2,2)*cellE(p%ic)%Ey_ijm1 + tmp(3,1)*tmp(2,2)*cellE(p%ic)%Ey_ip1jm1 + tmp(1,1)*tmp(2,2)*cellE(p%ic)%Ey_ip2jm1 +\
            tmp(2,1)*tmp(3,2)*cellE(p%ic)%Ey_im1jp1 + tmp(4,1)*tmp(3,2)*cellE(p%ic)%Ey_ijp1 + tmp(3,1)*tmp(3,2)*cellE(p%ic)%Ey_ip1jp1 + tmp(1,1)*tmp(3,2)*cellE(p%ic)%Ey_ip2jp1 +\
            tmp(2,1)*tmp(1,2)*cellE(p%ic)%Ey_im1jp2 + tmp(4,1)*tmp(1,2)*cellE(p%ic)%Ey_ijp2 + tmp(3,1)*tmp(1,2)*cellE(p%ic)%Ey_ip1jp2 + tmp(1,1)*tmp(1,2)*cellE(p%ic)%Ey_ip2jp2; \
    part_Ey = part_Ey/36.0_f64; \
    exit; \
 end do

#endif
