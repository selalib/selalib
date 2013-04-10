/**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
**************************************************************/
 
/*--------------------------------------------------------------
! file : ctime.c
! date : 08/30/2007
!  used to compute the CPU time
!  (without problem of negative error when the time is too long)
!---------------------------------------------------------------*/

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#define MSRATE 1000000LL
#define MAXSEC 1000000000LL

/*
 --------------------------------------------------------------------------  
   subroutine CLOCK_GET
 --------------------------------------------------------------------------  
*/
void
#ifndef NOUNDER
  clockget_ (long long int *dt) {
#else
  clockget (long long int *dt) {
#endif
  struct timeval tv2;
  unsigned long long int tvtmp;

  gettimeofday(&tv2, (struct timezone*)0);
  tvtmp = (((tv2).tv_sec)%MAXSEC  + 0.000001*(tv2).tv_usec) * MSRATE ;
  *dt   = tvtmp;
}


/*
 --------------------------------------------------------------------------  
   subroutine CLOCK_RATE
 --------------------------------------------------------------------------  
*/
void
#ifndef NOUNDER
  clockrate_ (long long int *rate) {
#else
  clockrate (long long int *rate) {
#endif
    *rate = MSRATE;
}


/*
 --------------------------------------------------------------------------  
   subroutine CLOCKMAXPERIOD
 --------------------------------------------------------------------------  
*/
void
#ifndef NOUNDER
  clockmaxperiod_ (long long int *periodm) {
#else
  clockmaxperiod (long long int *periodm) {
#endif
    *periodm = MSRATE * MAXSEC;
}

