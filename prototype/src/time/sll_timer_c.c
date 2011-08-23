#include <time.h>
#include <stdio.h>

#include "utils.h"

/* eventually put these in their own proper place */
#define STRING(x)      #x
#define XSTRNG(x) STRNG(x)

#define MALLOC(t,n)  ((t*)robust_malloc( sizeof(t)*(n),       \
		      "MALLOC("#t","#n") at "__FILE__"("XSTRNG(__LINE__)")"))

void *robust_malloc( size_t sz, char *descr ){
  void *p;
  if( sz==0 ) return NULL;  /* null memory request */
  p = malloc(sz);           /* sz = sizeof(t)*(n)   */
  if( p==NULL ) {
    if( descr != NULL ) {
      fprintf( stderr, "%s failed, exiting program (sz = %lu)\n",
	       descr, (unsigned long int) sz );
    } else {
      fprintf( stderr, "undescribed malloc failed, exiting program"
	       "(sz=%lu)\n", (unsigned long int) sz );
    }    
    exit(1);
  }
  return p;
}

struct timespec *set_time_mark_C(void) {
  struct timespec * time_ptr;
  time_ptr = MALLOC(struct timespec,1);
  Try(clock_gettime(CLOCK_MONOTONIC, time_ptr));
  return time_ptr;
}

struct timespec *reset_time_mark_C( struct timespec* mark ) {
  Try(clock_gettime(CLOCK_MONOTONIC, mark));
  return mark;
}

double time_elapsed_since_C(struct timespec* t0 ) {
  struct timespec t1;
  Try(clock_gettime(CLOCK_MONOTONIC, &t1));
  return timerdiff( &t1,t0 );
}

