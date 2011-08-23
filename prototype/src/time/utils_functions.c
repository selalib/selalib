#include "utils.h"

void MYerror(char * where,int line,char *what) {
  fprintf(stderr,SETCOLOR_FAILURE"%s,%d:%s" SETCOLOR_NORMAL,
	  where,line,what);
  fflush(stderr);
  fflush(stdout);
  fsync(STDERR_FILENO);
  fsync(STDOUT_FILENO);
  ++global_error;
}

int MYperror(char * where,int line,char *what) {
  fprintf(stderr,SETCOLOR_FAILURE"%s,%d:",where,line);
  perror(what);
  fprintf(stderr,SETCOLOR_NORMAL);
  fflush(stderr);
  fflush(stdout);
  fsync(STDERR_FILENO);
  fsync(STDOUT_FILENO);
  ++global_error;
  by_now();
}

int MY_c_perror(int e,char * where,int line,char *what) {
  char *fmt;
  fmt = e ? "expected" :SETCOLOR_FAILURE"UNEXPECTED";
  if ( e != errno || detail) {
    fprintf(stderr,"%s %s,%d:",fmt,where,line);
    perror(what);
  }
  if (e != errno){
    errno = e;
    perror(SETCOLOR_FAILURE"Expected");
    global_error++;
    fprintf(stderr,SETCOLOR_NORMAL);
  }
  fflush(stderr);
  fflush(stdout);
  fsync(STDERR_FILENO);
  fsync(STDOUT_FILENO);
  return e ? 0:-1;
}

int No_error(int e,char * where,int line,char *what) {
  char *fmt;
  fmt = e ? SETCOLOR_FAILURE"ERROR expected but not found" : "Cool";
  if (e || detail) {
    fprintf(stderr,"%s %s,%d:%s\n",
	    fmt,where,line,what);
  }
  if (e){
    errno = 0;
    perror("Found   ");
    errno = e;
    perror("Expected");
    global_error++;
    fprintf(stderr,SETCOLOR_NORMAL);
  }
  fflush(stderr);
  fflush(stdout);
  fsync(STDERR_FILENO);
  fsync(STDOUT_FILENO);
  return e ? -1 : 0;
}

void printfx(char *s,int d, int x){
  printf(s,d/x,d%x);
}

void printf2(char *s,int d, int e){ 
  printf(s,d/1000,d%1000,e/1000,e%1000);
}


