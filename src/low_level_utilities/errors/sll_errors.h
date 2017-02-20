#ifndef _SLL_ERRORS_H_
#define _SLL_ERRORS_H_

#ifdef __PGI

#define SLL_ERROR(fun,msg)   stop msg
#define SLL_WARNING(fun,msg) print*, msg

#else

#define SLL_ERROR(fun,msg)   call sll_s_error_handler  (__FILE__,__LINE__,fun,msg)
#define SLL_WARNING(fun,msg) call sll_s_warning_handler(__FILE__,__LINE__,fun,msg)

#endif

use sll_m_errors

#endif
