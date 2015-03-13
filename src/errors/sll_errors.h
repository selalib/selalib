#ifndef _SLL_ERRORS_H_
#define _SLL_ERRORS_H_

#define SLL_ERROR(fun,msg)   call sll_error_handler  (__FILE__,__LINE__,fun,msg)
#define SLL_WARNING(fun,msg) call sll_warning_handler(__FILE__,__LINE__,fun,msg)

use sll_errors

#endif
