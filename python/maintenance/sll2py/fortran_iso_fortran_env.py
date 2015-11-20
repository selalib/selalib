named_constants = { \
'ATOMIC_INT_KIND', #   Default-kind integer constant to be used as kind parameter when defining integer variables used in atomic operations. (Fortran 2008 or later.)
'ATOMIC_LOGICAL_KIND', # Default-kind integer constant to be used as kind parameter when defining logical variables used in atomic operations. (Fortran 2008 or later.)
'CHARACTER_KINDS', #    Default-kind integer constant array of rank one containing the supported kind parameters of the CHARACTER type. (Fortran 2008 or later.)
'CHARACTER_STORAGE_SIZE', #    Size in bits of the character storage unit.
'ERROR_UNIT', #    Identifies the preconnected unit used for error reporting.
'FILE_STORAGE_SIZE', #   Size in bits of the file-storage unit.
'INPUT_UNIT', #    Identifies the preconnected unit identified by the asterisk (*) in READ statement.
'INT8', 
'INT16', 
'INT32', 
'INT64', #    Kind type parameters to specify an INTEGER type with a storage size of 16, 32, and 64 bits. It is negative if a target platform does not support the particular kind. (Fortran 2008 or later.)
'INTEGER_KINDS', #    Default-kind integer constant array of rank one containing the supported kind parameters of the INTEGER type. (Fortran 2008 or later.)
'IOSTAT_END', # The value assigned to the variable passed to the IOSTAT= specifier of an input/output statement if an end-of-file condition occurred.
'IOSTAT_EOR', #    The value assigned to the variable passed to the IOSTAT= specifier of an input/output statement if an end-of-record condition occurred.
'IOSTAT_INQUIRE_INTERNAL_UNIT', #    Scalar default-integer constant, used by INQUIRE for the IOSTAT= specifier to denote an that a unit number identifies an internal unit. (Fortran 2008 or later.)
'NUMERIC_STORAGE_SIZE', # The size in bits of the numeric storage unit.
'LOGICAL_KINDS', #    Default-kind integer constant array of rank one containing the supported kind parameters of the LOGICAL type. (Fortran 2008 or later.)
'OUTPUT_UNIT', #    Identifies the preconnected unit identified by the asterisk (*) in WRITE statement.
'REAL32', 
'REAL64', 
'REAL128', #    Kind type parameters to specify a REAL type with a storage size of 32, 64, and 128 bits. It is negative if a target platform does not support the particular kind. (Fortran 2008 or later.)
'REAL_KINDS', #    Default-kind integer constant array of rank one containing the supported kind parameters of the REAL type. (Fortran 2008 or later.)
'STAT_LOCKED',  #   Scalar default-integer constant used as STAT= return value by LOCK to denote that the lock variable is locked by the executing image. (Fortran 2008 or later.)
'STAT_LOCKED_OTHER_IMAGE', #    Scalar default-integer constant used as STAT= return value by UNLOCK to denote that the lock variable is locked by another image. (Fortran 2008 or later.)
'STAT_STOPPED_IMAGE', #    Positive, scalar default-integer constant used as STAT= return value if the argument in the statement requires synchronisation with an image, which has initiated the termination of the execution. (Fortran 2008 or later.)
'STAT_FAILED_IMAGE', #    Positive, scalar default-integer constant used as STAT= return value if the argument in the statement requires communication with an image, which has is in the failed state. (TS 18508 or later.)
'STAT_UNLOCKED', #    Scalar default-integer constant used as STAT= return value by UNLOCK to denote that the lock variable is unlocked. (Fortran 2008 or later.) 
}

derived_types = {  \
'LOCK_TYPE', #    Derived type with private components to be use with the LOCK and UNLOCK statement. A variable of its type has to be always declared as coarray and may not appear in a variable-definition context. (Fortran 2008 or later.) 
}

intrinsic_procedures = { \
'compiler_options',
'compiler_version'
}
