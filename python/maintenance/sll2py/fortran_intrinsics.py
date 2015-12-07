intrinsic_types = { \
'integer',
'real',
'double precision',
'complex',
'character',
'logical',
}

intrinsic_procedures = { \
'abort', # Terminates program and produces core dump (GNU)
'abs', # Absolute value
'achar', # Character in ASCII collating sequence
'acos', # Arccosine function
'acosh', # Inverse hyperbolic cosine function
'adjustl', # Left adjust a string
'adjustr', # Right adjust a string
'aimag', # Imaginary part of complex number
'aint', # Truncate to a whole number
'all', # All values in MASK along DIM are true
'allocated', # Status of an allocatable entity
'anint', # Nearest whole number
'any', # Any value in MASK along DIM is true
'asin', # Arcsine function
'asinh', # Inverse hyperbolic sine function
'associated', # Status of a pointer or pointer/target pair
'atan', # Arctangent function
'atan2', # Arctangent function
'atanh', # Inverse hyperbolic tangent function
'bessel_j0', # Bessel function of the first kind of order 0
'bessel_j1', # Bessel function of the first kind of order 1
'bessel_jn', # Bessel function of the first kind
'bessel_y0', # Bessel function of the second kind of order 0
'bessel_y1', # Bessel function of the second kind of order 1
'bessel_yn', # Bessel function of the second kind
'bge', # Bitwise greater than or equal to
'bgt', # Bitwise greater than
'bit_size', # Bit size inquiry function
'ble', # Bitwise less than or equal to
'blt', # Bitwise less than
'btest', # Bit test function
#------------------------------------------------------------------------------
# Following symbols come from module 'iso_c_binding', they should not be here!
#------------------------------------------------------------------------------
# 'c_associated',    # Status of a C pointer
# 'c_funloc',        # Obtain the C address of a procedure
# 'c_f_procpointer', # Convert C into Fortran procedure pointer
# 'c_f_pointer',     # Convert C into Fortran pointer
# 'c_loc',           # Obtain the C address of an object
# 'c_sizeof',        # Size in bytes of an expression
#------------------------------------------------------------------------------
'ceiling', # Integer ceiling function
'char', # Character conversion function
'cmplx', # Complex conversion function
'command_argument_count', # Get number of command line arguments
'conjg', # Complex conjugate function
'cos', # Cosine function
'cosh', # Hyperbolic cosine function
'count', # Count function
'cpu_time', # CPU elapsed time in seconds
'cshift', # Circular shift elements of an array
'date_and_time', # Date and time subroutine
'dble', # Double conversion function
'digits', # Significant digits function
'dim', # Positive difference
'dot_product', # Dot product function
'dprod', # Double product function
'dshiftl', # Combined left shift
'dshiftr', # Combined right shift
'eoshift', # End-off shift elements of an array
'epsilon', # Epsilon function
'erf', # Error function
'erfc', # Complementary error function
'erfc_scaled', # Error function
'execute_command_line', # Execute a shell command
'exit', # Terminates program with given status (GNU)
'exp', # Exponential function
'exponent', # Exponent function
'extends_type_of', # Type extension inquiry
'float', # Convert integer to default real
'floor', # Integer floor function
'fraction', # Fractional part of the model representation
'gamma', # Gamma function
'get_command', # Get the entire command line
'get_command_argument', # Get command line arguments
'get_environment_variable', # Get an environmental variable
'huge', # Largest number of a kind
'hypot', # Euclidean distance function
'iachar', # Code in ASCII collating sequence
'iand', # Bitwise logical and
'ibclr', # Clear bit
'ibits', # Bit extraction
'ibset', # Set bit
'ichar', # Character-to-integer conversion function
'ieor', # Bitwise logical exclusive or
'index', # Position of a substring within a string
'int', # Convert to integer type
'ior', # Bitwise logical inclusive or
'is_iostat_end', # Test for end-of-file value
'is_iostat_eor', # Test for end-of-record value
'ishft', # Shift bits
'ishftc', # Shift bits circularly
'kind', # Kind of an entity
'lbound', # Lower dimension bounds of an array
'leadz', # Number of leading zero bits of an integer
'len', # Length of a character entity
'len_trim', # Length of a character entity without trailing blank characters
'lge', # Lexical greater than or equal
'lgt', # Lexical greater than
'lle', # Lexical less than or equal
'llt', # Lexical less than
'log', # Logarithm function
'log10', # Base 10 logarithm function
'log_gamma', # Logarithm of the Gamma function
'logical', # Convert to logical type
'maskl', # Left justified mask
'maskr', # Right justified mask
'matmul', # matrix multiplication
'max', # Maximum value of an argument list
'maxexponent', # Maximum exponent of a real kind
'maxloc', # Location of the maximum value within an array
'maxval', # Maximum value of an array
'merge', # Merge variables
'merge_bits', # Merge of bits under mask
'min', # Minimum value of an argument list
'minexponent', # Minimum exponent of a real kind
'minloc', # Location of the minimum value within an array
'minval', # Minimum value of an array
'mod', # Remainder function
'modulo', # Modulo function
'move_alloc', # Move allocation from one object to another
'mvbits', # Move bits from one integer to another
'nearest', # Nearest representable number
'new_line', # New line character
'nint', # Nearest whole number
'not', # Logical negation
'norm2', # Euclidean vector norm
'null', # Function that returns an disassociated pointer
'pack', # Pack an array into an array of rank one
'parity', # Reduction with exclusive or
'popcnt', # Number of bits set
'poppar', # Parity of the number of bits set
'precision', # Decimal precision of a real kind
'present', # Determine whether an optional dummy argument is specified
'product', # Product of array elements
'radix', # Base of a model number
'random_number', # Pseudo-random number
'random_seed', # Initialize a pseudo-random number sequence
'range', # Decimal exponent range of a real kind
'real', # Convert to real type
'repeat', # Repeated string concatenation
'reshape', # Function to reshape an array
'rrspacing', # Reciprocal of the relative spacing
'same_type_as', # Dynamic type inquiry
'scale', # Scale a real value
'scan', # Scan a string for the presence of a set of characters
'selected_char_kind', # Choose character kind
'selected_int_kind', # Choose integer kind
'selected_real_kind', # Choose real kind
'set_exponent', # Set the exponent of the model
'shape', # Determine the shape of an array
'shifta', # Right shift with fill
'shiftl', # Left shift
'shiftr', # Right shift
'sign', # Sign copying function
'sin', # Sine function
'sinh', # Hyperbolic sine function
'size', # Determine the size of an array
'sleep', # Causes the process to pause for Seconds seconds (GCC and ifort)
'sngl', # Convert double precision real to default real
'spacing', # Smallest distance between two numbers of a given type
'spread', # Add a dimension to an array
'sqrt', # Square-root function
'storage_size', # Storage size of argument A in bits (Fortran 2008 and later)
'sum', # Sum of array elements
'system', # Passes the command COMMAND to a shell (GCC and ifort)
          # F2008 'execute_command_line' exists (GCC 4.7+ or ifort 15+ needed)
'system_clock', # Time function
'tan', # Tangent function
'tanh', # Hyperbolic tangent function
'tiny', # Smallest positive number of a real kind
'trailz', # Number of trailing zero bits of an integer
'transfer', # Transfer bit patterns
'transpose', # Transpose an array of rank two
'trim', # Remove trailing blank characters of a string
'ubound', # Upper dimension bounds of an array
'unpack', # Store the elements of a vector in an array of higher rank
'verify', # Scan a string for the absence of a set of characters
}

logical_operators = { \
'.not.',
'.and.',
'.or.',
'.eqv.',
'.neqv.',
}

logical_constants = { \
'.true.',
'.false.',
}

relational_operators = { \
'.lt.',
'.le.',
'.eq.',
'.ne.',
'.gt.',
'.ge.',
}
