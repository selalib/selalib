derived_types = { \
'C_PTR',    # Pointer to variable
'C_FUNPTR', # Pointer to function
}

intrinsic_procedures = { \
'C_ASSOCIATED', 
'C_F_POINTER',
'C_F_PROCPOINTER',
'C_FUNLOC',
'C_LOC',
'C_SIZEOF',
}

named_constants = { \
'C_int', 
'C_SHORT', # 	short int
'C_LONG', # 	long int
'C_LONG_LONG', # 	long long int
'C_SIGNED_CHAR', # 	signed char/unsigned char
'C_SIZE_T', # 	size_t
'C_INT8_T', # 	int8_t
'C_INT16_T', # 	int16_t
'C_INT32_T', # 	int32_t
'C_INT64_T ', #	int64_t
'C_INT128_T', # 	int128_t 	Ext.
'C_INT_LEAST8_T', # 	int_least8_t
'C_INT_LEAST16_T', # 	int_least16_t
'C_INT_LEAST32_T', # 	int_least32_t
'C_INT_LEAST64_T', # 	int_least64_t
'C_INT_LEAST128_T', #	int_least128_t 	Ext.
'C_INT_FAST8_T', # 	int_fast8_t
'C_INT_FAST16_T', # 	int_fast16_t
'C_INT_FAST32_T', # 	int_fast32_t
'C_INT_FAST64_T', # 	int_fast64_t
'C_INT_FAST128_T', # 	int_fast128_t 	Ext.
'C_INTMAX_T', # 	intmax_t
'C_INTPTR_T', #	intptr_t
'C_PTRDIFF_T', # 	intptr_t 	TS 29113
'C_FLOAT', #	float
'C_DOUBLE', # 	double
'C_LONG_DOUBLE', # 	long double
'C_FLOAT128', # 	__float128 	Ext.
'C_FLOAT_COMPLEX', # 	float _Complex
'C_DOUBLE_COMPLEX', #	double _Complex
'C_LONG_DOUBLE_COMPLEX', #	long double _Complex
'C_FLOAT128_COMPLEX', # 	__float128 _Complex 	Ext.
'C_BOOL', #	_Bool
'C_CHAR', #	char
'C_NULL_PTR', #  	C_PTR
'C_NULL_FUNPTR', #  	C_FUNPTR
}

parameters = { \
'C_NULL_CHAR', # 	null character 	'\0'
'C_ALERT', # 	alert 	'\a'
'C_BACKSPACE', # 	backspace 	'\b'
'C_FORM_FEED', # 	form feed 	'\f'
'C_NEW_LINE', # 	new line 	'\n'
'C_CARRIAGE_RETURN', # 	carriage return 	'\r'
'C_HORIZONTAL_TAB', # 	horizontal tab 	'\t'
'C_VERTICAL_TAB', # 	vertical tab 	'\v' 
}

