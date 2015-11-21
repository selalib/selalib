data_types = [ \
'S',  # real
'D',  # double precision
'C',  # complex
'Z',  # complex*16 or double complex
]

level1 = [ \
'xROTG' ,  # setup Givens rotation
'xROTMG',  # setup modified Givens rotation
'xROT'  ,  # apply Givens rotation
'xROTM' ,  # apply modified Givens rotation
'xSWAP' ,  # swap x and y
'xSCAL' ,  # x = a*x
'xCOPY' ,  # copy x into y
'xAXPY' ,  # y = a*x + y
'xDOT'  ,  # dot product
'xDSDOT',  # dot product with extended precision accumulation
'xNRM2' ,  # Euclidean norm
'xCNRM2',  # Euclidean norm
'xASUM' ,  # sum of absolute values
'IxAMAX',  # index of max abs value
]

level2 = [ \
'xGEMV',  # matrix vector multiply
'xGBMV',  # banded matrix vector multiply
'xSYMV',  # symmetric matrix vector multiply
'xSBMV',  # symmetric banded matrix vector multiply
'xSPMV',  # symmetric packed matrix vector multiply
'xTRMV',  # triangular matrix vector multiply
'xTBMV',  # triangular banded matrix vector multiply
'xTPMV',  # triangular packed matrix vector multiply
'xTRSV',  # solving triangular matrix problems
'xTBSV',  # solving triangular banded matrix problems
'xTPSV',  # solving triangular packed matrix problems
'xGER' ,  # performs the rank 1 operation A := alpha*x*y' + A
'xSYR' ,  # performs the symmetric rank 1 operation A := alpha*x*x' + A
'xSPR' ,  # symmetric packed rank 1 operation A := alpha*x*x' + A
'xSYR2',  # performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
'xSPR2',  # performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
]

level3 = [ \
'xGEMM',  # matrix matrix multiply
'xSYMM',  # symmetric matrix matrix multiply
'xHEMM',  # hermitian matrix matrix multiply
'xSYRK',  # symmetric rank-k update to a matrix
'xHERK',  # hermitian rank-k update to a matrix
'xSYR2K', # symmetric rank-2k update to a matrix
'xHER2K', # hermitian rank-2k update to a matrix
'xTRMM',  # triangular matrix matrix multiply
'xTRSM',  # solving triangular matrix with multiple right hand sides
]

# Real routines
s_level1 = [r.replace('x','S') for r in level1]
s_level2 = [r.replace('x','S') for r in level2]
s_level3 = [r.replace('x','S') for r in level3 if 'HE' not in r]

# Double precision routines
d_level1 = [r.replace('x','D') for r in level1]
d_level2 = [r.replace('x','D') for r in level2]
d_level3 = [r.replace('x','D') for r in level3 if 'HE' not in r]

# Complex routines
c_level1 = [r.replace('x','C') for r in level1]
c_level2 = [r.replace('x','C') for r in level2]
c_level3 = [r.replace('x','C') for r in level3]

# Complex*16 (or double complex) routines
z_level1 = [r.replace('x','Z') for r in level1]
z_level2 = [r.replace('x','Z') for r in level2]
z_level3 = [r.replace('x','Z') for r in level3]

# All routines
all_routines = \
        s_level1 + s_level2 + s_level3 + \
        d_level1 + d_level2 + d_level3 + \
        c_level1 + c_level2 + c_level3 + \
        z_level1 + z_level2 + z_level3
