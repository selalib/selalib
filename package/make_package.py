# coding: utf8
"""
Repackage contents of all Selalib static libraries into a single archive
named 'libselalib.a'.

Modules required
----------------
  * Built-in  : os, shutil

"""
#
# Author: Yaman Güçlü, Jan 2016 - IPP Garching
#
# Last revision: 13 Jan 2016
#
from __future__ import print_function
import os
import shutil

#===============================================================================

# Identify path for package directory, and create new directory if missing
# TODO: PACKAGE_DIR should be passed as argument from CMake to Python
#
PWD = os.path.abspath( os.path.curdir )
BUILD_DIR   = PWD
LIBS_ROOT   = os.path.join( BUILD_DIR, "src" )
PACKAGE_DIR = os.path.join( BUILD_DIR, "package" )

if not os.path.exists( LIBS_ROOT ):
    print( "ERROR: 'LIBS_ROOT' path does not exist: %s" % LIBS_ROOT )
    raise SystemExit()
elif not os.path.isdir( LIBS_ROOT ):
    print( "ERROR: 'LIBS_ROOT' is not a directory: %s" % LIBS_ROOT )
    raise SystemExit()

#if not os.path.isdir( PACKAGE_DIR ):
#    os.mkdir( PACKAGE_DIR )

# Collect the paths to all library files in the build directory
list_lib = []
for dirpath, dirnames, filenames in os.walk( LIBS_ROOT ):
    for fname in filenames:
        root, ext = os.path.splitext( fname )
        if root.startswith( "libsll_" ) and ext == ".a":
            fpath = os.path.join( dirpath, fname )
            list_lib.append( fpath )
list_lib.sort()

# Extract all object files into a temporary directory
TMP_DIR = os.path.join( PACKAGE_DIR, "tmp" )

if os.path.isdir( TMP_DIR ):
    shutil.rmtree( TMP_DIR )

os.mkdir( TMP_DIR )
os.chdir( TMP_DIR )

for fpath in list_lib:
    command = "ar -x " + fpath
    exit_status = os.system( command )
    assert exit_status == 0

# Repackage all object files into a single archive (static library)
command = "ar -qc libselalib.a *.o"
exit_status = os.system( command )
assert exit_status == 0

# Move library to proper build directory, and remove temporary dir
shutil.copy( "libselalib.a", PACKAGE_DIR )
shutil.rmtree( TMP_DIR )
os.chdir( BUILD_DIR )
