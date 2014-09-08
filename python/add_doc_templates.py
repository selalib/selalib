import os

root = "../prototype/src/"

def write_doc_file(library):
     _f = open(root+library+"/%s_doc.F90"%(library), "w")
     _f.write("! This file is read by doxygen software\n")
     _f.write("! Change it to match with your library\n")
     _f.write("! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks\n")
     _f.write("! Type 'make doc' in build directory\n")
     _f.write("! To check the results, open : \n")
     _f.write("! selalib/prototype/documentation/build/html/doxygen/html/namespaces.html \n")
     _f.write("! The following lines will be read by doxygen to generate documentation:\n\n\n")
     _f.write("!> @namespace sll_%s \n"%(library))
     _f.write("!> @brief \n")
     _f.write("!> Description of %s library (72 characters)\n"%(library))
     _f.write("!> @author Selalib team \n")
     _f.write("!> You can add a contact, do not put your email to prevent spam.")
     _f.write("!> @details\n")
     _f.write("!> Long description of  %s, you can add some references or math equations.\n"%(library))
     _f.write("!>\n")
     _f.write("!> <b> Headers file available </b>\n")
     _f.write("!>  - sll_%s.h\n"%(library))
     _f.write("!>\n")
     _f.write("!> <b> Modules available </b>\n")
     _f.write("!>  List fortran module available\n")
     _f.write("!>  - sll_%s\n"%(library))
     _f.write("!>\n")
     _f.write("!> <b> How to use it </b>\n")
     _f.write("!> - Header file : \code #include 'sll_%s.h' \endcode\n"%(library))
     _f.write("!> - Link with   <code>-lsll_%s</code>\n")
     _f.write("!> - Add <code> use sll_%s </code>\n"%(library))
     _f.write("!>\n")
     _f.write("!> <b> Examples </b>\n")
     _f.write("!> -Add some fortran lines to explain how ti use the library\n")
     _f.write("!> \code\n")
     _f.write("!> call initialize(my_type, arg_1, arg_2)\n")
     _f.write("!> call solve(my_type, your_result)\n")
     _f.write("!> \endcode\n")
     _f.write("!>\n")
     _f.close()



#for dir,subdir,name in os.walk(root):
#    if os.path.join(dir, name):
#        print dir

for dir in os.walk(root).next()[1]:
    if (dir != "CMakeModules" and dir != "package"):
        print root+dir
        mylib = dir.lower()
        write_doc_file(mylib)  