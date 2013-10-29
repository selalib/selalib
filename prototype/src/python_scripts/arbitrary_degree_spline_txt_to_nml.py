#!/usr/bin/env python3
"""
SYNOPSIS

    arbitrary_degree_spline_txt_to_nml [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    This script translates a .txt input file with information useful
    to specify a 2D arbitrary degree spline into a namelist format readable
    by Selalib's arbitrary degree spline 2D interpolator. To use, call the
    script with the .txt file as the argument. The output will be a .nml file
    with the same name as the original file. The original file should have 
    contents formatted as:

# degree
2, 1
# shape
9, 2
# rational
1
# knots
 0.    0.    0.    0.25  0.25  0.5   0.5   0.75  0.75  1.    1.    1.  
# knots
 0.  0.  1.  1.
# points
 0.  -0.5  0. 
 0. -1.  0.
-0.5 -0.5  0. 
-1. -1.  0.
-0.5  0.   0. 
-1.  0.  0.
-0.5  0.5  0. 
-1.  1.  0.
 0.   0.5  0. 
 0.  1.  0.
 0.5  0.5  0. 
 1.  1.  0.
 0.5  0.   0. 
 1.  0.  0.
 0.5 -0.5  0. 
 1. -1.  0.
 0.  -0.5  0. 
 0. -1.  0.
# weights
1.0
1.0
0.707106781187
0.707106781187
1.0
1.0
0.707106781187
0.707106781187
1.0
1.0
0.707106781187
0.707106781187
1.0
1.0
0.707106781187
0.707106781187
1.0
1.0

EXAMPLES

    ./arbitrary_degree_spline_txt_to_nml.py input_file.txt

EXIT STATUS

    TODO:

AUTHOR

    Edwin CHACON-GOLCHER <golcher@math.unistra.fr>

LICENSE

    Same as Selalib's...

VERSION

    1.0
"""

import sys, os, traceback, optparse
import time
import re
#from pexpect import run, spawn

def main ():

    global options, args
    # TODO: Do something more interesting here...
    inputname      = "" # empty filename at the start
    inputfilename  = "" # name of the read-only file (could differ from input)
    readfilename   = ""
    currently_reading = ""
    num_slots = 0
    knots1_is_written = False
    x1 = []
    x2 = []
    wgts = []
    size_x1 = 0
    size_x2 = 0
    point_counter = 0
    logical_mesh_locations1 = []
    logical_mesh_locations2 = []
    nc1 = 0  # number cells for logical mesh
    nc2 = 0  # number cells for logical mesh
    print('number of arguments passed')
    print(len(args))
    # in case that no arguments were provided, request a filename from the user
    if (len(args) == 0):
        print( "Enter the name of the file to convert to namelist format:")
        while 1:
            next = sys.stdin.read(1)
            if next == "\n":
                break
            inputname += next
    elif (len(args) > 1): # 2 or more arguments given
        print( "Please enter one argument only. Usage: ")
        print( "user$ ./arb_deg_txt_to_namelist.py filename.txt")
        sys.exit()
    else:                   # exactly one argument given
        inputname = args[0]
    # check whether the user has given the .txt extension or not, and create 
    # the name of the output file. Echo to screen the names of the files to be
    # read and written.
    numdots = inputname.count('.')
    if numdots == 0: # no extension, thus just add the extension to files
        outputname   = inputname + ".nml"   # output file name
        readfilename = inputname + ".txt" 
    elif numdots == 1: # there is an extension
        dotposition  = inputname.find('.')
        #   print inputname[dotposition:]
        if inputname[dotposition:] == ".txt": # it has the right extension
            readfilename = inputname           # open file with name as given
            outputname   = inputname[:dotposition] + ".nml" # create output name
        else:
            print( "Wrong extension. Only .txt files are allowed.")
            sys.exit()
    else:
        print( "That is a very weird-looking filename!")
        outputname   = inputname + ".nml"   # just add the extension
        readfilename = inputname + ".txt"  

    print( "The file to be processed is: {0} ".format(readfilename))
    print( "Converting {0} to {1}".format(readfilename, outputname))

    with open(readfilename,'r') as readfile, open(outputname,'w') as writefile:
        now  = time.localtime()
        date = str(now[1]) + "/" + str(now[2]) + "/" + str(now[0])
        mytime = str(now[3]) + ":" + str(now[4]) + ":" + str(now[5]) + "\n"
        writefile.write("! Input namelist file intended to initialize a ")
        writefile.write("2D coordinate transformation.\n")
        writefile.write("! This should be done with a call like: \n\n")
        writefile.write("!       call T%initialize_from_file(filename)\n\n")
        writefile.write("! Generated by "+sys.argv[0]+"\n")
        writefile.write("! on: " + date + "," + mytime)
        writefile.write("! Original input file: " + "\n" + "! "+ readfilename + "\n\n\n")

        linelist = readfile.readlines()

        for line in linelist:
            linetemp = line.split() # split in constituent tokens
            if len(linetemp) > 0:
                if currently_reading == "":   # seeking which slot to fill
                    if linetemp[0] == "#":
                        if linetemp[1] == "degree":  
                            num_slots += 1           # add 1 to slot count
                            currently_reading = "degree"
                            writefile.write("&degree\n")
                            continue
                        elif linetemp[1] == "shape":
                            num_slots += 1           # add 1 to slot count
                            currently_reading = "shape"
                            writefile.write("&shape\n")
                            continue
                        elif linetemp[1] == "rational":
                            num_slots += 1           # add 1 to slot count
                            currently_reading = "rational"
                            writefile.write("&rational\n")
                            continue
                        elif linetemp[1] == "knots":
                            if knots1_is_written == False:
                                num_slots += 1           # add 1 to slot count
                                currently_reading = "knots"
                                writefile.write("&knots_1\n")
                                continue
                            elif knots1_is_written == True:
                                num_slots += 1           # add 1 to slot count
                                currently_reading = "knots"
                                writefile.write("&knots_2\n")
                                continue
                        elif linetemp[1] == "points":
                            num_slots += 1           # add 1 to slot count
                            currently_reading = "points"
                            writefile.write("&control_points\n")             
                            continue
                        elif linetemp[1] == "weights":
                            num_slots += 1           # add 1 to slot count
                            currently_reading = "weights"
                            # this is never seen...
                            continue
                        else:
                            print('It seems there is an input file problem: ')
                            print( 'untagged data present?')
                elif currently_reading == "degree":
                    # first entry comes with a comma in original file
                    comma_position = linetemp[0].find(',')
                    tmp = linetemp[0][:comma_position]
                    writefile.write("    spline_deg1 = " + tmp +"\n")
                    writefile.write("    spline_deg2 = " + linetemp[1] +"\n")
                    writefile.write("/" + "\n\n")
                    currently_reading = ""
                    continue
                elif currently_reading == "shape":
                    comma_position = linetemp[0].find(',')
                    tmp = linetemp[0][:comma_position]
                    writefile.write("    num_pts1 = " + tmp +"\n")
                    writefile.write("    num_pts2 = " + linetemp[1] +"\n")
                    writefile.write("/" + "\n\n")
                    size_x1 = int(tmp)
                    size_x2 = int(linetemp[1])
                    currently_reading = ""
                    continue
                elif currently_reading == "rational":
                    writefile.write("    is_rational = "+ linetemp[0] + "\n")
                    writefile.write("/" + "\n\n")
                    currently_reading = ""
                    continue
                elif currently_reading == "knots":
                    if knots1_is_written == False:
                        logical_mesh_locations1 = remove_duplicates(linetemp)
                        tmp = " ".join(linetemp) # one space between elems
                        writefile.write("    knots1 = "+ tmp + "\n")
                        writefile.write("/" + "\n\n")
                        currently_reading = ""
                        knots1_is_written = True
                        continue
                    elif knots1_is_written == True:
                        logical_mesh_locations2 = remove_duplicates(linetemp)
                        tmp = " ".join(linetemp) # one space between elems
                        writefile.write("    knots2 = "+ tmp + "\n")
                        writefile.write("/" + "\n\n")
                        currently_reading = ""
                        continue
                elif currently_reading == "points":
                    if linetemp[0] == "#": # finished list of points
                        writefile.write("    control_pts1 = "+" ".join(x1)+"\n")
                        writefile.write("    control_pts2 = "+" ".join(x2)+"\n")
                        writefile.write("/" + "\n\n")
                        # and prepare for entering data in weights, since
                        # that's where we are at this point in the file...
                        writefile.write("&pt_weights\n")
                        currently_reading = "weights"
                        if point_counter != size_x1*size_x2 :
                            print("Warning: the number of points read is not ")
                            print("the same as the number expected.\n")
                        continue
                    else:  # still reading points
                        x1.append(linetemp[0])
                        x2.append(linetemp[1])
                        # we discard the third field which represents 'z'
                        point_counter += 1
                elif currently_reading == "weights":
                    wgts.append(linetemp[0])
                    # add a test to verify the right number of weights
                      #  if point_counter != size_x1*size_x2 :
                      #      print("Warning: the number of points read is not ")
                      #      print("the same as the number expected.\n")
                      #  continue
        writefile.write("    weights = "+" ".join(wgts)+"\n")
        writefile.write("/" + "\n\n")
        currently_reading = ""
        # add information relevant to the construction of the logical mesh.
        writefile.write("&logical_mesh_2d\n")
        nc1 = len(logical_mesh_locations1) - 1
        print('for the logical mesh implicit in the given transformation:')
        print('number of cells in dimension 1: ', nc1)
        writefile.write("    number_cells1 = " + str(nc1) + "\n")
        nc2 = len(logical_mesh_locations2) - 1
        print('number of cells in dimension 2: ', nc2)
        writefile.write("    number_cells2 = " + str(nc2) + "\n")
        writefile.write("/" + "\n\n")

def remove_duplicates( lst ):
    no_dups = []
    for i in lst:
        if i not in no_dups:
            no_dups.append(i)
    return no_dups


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='1.0')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        (options, args) = parser.parse_args()
        #if len(args) < 1:
        #    parser.error ('missing argument')
        if options.verbose: print( time.asctime())
        main()
        if options.verbose: print(* time.asctime())
        if options.verbose: print( 'execution time in seconds:')
        if options.verbose: print( (time.time() - start_time))
        sys.exit(0)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print( 'ERROR, UNEXPECTED EXCEPTION')
        print( str(e))
        traceback.print_exc()
        os._exit(1)


