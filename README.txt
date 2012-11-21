
This Pipol package expose as a simple example the automation of the
chain :

     1. architecture reservation, OS deployement, OS configuration
     2. svn update of the sources from gforge 
     3. compilation
     4. test and visualisation on a dashboard
     5. packages generation
     6. packages upload on the gforge

For details on Pipol automations, see :
http://pipol.inrialpes.fr/docs/user/user-manual.html

Under YourPipolHomeDir/.pipol (see user manual):

 - some examples of operating system configuration scripts rc and
   rc.<system>

 - under nightly/ a nightly build script example with comments


As this example is not a gforge project, svn update and packages
upload are commented in

YourPipolHomeDir/.pipol/nightly/HelloWorldBuild


Quick Start:
-----------

I) manual Hello World compilation, testing and packaging

mkdir $PIPOL_WDIR/build
cd $PIPOL_WDIR/build
cmake <the path of this directory>
make Experimental
(the test result goes to http://cdash.inria.fr/CDash/index.php?project=Pipol)

#package generation (with cmake >= 2.6)
cpack -G <the cpack generator for your system : DEB, RPM, ... >


II) batch build with a script : 

in YourPipolHomeDir/.pipol/nightly/HelloWorldBuild we assume this
package as been unpacked in your pipol home dir under $HOME/src

If it is not the case modify the script
YourPipolHomeDir/.pipol/nightly/HelloWorldBuild

Then you can have a test :

pipol-sub esn amd64-linux-debian-testing.dd none 2:00 $PWD/YourPipolHomeDir/.pipol/nightly/HelloWorldBuild 

See the result on :

http://cdash.inria.fr/CDash/index.php?project=Pipol

# example of system configuration :

1) get cmake-2.6.0.tar.gz in $HOME/.pipol/src 
cd $HOME/src
wget http://www.cmake.org/files/v2.6/cmake-2.6.0.tar.gz
cp PipolHelloWorld-1.0.0-Source/YourPipolHomeDir/.pipol/rc $HOME/.pipol/rc

The previous pipol-sub command should correctly create the packages

III) nightly build

As this example would lead to many reservations at night, this should
be setup only for usefull build. So please at the beginning use only a
limited set of #PIPOL lines and extend them as needed.

look at PipolHelloWorld-1.2-Source/YourPipolHomeDir/.pipol/nightly

and create your own $HOME/.pipol/nightly/ProjectBuild
