.. Description of the prototype
.. module:: prototype
.. _prototype-page:

=================
Prototype
=================

The ``prototype`` is the first program which uses selalib


Documentation for Git is available `here <http://git-scm.com/>`_.

First time set up
-----------------

Set up your developer info::

 git config --global user.name "Your NAME"
 git config --global user.email "you@somewhere"
 git config --global core.editor emacs (or other) 
 git config --global merge.tool opendiff (or other)
 
Check configuration with::

 git config --list

For help::

 git help <command-you-want-help-about>

Get the prototype
-----------------
Developer GIT Access via SSH

`Edit keys on your INRIA gforge account <https://gforge.inria.fr/account/editsshkeys.php>`_ and follow instructions.
Only project developers can access the GIT tree via this method. Substitute with the proper values::

 git clone git+ssh://YOUR_INRIA_LOGIN@scm.gforge.inria.fr//gitroot//selalib/selalib.git
 cd selalib/

Display all branches with::

 git branch -a

Choose the remote prototype branch:: 

 git branch prototype-devel
 git checkout prototype-devel
 git merge origin/prototype-devel

More information available on document `An overview of GIT A short tutorial for SELALIB developers <https://gforge.inria.fr/docman/view.php/3042/7642/selalib_coding_guidelines.pdf>`_ by Edwin Chacon-Golcher.

Build the prototype
-------------------

Go to prototype/src/ directory, build the libraries and the program with ::

 cd prototype/src
 export HDF5INCPATH=/usr/include
 export HDF5LIBPATH=/usr/lib 
 export HDF5MODPATH=/usr/include 
 export FFTW3INCLUDE=/usr/include 
 make

You can set these variables in your shell profile (csh example) ::
 
 setenv HDF5INCPATH "/opt/local/include"
 setenv HDF5MODPATH "/opt/local/include"
 setenv HDF5LIBPATH "/opt/local/lib"
 setenv FFTW3INCLUDE "/opt/local/lib"

If library location is unusual, you probably have to set your runtime path with::

 setenv LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5LIBPATH}

Commit your modifications in a new branch
-----------------------------------------

Choose a good name for your new-branch::

 git status
 git add some_file.F90
 git status
 git commit -m 'New some_file that does something'
 git push origin prototype-devel:new-branch
