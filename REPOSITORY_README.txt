Some notes on the organization of the Selalib repository:

As of the initial commit, we intend the repo to have a file
structure like this:
                      ________________
                     |    Selalib     |
                     |________________|
           __________________|___________________________________________
    ______|________   _______|_________   ______|________   ___|___   __|____
   |   prototype   | |     selalib     | | gysela_public | | slv2d | | pic2d |
   |_______________| |_________________| |_______________| |_______| |_______|

The 'prototype' directory will contain the files related with the preliminary
testing and development of an architecture for Selalib. During the initial
stages, we expect most development to occur here. 

'selalib' will include the more mature set of capabilities that can be 
regarded as a library in itself. Initially it will only include some low-level 
utilities and will include higher-level functionality as time progresses.

'gysela_public' will contain the baseline public version of Gysela, whose 
capabilities we intend to keep compatible with our other developments. However,
since the public version of Gysela is derived from the proprietary version
owned by the CEA, only the CEA personnel will work updating this part of the
repository.

'slv2d' is a code written by Eric Sonnendrucker that contains all functions
we need to have in selalib. This part of the code will be used for testing purpose

'pic2d' is a code written by Regine Bartheleme that contains a test case using PIC 
method. This part of the code will be used for testing purpose of PIC features of selalib.

AVAILABLE BRANCHES:

Our Selalib repository has to hold information that while related, is not 
really unified. For example, the public version of Gysela, the prototype, and 
the modules of the Selalib library. Maybe even more. At the same time, we have 
only one repository to store all of this information. It would be convenient 
if we would use Git to apply some level of separation among these entities.

To manage the work for the prototype and selalib components, we have chosen to 
have for each of them -main and -devel branches. A -main branch should contain 
the most trusted, least frequently updated set of code and documentation files. In a more advanced state of the project this would represent something like an 
official release. Merging into a -main branch should be the result of 
completing some integrity tests, such as whole system builds, test runs and 
when available, non-regression tests. Non-code files, such as documentation, 
should be subject to some preliminary review. Only a few individuals should be 
authorized to merge with the -main branches and do so under certain 
pre-defined guidelines.

The -devel branches should contain the files that have yet to be put through 
the integrated tests necessary for merging into the -main branches. To be in a 
-devel branch, code should aim at being accompanied by the corresponding unit 
test. For some really preliminary development that needs to be shared, coders 
are encouraged to make their own short-lived topic branches, that can be 
merged and then deleted when the work is done, much like the local workflow 
that one uses with Git.

The independent branch gysela-public, is meant to exist so that the CEA can 
upload the public version of Gysela directly into its own branch. The mergings 
of branches are supposed to flow always in the direction of: 

                       devel ---> main ---> master

and

                       gysela-public ---> master

Therefore, at any given moment, master should have the most comprehensive 
collection of trusted files for the whole project. Anything below that 
represents a sub-project and a given stage of development.

Thus, the available branches are:
master
gysela-public
prototype-main
prototype-devel
selalib-main
selalib-devel


Downloading the repository implies downloading all the project and all the 
branches. Specifics as to how to deal with the branching model in Git are
further explained in the 'An Overview of Git: A short tutorial for Selalib 
developers'. 
