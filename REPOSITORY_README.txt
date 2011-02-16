Some notes on the organization of the Selalib repository:

February 16, 2011: As othe initial commit, we intend the repo to have a
structure like this:
                        ________________
                       |   Selalib      |
                       |________________|
           ____________________|__________________
    ______|________     _______|_________   ______|________
   |   prototype   |   | gysela_public   | |    selalib    |
   |_______________|   |_________________| |_______________|

The prototype directory will contain the files related with the preliminary
testing and development of an architecture for Selalib. During the initial
stages, we expect most development to occur here. gysela_public 
contains the baseline public version of gysela, whose capabilities we intend
to keep compatible with our other developments. Selalib will include the more
mature set of capabilities that can be regarded as a library in itself. 
Initially it will only include some low-level utilities and will include 
higher-level functionality as time progresses.

Downloading the repository implies downloading all the project: prototype,
gysela_public and selalib. 
