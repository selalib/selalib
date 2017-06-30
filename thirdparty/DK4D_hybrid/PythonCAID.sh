#!/usr/bin/env bash

source ~/bin/config_gnu_django.sh
module rm epd
module load python

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PIGASUS_DIR/lib/python2.7/site-packages/pigasus:$PIGASUS_DIR/lib
