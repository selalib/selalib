:orphan:

Selalib's releases 
==================


Selalib is available as an archive 

* `tar gz file <selalib-0.5.0.tar.gz>`_

#* as a `deb package <selalib_1.0.0_amd64.deb>`_
#* as a `rpm package <selalib-1.0.0.fc18.x86_64.rpm>`_
#
#Install the package on Ubuntu 12.10 with::
#
#    sudo dpkg -i selalib_1.0.0_amd64.deb
#
#or on fedora core 18 ::
#
#    sudo rpm -i selalib-1.0.0.fc18.x86_64.rpm

Use the library::

    gfortran -I/usr/lib/selalib/include landau.F90 -L/usr/lib/selalib/lib -lselalib -ldfftpack

    ./a.out | gnuplot


.. note:: A detailed installation and setup guide is in preparation.
