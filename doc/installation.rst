Installation
************

In the Sali Lab
===============

If you are working in the Sali lab, you don't need to build and install
PCSS - it is already set up for you as a module. Just run
``module load pcss`` to load it.

Dependencies
============

All dependencies listed below are expected to be found in standard
system paths. This may require setting ``PATH`` and/or
``LD_LIBRARY_PATH`` environment variables, or modifying the global parameter
file. Note that Linux is the only platform on which PCSS has been tested.

* `Perl <https://www.perl.org/>_`.

* `SVMlight <http://svmlight.joachims.org/>_`.

In the Sali lab, running 
``module load svm_light``
will get all of these dependencies.

Building
========

Use ``make install`` to install the library.
In most cases you will need to tell ``make`` where to install (if running on
a Linux cluster, PCSS will need to be installed on a network-accessible
filesystem), with something like
``make PREFIX=/shared/pcss install``. See
``Makefile.include`` for all make variables that can be configured.
