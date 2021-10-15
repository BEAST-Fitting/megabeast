#####################
Running the MegaBEAST
#####################

**********
Background
**********

The MegaBEAST is specifically setup to use the output of BEAST
runs.  Currently it can be run on a single set output from a single
BEAST run.

Plans and prototype code exist for running the MegaBEAST on
BEAST results ordered spatially to have
all the stars in a region split in to smaller pixels.
More information on the BEAST spatial reordering can be found in the
`BEAST documentation <http://beast.readthedocs.io/en/latest/workflow.html#post-processing>`_.

*******
Running
*******

The MegaBEAST should be installed (see :ref:`mb-install`.)

The MegaBEAST is run from the commandline using a :ref:`mb-settings`.

For fitting a single stellar population with a foreground screen of dust:

.. code-block:: shell

  $ mb_fit_single mb_settings.txt

*****************
MegaBEAST outputs
*****************

N/A
