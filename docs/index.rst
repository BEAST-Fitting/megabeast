*********
MegaBEAST
*********

``megabeast`` is a Bayesian model that fits ensembles of BEAST results for
single stars to derive stellar population and dust/star geometry parameters.
Thus, the combination of the BEAST and megaBEAST is a Hierarchical Bayesian
model.  The goal is to create maps of stellar and dust parameters of galaxies
(e.g., A(V) or stellar age).

User Documentation
==================

.. toctree::
   :maxdepth: 2

   Motivation and Goals <megabeast/goals.rst>
   Settings file <megabeast/settings.rst>

   How to run <megabeast/running.rst>
   Creating simulated data <megabeast/simulate.rst>

Installation
============

.. toctree::
  :maxdepth: 2

  How to install <megabeast/install.rst>

Reporting Issues
================

If you have found a bug in ``megabeast`` please report it by creating a
new issue on the ``megabeast`` `GitHub issue tracker
<https://github.com/BEAST-Fitting/megabeast/issues>`_.

Please include an example that demonstrates the issue sufficiently so that
the developers can reproduce and fix the problem. You may also be asked to
provide information about your operating system and a full Python
stack trace.  The developers will walk you through obtaining a stack
trace if it is necessary.

Contributing
============

users.  We accept contributions at all levels, spanning the gamut from fixing a
typo in the documentation to developing a major new feature. We welcome
contributors who will abide by the `Python Software Foundation Code of Conduct
<https://www.python.org/psf/conduct/>`_.

``dust_extinction`` follows the same workflow and coding guidelines as
`Astropy`_.  The following pages will help you get started with contributing
fixes, code, or documentation (no git or GitHub experience necessary):

* `How to make a code contribution <https://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_

* `Coding Guidelines <https://docs.astropy.org/en/stable/development/codeguide.html>`_

* `Developer Documentation <https://docs.astropy.org/en/stable/#developer-documentation>`_


For the complete list of contributors please see the `megabeast
contributors page on Github
<https://github.com/BEAST-Fitting/megabeast/graphs/contributors>`_.

Reference API
=============
.. toctree::
   :maxdepth: 1

   megabeast/core_api.rst
   megabeast/tools_api.rst
