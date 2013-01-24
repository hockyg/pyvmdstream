pyvmdstream
===========

A python library for interacting with VMD via TCP

Info
-------
Started January, 2013 by Glen Hocky

Based on the C++ vmdstream library by Aaron S. Keys - [https://github.com/askeys/vmdstream](https://github.com/askeys/vmdstream)

Authors
-------
Glen Hocky

Dependencies
------------
This library requires python and VMD to be previously installed. It has been tested to work with python2.7 and VMD 1.9.1. It also requires the python libraries numpy and optionally matplotlib.

Quick Instructions
------------------
### Before starting ###

 * VMD should be in PATH, meaning that if you type vmd on the command line, it will open.
 * python libraries, such as numpy, should be in PYTHONPATH

### To test ###

 * go to the directory with pyvmdstream.py, and type "python examples/draw_example.py"

### To use ###

 1. Include pyvmdstream.py in the same directory as your python script, or
 2. Add the directory where pyvmdstream.py is located to PYTHONPATH

Further instructions on usage may be found in:

**examples/draw_example.py**

### To get help ###

 * Type pydoc PATH-TO-LIBRARY/pyvmdstream.py

Contents
--------
* README.md

    This file, describing the overall contents of this project
* pyvmdstream.py

    The python library itself
* examples/

    Example usage files, which also serve as usage *instructions*

