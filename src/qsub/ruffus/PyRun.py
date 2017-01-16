#!/usr/bin/env python
#make executable in bash chmod +x PyRun

import sys
import inspect
import importlib
import os

if __name__ == "__main__":
    cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
  
    if cmd_folder not in sys.path:
        sys.path.insert(0, cmd_folder)

    # get the second argument from the command line      
    methodname = sys.argv[1]

    # split this into module, class and function name
    modulename, classname, funcname = methodname.split(".")
    
    # get pointers to the objects based on the string names
    themodule = importlib.import_module(modulename)
    theclass = getattr(themodule, classname)
    thefunc = getattr(theclass, funcname)

    # pass all the parameters from the third until the end of 
    # what the function needs & ignore the rest
    args = inspect.getargspec(thefunc)
    z = len(args[0]) + 2
    params=sys.argv[2:z]
    thefunc(*params)
