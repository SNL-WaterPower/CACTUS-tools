#!/usr/bin/env python
# pyCactusCloneCase.py
# =====================
# Copy a CACTUS case directory WITHOUT output data.

import os
import sys
import shutil
import fnmatch
import argparse

def ignore_list(path, files):
    # specify file extensions to ignore
    ignore_patterns = ['*.csv', '*.pvd', '*.vts', '*.vtu', '*.out', '*.tp']
    
    ignored_files = []
    
    for ignore_pattern in ignore_patterns:
        ignored_files = ignored_files + fnmatch.filter(files,ignore_pattern)

    return ignored_files


def clone_case_inputfiles(case_path,output_path):
    # copy the files (using the ignore list function above)
    shutil.copytree(case_path, output_path, ignore=ignore_list)

    # remove empty directories (topdown)
    # http://stackoverflow.com/questions/18165343/python-empty-dirs-subdirs-after-a-shutil-copytree-function
    for root, dirs, files in os.walk(output_path, topdown=True):
        for dn in dirs:
            pth = os.path.join(root, dn)
            try:
                os.rmdir(pth)
            except OSError:
                pass

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(description="""Copy a CACTUS case directory WITHOUT output data.'""")

    parser.add_argument("case_path", help="path to original case.", type=str)
    parser.add_argument("output_path", help="path to copied case directory.", type=str)
    
    
    args = parser.parse_args()

    # read in command-line arguments
    case_path   = args.case_path
    output_path = args.output_path

    # make directories
    if os.path.exists(output_path):
        print 'Error: The path already exists!' + output_path
        sys.exit()

    # call the copy function
    clone_case_inputfiles(case_path,output_path)