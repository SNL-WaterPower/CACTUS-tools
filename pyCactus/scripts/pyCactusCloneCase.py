#!/usr/bin/env python
"""Copy a CACTUS case directory without output data.

Copy all files except those matching:
    '*.csv',
    '*.pvd',
    '*.vts',
    '*.vtu',
    '*.out',
    '*.tp',
    'output'
"""

import os
import sys
import shutil
import fnmatch
import argparse

def ignore_list(path, files):
    # specify file extensions to ignore
    ignore_patterns = ['*.csv',
                       '*.pvd',
                       '*.vts',
                       '*.vtu',
                       '*.out',
                       '*.tp',
                       'output']
    
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


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".

    https://code.activestate.com/recipes/577058/
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description="""
Copy a CACTUS case directory without output data.

Copy all files except those matching:
    '*.csv',
    '*.pvd',
    '*.vts',
    '*.vtu',
    '*.out',
    '*.tp',
    'output'
""",formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("case_path",
                        help="path to original case.",
                        type=str)
    parser.add_argument("output_path",
                        help="path to copied case directory.",
                        type=str)
    parser.add_argument("--clean",
                        help="delete existing directory if it exists -- danger!",
                        action="store_true",
                        default=False)
    
    args = parser.parse_args()

    # read in command-line arguments
    case_path   = args.case_path
    output_path = args.output_path
    clean = args.clean

    # make directories
    if os.path.exists(output_path):
        if clean:
            query_yes_no("Are you sure you wish to delete %s?" % (output_path))
            try:
                shutil.rmtree(output_path)
            except RuntimeError:
                print "There was a problem deleting %s" % (output_path)
                raise
        else:
            raise RuntimeError('The path %s already exists!' % (output_path))

    # call the copy function
    clone_case_inputfiles(case_path,output_path)

if __name__ == "__main__":
    main()
