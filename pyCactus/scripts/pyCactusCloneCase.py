#!/usr/bin/env python
"""Copy a CACTUS case directory without output data."""

import os
import sys
import shutil
import fnmatch
import argparse


def find_files(directory, pattern='*'):
    """Find files recursively in a directory matching a glob-style pattern.

    e.g., find_files('.', 'output/element/*.csv')
    """
    if not os.path.exists(directory):
        raise ValueError("Directory not found {}".format(directory))

    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            full_path = os.path.join(root, filename)
            if fnmatch.filter([full_path], pattern):
                matches.append(os.path.join(root, filename))
    return matches


def clone_case(case_path, output_path, ignore_patterns):
    """Clone the case, ignoring glob-style patterns."""
    # find all the ignored files
    ignored_files = []
    ignored_files_dict = {}
    for ignore_pattern in ignore_patterns:
        found_files = find_files(case_path,
                                 pattern=os.path.join(case_path,
                                                      ignore_pattern))
        ignored_files += found_files
        ignored_files_dict[ignore_pattern] = len(found_files)

    def __ignore_this(path, files):
        matches = []
        for file in files:
            if os.path.join(path, file) in ignored_files:
                matches = matches + [file]
        return matches

    # copy all files except those matching the pattern
    shutil.copytree(case_path, output_path,
                    ignore=__ignore_this)

    # remove empty directories (topdown)
    # http://stackoverflow.com/questions/18165343/python-empty-dirs-subdirs-after-a-shutil-copytree-function
    for root, dirs, files in os.walk(output_path, topdown=True):
        for dn in dirs:
            pth = os.path.join(root, dn)
            try:
                os.rmdir(pth)
            except OSError:
                pass

    print 'Copied%s to %s' % (case_path, output_path)

    if len(ignored_files) > 0:
        print 'Files ignored:'
        for ignore_pattern, num_found_files in ignored_files_dict.iteritems():
            print ignore_pattern + '\t' + str(num_found_files)


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
    """Main."""
    # parse command line arguments
    parser = argparse.ArgumentParser(description="""
Copy a CACTUS case directory without output data.

By default, no output data is copied (file matching the following patterns are
ignored).

    '*.csv',
    '*.pvd',
    '*.vts',
    '*.vtu',
    '*.out',
    '*.tp',
    'output'

Flags may be used to copy output data (see help).""", formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("case_path",
                        help="path to original case.",
                        type=str)
    parser.add_argument("output_path",
                        help="path to copied case directory.",
                        type=str)
    parser.add_argument("--clean",
                        help="delete existing directory if it exists!",
                        action="store_true",
                        default=False)
    parser.add_argument("--copy-output",
                        help="copy output data.",
                        action="store_true",
                        default=False)
    parser.add_argument("--no-field",
                        help="ignore field data (only used if --copy-output).",
                        action="store_true",
                        default=False)
    parser.add_argument("--no-wakeelem",
                        help="ignore wake element data (only used if --copy-output).",
                        action="store_true",
                        default=False)
    parser.add_argument("--no-wall",
                        help="ignore field data (only used if --copy-output).",
                        action="store_true",
                        default=False)

    args = parser.parse_args()

    # read in command-line arguments
    case_path   = os.path.abspath(args.case_path)
    output_path = os.path.abspath(args.output_path)
    copy_output = args.copy_output
    no_field  = args.no_field
    no_wakelem = args.no_wakeelem
    no_wall = args.no_wall
    clean = args.clean

    # define patterns
    field_data_pattern = ['output/field/*']
    wakeelem_data_pattern = ['output/element/*']
    wall_data_pattern = ['output/wall*/']

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

    # set the ignore patterns
    if not copy_output:
        ignore_patterns = ['*.csv',  # csv files
                           '*.vtu',  # vtk unstructured
                           '*.vts',  # vtk structured
                           '*.pvd',  # vtk collection
                           '*.tp',   # tecplot file
                           '*.out',  # slurm output
                           'output/*']  # anything in the 'output' dir

    elif copy_output:
        ignore_patterns = []
        if no_field:
            ignore_patterns += field_data_pattern

        if no_wakelem:
            ignore_patterns += wakeelem_data_pattern

        if no_wall:
            ignore_patterns += wall_data_pattern

    # call the copy function
    clone_case(case_path, output_path, ignore_patterns)

if __name__ == "__main__":
    main()
