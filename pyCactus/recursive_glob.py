# recursive_glob.py

import os
import fnmatch

def recursive_glob(rootdir='.', pattern='*'):
	""" A function to search recursively for files matching a specified pattern.
		Adapted from http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python """

	matches = []
	for root, dirnames, filenames in os.walk(rootdir):
	  for filename in fnmatch.filter(filenames, pattern):
		  matches.append(os.path.join(root, filename))

	return matches