#!/usr/bin/python
# test_pyCactus.py : tests functionality of pyCactus.py

from pyCactus import CactusRun
import matplotlib.pyplot as plt

run_directory = './data'
case_name = 'TestHAWT'

run = CactusRun(run_directory,case_name)

columns, headers = run.get_named_subset('time_Cp')

plt.plot(columns[0], columns[1])
plt.show()