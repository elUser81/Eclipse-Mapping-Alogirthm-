# Eclipse Mapping Alogirthm 

This script is designed to take flux vs time data for a Cataclysmic Variable
system and out put a map of the accretion disk around the primary. 
The code closely follows Baptista, Steiner (1993). See my thesis for more details
and references. 

To use this code, git clone this repository to your local machine, then create 
another python file and import eclipse_mapping_algorithm.py

Before you run the script, create a folder containing a parameter file and lightcurve data.

The basic call for the eclipse map script looks like:

    from eclipse_mapping.eclipse_mapping_algorithm import Processes, Binary

    binary_data_dir = 'path/to/folder/containing/flux data/and/parameter file'

    b = Processes.load_binary(binary_data_dir)
    Processes.build_compute_save(binary_data_dir)
    

Depending on your computer specs, the solver will take 4-6 hours for a 30 by 30 resolution.
I would recommend running the script at night before bed and checking in the morning!

---Note that the light curve data in Examples/Example Inputs is normalized over 
the maximum flux value in the data.
