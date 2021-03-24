# BELLA: a method to design blended multi-panel composite laminates with many plies and panels

--------------------------------------------------------------------------

The files correspond to the PhD Thesis of Noémie Fedon.

BELLA is a method for optimising the layout of composite laminate strucutres when panel 
thicknesses are fixed. The algorithm enforoces lay-up and ply-drop design guidelines and 
optimise panel convergence towards lamination-paramter targets.

--------------------------------------------------------------------------
Requirements:

1. A python IDE running python 3.7+

2. The following libraries should accessible:

	- matplotlib
	- pandas
	- numpy

---------------------------------------------------------------------------
In order to use it:

1. clone or download the repository in your desired folder.

2. Set up your environment with the appropriate libraries.

3. Change the settings and run one of the files used to test LAYLA: 
run_BELLA.py, run_BELLA_from_input_file.py, run_BELLA_from_input_file_horseshoe.py 

Refer to the documentation for more details.
--------------------------------------------------------------------------
Folder Structure

- src and subfolders contain the files needed to run the code

- input-files contains the files storing the input-files used for testing LAYLA

- results contains the results and analyses generated for the thesis.

- FXI-results contains the results generated using the evolutionary algorithm of 
François-Xavier Irisarri.

- run_BELLA.py is used for to run BELLA without using an input-file.

- run_BELLA_from_input_file.py is used for testing BELLA based on input-files.

- run_BELLA_from_input_file_horseshoe.py is used for testing BELLA based on input-files 
specific to the benchmark problem of composite-laminate design.

--------------------------------------------------------------------------
Version 1.0.0

--------------------------------------------------------------------------
License

This project is licensed under the MIT License. See the LICENSE for details.

--------------------------------------------------------------------------
Acknowledgments

- Terence Macquart, Paul Weaver and Alberto Pirrera, my PhD supervisors.

--------------------------------------------------------------------------
Author:

Noémie Fedon

For any questions or feedback don't hesitate to contact me at: noemiefedon@gmail.com
or through github at https://github.com/noemiefedon/BELLA
