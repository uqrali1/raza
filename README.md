# Rapid 3D z-crossings Algorithm (RAZA) FOR DISTRIBUTION			       

RAZA provides a robust, accurate, intuitive, fast, and generally applicable segmentation algorithm capable of detecting organelles, membranes, macromolecular as- semblies and extrinsic membrane protein domains. It defines each continuous contour within a tomogram as a discrete object and extracts a set of 3D structural fingerprints (major, middle and minor axes, surface area and volume), enabling selective, semi-automated segmentation and object extraction. 


If you use RAZA software in your research, please cite the following paper: 
Rubbiya A. Ali, Ahmed M. Mehdi, Rosalba Rothnagela, Nicholas A. Hamilton,Christoph Gerlec, Michael J. Landsberga, Ben Hankamer "RAZA: A Rapid 3D z-crossings Algorithm to segment electron tomograms and extract organelles and macromolecules" (2017), Journal of Structural Biology 200 (2):73-86
https://pubmed.ncbi.nlm.nih.gov/29032142/


> ### Note: Please note that input file ‘GroEL.mrc’ and output files generated at each step are already provided with dataset. If you want to re-run this test data please use following commands to regenerate these results.

We assume you have successful installation of linux and RAZA source code on your computer, please run the following to run example dataset.

 Run following command that uses reference file (already provided with test dataset) and input MRC file to use particle selection feature of RAZA:

$ raza -s 1 -z -1 -r 50 -c 1 GroEL_16bit_scalled.mrc GroEL_16bit_scalled_s1_z-1_r50_c1.mrc


The software package provides source code, installation instructions, a test dataset and tutorial. Please feel free to contact at the following address:

### r.ali1@uq.edu.au


## complete userguide

> Please note the software should be installed successfully to run following command.
#### 1.Open terminal and run following command:
> $ raza

#### 2.	How to make RAZA work on your dataset?
>Current version of RAZA is a command line tool, and does not require the parallelization or GPU based cluster. It is autonomous software, which can be further combined with other image visualization packages. RAZA works on 3D MRC files only. If you have any other format, please convert data to MRC format. Here, we are using our test dataset ‘GroEL.mrc’ to run RAZA. The details of this dataset have been provided in manuscript.

Following commands can be used in this software package;
##### I.	tomogramNormalization
input: MRC file
output: Normalized output file (MRC file format). The output file can be visualized in IMOD.
##### II.	raza_log
input: MRC file
output: Laplacian of Gaussian of input file (Non-binary MRC file format). Highlighted edges can be visualized using IMOD.
##### III.	raza_logZcross
input: MRC file
output: Laplacian of Gaussian of input file followed by z-crossings (binary MRC file format). The output file can be visualized in IMOD.
##### IV.	raza
input: Input volume (MRC file) and reference/template object details (text file)
output: Particle selection after performing raza_logZcross
(the output file has a binary MRC file format showing selected particles based on reference structural fingerprints). The output file can be visualized in IMOD.

Once RAZA has been compiled successfully, please do following:
a.	Make a directory ‘Test’, and copy input file ‘GroEL.mrc’
> $ cd Desktop
> $ mkdir Test
> $ cd Test/
> $ ls
> GroEL.mrc
	
> #### Steps (b-e) detailed below are used for optimization purpose only. You may not need to run steps (b-e). You can jump to step (f), if you wish to run RAZA directly on your dataset. 

b.	In the header of input file check whether the file is 8-bit, 16-bit or 32-bit format. If it is not 16-bit, then convert it to 16-bit by using ‘newstack’
> $ newstack -mode 1 GroEL.mrc GroEL_16bit.mrc
> $ header GroEL_16bit.mrc
 Number of columns, rows, sections .....     398     520      95
 Map mode ..............................    1(16-bit integer)           
 
c.	Data is standardized by normalizing the intensity range of the input volume to -32767 to 32767 using ‘tomgramNormalization’ program. Normalization of data is highly recommended. However, RAZA can be run without normalization procedure if you wish not to normalize your data.
> $ tomogramNormalization GroEL_16bit.mrc GroEL_16bit_scaled.mrc

An output file “GroEL_16bit_scaled.mrc” has been generated. You will use this normalized MRC file for further processing. 

d.	Now apply Laplacian of Gaussian (LoG). To find out usage details please enter: ‘raza_log’ and hit enter. 

> #### Recommended values for sigma (σ) are 0.5-3.0.
> $ raza_log 

raza_log by Rubbiya Akram Ali 
Usage: raza_log [options] <input file> <output file>
Options:
	-s #	Sigma (σ) of LoG filter (default 0.5) 

> #### The above step is not compulsory before running RAZA. It helps you to test different σ values on your dataset and visualize your output using IMOD.

e.	Now use the optimized σ values to run RAZA LoG filter. You may skip step (e), if you would like to run RAZA LoG filter directly without performing optimization step.
> $ raza_log -s 1 GroEL_16bit_scaled.mrc GroEL_16bit_scaled _s1.mrc

f.	Now apply z-crossings. To find out usage details please enter: ‘raza_logZcross’ and hit enter. 
> $ raza_logZcross 

raza_logZcross by Rubbiya Akram Ali 
Usage: raza_logZcross [options] <input file> <output file>
Options:
	-z #	Z Crossing value, threshold for gradients (default 1.00)
	-s #	Sigma (σ) of LoG filter (default 0.00) 

g.	Now analyze intensity values and select z-crossings value. To apply z-crossings run following command:
> $ raza_logZcross -s 1 -z -1 GroEL_16bit_ scaled.mrc GroEL_16bit_ scaled _s1_z-1.mrc

4.	Particle selection using RAZA
In order to use the particle selection feature of RAZA, Steps (a-g) can help user to optimize σ and z-crossings values.

> ### Selection of reference objects
Before using particle selection feature of RAZA, we first generate a text file containing details of number of contours, of reference object(s). 
Use following commands to select structural fingerprints of reference objects: The functions/methods used below for selecting objects as reference have not been constructed by us, rather they are part of IMOD package.
 
> $ 3dmod  GroEL_16bit_scalled.mrc

In IMOD, click on 
	3dmod Zap Window 
	Go to ‘Specials’
	Drawing tools 
	Select Drawing Mode: Sculpt
	Then make a rough circle on edges of your object
	Once you are done with one object, go to ‘Edit’ 
	Object -> New -> Done

then click on “new object” button (it’ll give you a new color). For each reference particle selection, you need to create new object with new color. Now press ‘S’ to save, and give ‘UserInputModel.mod’ as model name. Make sure model file is saved in the same parent directory which contains your input file. Then you need to convert your model file to x, y and z coordinates using following command.
> $ model2point -object -contour UserInputModel.mod UserInputModel_Model2Point.txt

The output text file UserInputModel_Model2Point.txt will be used as a reference structural fingerprints for particle selection. Once you have defined optimal values of both the z-crossings and σ parameters, check usage of ‘raza – particle picking tool’.
> $ raza

raza by Rubbiya Akram Ali 
Usage: raza [options] <input file> <output file>
Options:
	-z #	z-crossings value, threshold for gradients (default 1.00)
	-s #	Sigma (σ) of LoG filter (default 0.00) 
	-o #	Output only the given Z slice (numbered from 1)
	-r #	Tolerance rate (default 10)
	-c #	Condition for testing (1-7)
		  	1. All following condition together
		  	2. 3D Axes together
		  	3. 3D volume
		  	4. 3D surface area
		  	5. Major axis only
		  	6. Middle axis only
		  	7. Minor axis only

Run following command that uses reference file and input MRC file to use particle selection feature of RAZA:

> $ raza -s 1 -z -1 -r 50 -c 1 GroEL_16bit_scalled.mrc GroEL_16bit_scalled_s1_z-1_r50_c1.mrc

> It will generate four output files. 
> #### 1.	GroEL_16bit_scalled_s1_z-1_r50_c1.mrc
> #### 2.	model_file_major_minor.txt
> #### 3.	Final_subvolume_coordinates.txt
> #### 4.	StructuralDetails.txt
The output file is an MRC file that contains the selected particles based on reference volume. ‘model_file_major_minor.txt’ contains x, y and z coordinates of each major, middle and minor axes. ‘Final_subvolume_coordinates.txt’ contains the coordinates which can be used to extract objects using ‘trimvol’ (function) in IMOD package. ‘structuralDetails.txt’ contains structural details (length of major, middle and minor axes, surface area and volume) of each object. 
RAZA generates an MRC file as the output.  Post processing can therefore be conducted using MRC compatible image processing tools.
To link RAZA edge detection algorithms to post processing tools of the IMOD package. Typically, ‘imodauto’ is applied to convert the RAZA output to generate a stack of 2D contours that can be visualised in IMOD. The IMOD ‘Drawing tool’ is used for the removal of unwanted contours. ‘imodmesh’ is used to render the 3D surfaces. Alternatively, ‘isosurface’ can be used to generate a 3D surface directly, when specific contours do not have to be deleted.  


## Using singularity (>= 3.4)

To pull: 

	singularity pull --arch amd64 library://hoangnguyen177/default/raza:0.1

To run, cd to the input:

	singularity exec ${PATH_TO_RAZA_SIF} raza -s 1 -z -1 -r 50 -c 1 ./GroEL_16bit_scalled.mrc ./GroEL_16bit_scalled_s1_z-1_r50_c1.mrc

## Citation

Ali, R.A., Mehdi, A.M., Rothnagel, R., Hamilton, N.A., Gerle, C., Landsberg, M.J. and Hankamer, B., 2017. RAZA: a Rapid 3D z-crossings algorithm to segment electron tomograms and extract organelles and macromolecules. Journal of structural biology, 200(2), pp.73-86.
