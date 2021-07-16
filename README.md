# Rapid 3D z-crossings Algorithm (RAZA) FOR DISTRIBUTION			       

RAZA provides a robust, accurate, intuitive, fast, and generally applicable segmentation algorithm capable of detecting organelles, membranes, macromolecular as- semblies and extrinsic membrane protein domains. It defines each continuous contour within a tomogram as a discrete object and extracts a set of 3D structural fingerprints (major, middle and minor axes, surface area and volume), enabling selective, semi-automated segmentation and object extraction. 

For more details, please find the complete article here: https://pubmed.ncbi.nlm.nih.gov/29032142/


> ### Note: Please note that input file ‘GroEL.mrc’ and output files generated at each step are already provided with dataset. If you want to re-run this test data please use following commands to regenerate these results.

We assume you have successful installation of linux and RAZA source code on your computer, please run the following to run example dataset.

 Run following command that uses reference file (already provided with test dataset) and input MRC file to use particle selection feature of RAZA:

$ raza -s 1 -z -1 -r 50 -c 1 GroEL_16bit_scalled.mrc GroEL_16bit_scalled_s1_z-1_r50_c1.mrc


Please feel free to contact at the following address:

r.ali1@uq.edu.au
