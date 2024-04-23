###### DATA_AVAILABILITY ######

    This repository contains 2 file:
    
 
    (1) One ".cpp" file for the estimate of the rescue/extinction probability:
	
			C++ source code "rescue-NK.cpp"
	
			Requires the GSL scientific library:
		
			Ubuntu systems installation through terminal:
			
	        (C++)	~ sudo apt-get install build-essential
			(GSL)	~ sudo apt-get install gsl-bin
			(GSL)	~ sudo apt-get install libgsl-dev
			
			Ubuntu systems compilation through terminal:
		
			~ c++ -O3 rescue-NK.cpp -o [executable_name] -lm -lgsl -lgslcblas
	
			After compilation, we will get an executable.

                       To run the code, we must provide it with input data.
                        
                       The input data are:
                          
                       i- sequence size
                       ii- epistasis parameter K
                       iii - initial population size N_0
		       iv - mutation probability
		       v - number of independent runs (n_landscapes)
	               vi - the number of runs for each fitness landscape (conf_land). The number of independent fitness landscape is v/vi
		       vii - strength of selection, which we make equal to 1
                       viii - stress level, delta
                       ix - Maximum fitness W_max
                       x - range of phenotypic values, b 

                       Here is an example of how to run the code in an Ubuntu terminal

                      ./[executable_name] 12 11 10000 0.001 100000 100 1.0 0.3 1.5 6.0

                      The output of the code is created in an .DAT file. The name of the file is RESCUE-NK and also includes the values of parameter used.

                      The output contains 2 columns: the value of the initial drop in fitness (parameter), the probability of extinction

                      The output was made simpler than in our version, which includes several other measurements, some not used.

    
    (2) One ".py" file for the theoretical estimate of the rescue/extinction probability:

     			Python source code "rescue_NK-theory.py"

 			This simpler version generates the probability of extinction versus the stress level, delta.

    			Ubuntu systems compilation through terminal:
		
			~ python3 rescue_NK-theory.py

       			The code contains the needed libraries, and the input parameters can be set inside the script. They are:

     			i - the range of delta values
			ii - maximum reproductive value Wmax
   			iii - carry capacity K
      			iv - initial population size N0
	 		v - mutation rate mut
    			vi - genome size L
       			vii - range of phenotypic values b

   			The output of the code is created in an .txt file. The name of the file is "theory_NK_gen" and also includes the values of parameter used.

