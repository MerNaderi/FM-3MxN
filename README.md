# FM-3MxN
Three-way data clustering based on the mean-mixture of matrix-variate normal distributions

by: Mehrdad Naderi · Mostafa Tamandi · Elham Mirfarah · Wan-Lun Wang · Tsung-I Lin


Please copy the files to the "current working directory" of the R package.

The getwd() function shall determine an absolute pathname of the "current working directory". 


./Functions

	contains 

        (1) MxN-mix.EM.r: EM-based estimating program for the FM-MxN distributions;

        (2) MxESN-mix.EM.r: EM-based estimating program for the FM-MxESN and FM-MxSN distributions;

        (3) MeMxN-mix.EM.r: EM-based estimating program for the FM-MeMxN distributions;

        (4) MgMxN-mix.EM.r: EM-based estimating program for the FM-MgMxN distributions;

        (5) MwMxN-mix.EM.r: EM-based estimating program for the FM-MwMxN distributions;

        (6) MsMxN-mix.EM.r: EM-based estimating program for the FM-MsMxN distributions;

        (7) Prior-initial.r: function for generating initial values;


Notes: 
    
      1. The "moments", "kmed", and "cli" R packages are required;

      2. The "Insurance" data set is avalible in the R package "splm";

