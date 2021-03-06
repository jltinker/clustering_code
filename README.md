# clustering_code

command-line usage:

required libraries: in the initial commit I forgot that I call a library "cutil" in linking. This is a personal library, and it's been checked in as "lib" in my github. You just need to edit the path to your local libraries in the Makefile.

wp_LS_weight_omp galdat1 rand1.dat RR_file [covarfile] [collision_weight1] [collision_weight2] [rmin] [rmax] [nrbin] [njack_per_side] [pi_max]> wp.dat

galdat1 - galaxy input file. Format: Ascii
  1) ra [deg]
  2) dec [deg]
  3) redshift
  4) total weight of this target
  5) isurvey (just in case target sample combines two surveys, like BOSS+eBOSS)
  
rand1.dat - randoms. Format: Ascii
  1) ra [deg]
  2) dec [deg]
  
The redshifts of the randoms are drawn randomly from the data.

RR_file is a filename where the random counts are printed out. Since the same rand file can be used multiple times (and takes the majority of the CPU time) this can be read back in from this file. Currently the read-in is NOT IMPLEMENTED. If this file exists, it will be overwritten.

- covarfile - filename for outputting the covariance matrix
- collision_weight1 - if two galaxies are closer than 62 arcsecs, the pair is upweighted by this value (for survey==1)
- collision_weight2 - if two galaxies are closer than 62 arcsecs, the pair is upweighted by this value (for survey==2)
- rmin - minimum r_p value
- rmax - maximum r_r value
- nrbin - number of log(r_p) bins
- njack_per_side - for jackknife sampling: the same is divided into ra/dec partitions, for a total of (njack_per_side)^2 - - jackknife samples
- pi_max - value of pi out to which w_p is integrated.

Output--> standard out. Columns are:
  1)  mean w_p (pair-weighted over all pi)
  2)  w_p
  3) err w_p
  4) number of DD pairs
  5) mean w_p from all the jackknifes (for cross-check purposes)
  
---------------------------------------------------------------------------------

I have added code for doing the redshift-space multipoles. The command-line arguments are nearly identical:

ximulti_LS_weight galdat1 rand1.dat RR_file [covarfile] [collision_weight1] [collision_weight2] [rmin] [rmax] [nrbin] [njack_per_side]> ximulti.dat

Where all the inputs are the same, only now there is no pi_max (since no line-of-sight integration. The stdout of the code is:

 1) r [Mpc/h]
 2) xi_mono
 3) xi_mono_err
 4) xi_quad
 5) xi_quad_err

There are also two covariance matrix files that are outputted:

[covarfile]_mono and [covarfile]_quad


----------------------------------------------------------------------------------

I have now added cross-correlation for the wp code. Usage is:

wp_LS_weight galdat1 galdat2 rand1.dat RR_file [covarfile] [collision_weight1] [collision_weight2] [rmin] [rmax] [nrbin] [njack_per_side] [pi_max]> wp.dat

Here the usage is the same as the regular WP code, only now we have a galdata2 entry (expecting the same format as the first galdata1 entry). The output is the same as well. 

NB: The wp calculation is (D1D2-D1R-D2R-RR)/RR, which assumes that the two data samples have the same angular and redshift distribution.
