# clustering_code

command-line usage:

wp_LS_weight galdat1 rand1.dat RR_file [covarfile] [collision_weight1] [collision_weight2] [rmin] [rmax] [nrbin] [njack_per_side] [pi_max]> wp.dat

galdat1 - galaxy input file. Format:
  1) ra [deg]
  2) dec [deg]
  3) redshift
  4) total weight of this target
  5) isurvey (just in case target sample combines two surveys, like BOSS+eBOSS)
  
rand1.dat - randoms. Format:
  1) ra [deg]
  2) dec [deg]
  
The redshifts of the randoms are drawn randomly from the data.

RR_file is a filename where the random counts are printed out. Since the same rand file can be used multiple times (and takes the majority of the CPU time) this can be read back in from this file. Currently the read-in is NOT IMPLEMENTED. If this file exists, it will be overwritten.

covarfile - filename for outputting the covariance matrix
collision_weight1 - if two galaxies are closer than 62 arcsecs, the pair is upweighted by this value (for survey==1)
collision_weight2 - if two galaxies are closer than 62 arcsecs, the pair is upweighted by this value (for survey==2)
rmin - minimum r_p value
rmax - maximum r_r value
nrbin - number of log(r_p) bins
njack_per_side - for jackknife sampling: the same is divided into ra/dec partitions, for a total of (njack_per_side)^2 jackknife samples
pi_max - value of pi out to which w_p is integrated.
