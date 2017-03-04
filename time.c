#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
  return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

  /* note: on AIX and presumably many other 32bit systems, 
   * clock() has only a resolution of 10ms=0.01sec 
   */ 
}


/* returns the time difference between two measurements 
 * obtained with second(). The routine takes care of the 
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0,double t1)
{
  double dt;
  
  dt=t1-t0;

  if(dt<0)  /* overflow has occured */
    {
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
    }

  return dt;
}
