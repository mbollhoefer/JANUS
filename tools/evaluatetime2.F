#include <ilupack_fortran.h>
#include <ilupackmacros.h>

	FUNCTION EVALUATE_TIME(tu,ts) 
	IMPLICIT none 

	INTEGER, save :: cpt_init,cpt_max,freq
	REAL :: t1, tu, ts, sec, EVALUATE_TIME, EVALUATETIME
	REAL(kind(0d0)) :: cputm,eltm
	EXTERNAL EVALUATETIME

#if defined __WINDOWS__
		CALL CPU_TIME(t1)	
		sec=t1;	 
#else
		sec=EVALUATETIME(tu,ts)
#endif

	EVALUATE_TIME=sec
	RETURN 
	END FUNCTION
