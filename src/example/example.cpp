/* This is a stripped down version of v5.cpp. With most of the control
 * statements gone, this example should be helpfull for readers not familiar
 * with TMB to understand the core structure of the code.
 *
 * Note that the same effect of this example code can be achieved with v5.cpp,
 * and the inclusion of this file is for learning purposes only.
 *
 * Author/maintainer: Can Zhou [eidotog@gmail.com]
 * Date: May 29 2021
 * Version: 0.1
*/

#include <TMB.hpp>
#include "../growth.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
	using namespace density;
	using namespace Eigen;

	// Data
	// Length measurements
	DATA_VECTOR(l1);
	DATA_VECTOR(l2);

	// Time and date
	DATA_IVECTOR(t1);	//Month #
	DATA_VECTOR(d1);	//Day #/ # of days in that month
	DATA_IVECTOR(t2);
	DATA_VECTOR(d2);
	DATA_VECTOR(days);

	DATA_INTEGER(PAR1);

	
	PARAMETER(c);
	PARAMETER(k);
	PARAMETER(lognu);
	PARAMETER(logSig);

	Type nll = (Type) .0;

	// Local variables
	Type lt;
	int ds;
	int td;
	int tb;
	
	// Likelihood	
	for(int i=0; i < l1.size(); i++){
		lt = l1(i);
		tb = t1(i);
		td = t2(i) - t1(i);
		if(td == 0){		
			// Released and recaptured in the same month
			lt = grow::fish<Type>.growth_fn(lt, c, k, d2(i) - d1(i), PAR1, exp(lognu));
		} else {			
			// Released and recaptured in different months
			// First partial month
			lt = grow::fish<Type>.growth_fn(lt, c, k, 1 - d1(i), PAR1, exp(lognu));
			// Any full calendar month(s) in between
			for(int j = 1; j < td; j++){
				lt = grow::fish<Type>.growth_fn(lt, c, k, days(tb+j), PAR1, exp(lognu));
			}
			// Second partial month
			lt = grow::fish<Type>.growth_fn(lt, c, k, d2(i), PAR1, exp(lognu));
		}
		nll -= dnorm(l2(i), lt, exp(logSig), true);
	}

	return nll;
}
