/* Author/maintainer: Can Zhou [eidotog@gmail.com]
 * Date: Dec 1 2021
 * Version: 0.6
*/

#ifndef GROWTH_FUNCTIONS_XX
#define GROWTH_FUNCTIONS_XX

namespace grow
{
	template <class T>
	class g
	{
	public:
		// Growth function interface
		T growth_fn(T lt, T lx, T k, T delta_t, int type, T param);

		// (0) additive or (1) multiplicative effect of the environment on body length
		// for both cases, effect_size takes the full range of R
		// returns the expected length
		// deprecated since version 5
		T effect_tp(T lt, T effect_size, int type);

	private:
		// von Bertalanffy growth equation (Beverton, 1954)
		// Lx is your usual L infinity
		inline T von_bert(T lt, T lx, T k, T delta_t);

		// a growth model due to Geoffrey West et al (2001)
		// the original equation is on body weight
		// here a cubic relation between body length and body weight is assumed
		// lx is just a parameter, and not L infinity
		inline T west_g(T lt, T lx, T k, T delta_t);

		// General von Bertalanffy equation (1957 eqn 6)
		inline T g_von_bert(T lt, T lx, T k, T m, T delta_t);

		// Gompertz curve (Laird, 1964)
		inline T gomprz(T lt, T lx, T k, T delta_t);

		// Logistic curve (Verhulst, 1838)
		inline T logistic(T lt, T lx, T k, T delta_t);

		// Generalized logistic curve (Richards, 1959)
		inline T glogistic(T lt, T lx, T k, T v, T delta_t);
	};

}

#include "growth_imp.h"
#endif // GROWTH_FUNCTIONS_XX
