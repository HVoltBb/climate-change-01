/* Author/maintainer: Can Zhou [eidotog@gmail.com]
 * Date: Dec 1 2021
 * Version: 0.6
*/

#include <cmath>

namespace grow
{

	template <class T>
	T g<T>::growth_fn(T lt, T lx, T k, T delta_t, int type, T param)
	{
		T ans;
		switch (type)
		{
		case 0:
			ans = von_bert(lt, lx, k, delta_t);
			break;
		case 1:
			ans = west_g(lt, lx, k, delta_t);
			break;
		case 2:
			ans = gomprz(lt, lx, k, delta_t);
			break;
		case 3:
			ans = logistic(lt, lx, k, delta_t);
			break;
		case 4:
			ans = glogistic(lt, lx, k, param, delta_t);
			break;
		default:
			ans = g_von_bert(lt, lx, k, param, delta_t);
		};

		return ans;
	}

	template <class T>
	inline T g<T>::g_von_bert(T lt, T lx, T k, T m, T delta_t)
	{
		T ans = pow(lx - (lx - pow(lt, 1 - m)) * exp(-k * (1 - m) * delta_t), 1 / (1 - m));
		return ans;
	}

	template <class T>
	inline T g<T>::von_bert(T lt, T lx, T k, T delta_t)
	{
#ifdef _ZOOM
		return g_von_bert(lt, lx, k, .0, delta_t);
#else
		return lx - (lx - lt) * exp(-k * delta_t);
#endif
	}

	template <class T>
	inline T g<T>::west_g(T lt, T lx, T k, T delta_t)
	{
		return g_von_bert(lt, lx, k, 1 / 4.0, delta_t);
	}

	template <class T>
	inline T g<T>::gomprz(T lt, T lx, T k, T delta_t)
	{
		T ans = lx * exp(-log(lx / lt) * exp(-k * delta_t));
		return ans;
	}

	template <class T>
	inline T g<T>::logistic(T lt, T lx, T k, T delta_t)
	{
		T ans = lx / (1 + (lx - lt) / lt * exp(-k * delta_t));
		return ans;
	}

	template <class T>
	inline T g<T>::glogistic(T lt, T lx, T k, T v, T delta_t)
	{
		T ans = lx / pow(1 + (pow(lx / lt, v) - 1) * exp(-k * v * delta_t), 1 / v);
		return ans;
	}

	// Deprecated since version 5
	template <class T>
	T g<T>::effect_tp(T lt, T effect_size, int type)
	{
		T ans;
		switch (type)
		{
		case 0:
			ans = lt + effect_size;
			break;
		case 1:
			ans = lt * exp(effect_size);
		}
		return ans;
	}
};
