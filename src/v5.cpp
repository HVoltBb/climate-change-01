/* Author/maintainer: Can Zhou [eidotog@gmail.com]
 * Date: Dec 1 2021
 * Version: 0.6
 */

#include <TMB.hpp>
#include "growth.h"

// v0.4 added measurement error
// v0.5 removed l
// v0.6 removed c++14 warning, added parallel option

template <class Type>
Type objective_function<Type>::operator()()
{
	using namespace density;
	using namespace Eigen;
	grow::g<Type> shark;

	// Data
	DATA_VECTOR(l1);
	DATA_IVECTOR(t1); // Month #
	DATA_VECTOR(d1);  // Day #/ # of days in that month
	DATA_IVECTOR(t2);
	DATA_VECTOR(d2);

	// Higher resolution
	DATA_IVECTOR(d1_d); // Day #
	DATA_IVECTOR(d2_d);
	DATA_VECTOR(nao_d);
	DATA_MATRIX(X_d);
	DATA_IVECTOR(MNS_map);

	DATA_VECTOR(l2);
	DATA_IVECTOR(SEX);
	DATA_IVECTOR(DUPE);
	DATA_IVECTOR(MNS);
	DATA_INTEGER(PAR1);
	DATA_INTEGER(PAR2);
	DATA_SCALAR(PAR3);
	DATA_INTEGER(PAR4);
	DATA_INTEGER(PAR5);
	DATA_VECTOR(nao);
	DATA_VECTOR(days);

	// Reference
	DATA_INTEGER(REFSEX);
	DATA_SCALAR(REFLEN);
	DATA_SCALAR(REFTIM);
	DATA_INTEGER(REFMON);

	// Prediction
	DATA_VECTOR(days_p);
	DATA_INTEGER(lag_);

	// Cubic spline
	DATA_MATRIX(X);
	DATA_SPARSE_MATRIX(S);
	DATA_INTEGER(Sdim);
	DATA_MATRIX(prediction_design_matrix);
	DATA_INTEGER(Pdim);
	DATA_MATRIX(p_dm);

	// NAO effects
	PARAMETER(c);
	PARAMETER(k);
	PARAMETER(b);
	PARAMETER_VECTOR(gc);
	PARAMETER(logNaoc);
	PARAMETER_VECTOR(gk);
	PARAMETER(logNaok);

	// Observation error
	PARAMETER(logSig);
	DATA_INTEGER(PAR6);
	PARAMETER_VECTOR(e_o);

	// Sex effects
	PARAMETER_VECTOR(lsex);
	PARAMETER(logSexl);
	PARAMETER_VECTOR(ksex);
	PARAMETER(logSexk);

	// Month effects
	PARAMETER_VECTOR(theta_l);
	PARAMETER(logMonthl);

	// Individuality
	PARAMETER_VECTOR(indl);
	PARAMETER(logvIndl);

	// Generalized logistic
	PARAMETER(lognu);

#ifdef _OPENMP
	parallel_accumulator<Type> nll(this);
#else
	Type nll = .0;
#endif
	// Guarded individuality
	if (CppAD::Variable(logvIndl))
		nll -= sum(dnorm(indl, Type(.0), exp(logvIndl), true));

	// Guarded Sex effect
	if (CppAD::Variable(logSexl))
		nll -= sum(dnorm(lsex, Type(.0), exp(logSexl), true));
	if (CppAD::Variable(logSexk))
		nll -= sum(dnorm(ksex, Type(.0), exp(logSexk), true));

	// Guarded month effect
	if (CppAD::Variable(logMonthl))
		nll -= sum(dnorm(theta_l, Type(.0), exp(logMonthl), true));

	vector<Type> naoc = X * gc;
	vector<Type> naoc_d = X_d * gc;
	if (CppAD::Variable(logNaoc))
		nll -= 0.5 * Sdim * logNaoc - 0.5 * exp(logNaoc) * GMRF(S).Quadform(gc);

	if (CppAD::Variable(b))
	{
		naoc = nao * b;
		naoc_d = nao_d * b;
	}

	vector<Type> naok = X * gk;
	vector<Type> naok_d = X_d * gk;
	if (CppAD::Variable(logNaok))
		nll -= 0.5 * Sdim * logNaok - 0.5 * exp(logNaok) * GMRF(S).Quadform(gk);

	if (PAR6 != 0)
		nll -= sum(dnorm(e_o, Type(.0), exp(logSig), true));

	// Local variables
	Type lt;
	Type baseline;
	Type intrinsic_effect;
	int ds;
	int td;
	int tb;
	int id;
	int sx;
	int ms;

	// Likelihood
	for (int i = 0; i < l1.size(); i++)
	{
		ds = d2_d(i) - d1_d(i);
		lt = l1(i) + e_o(i);
		id = DUPE(i); // specimen id map
		sx = SEX(i);
		ms = MNS(i);
		if (ds < PAR4)
		{
			tb = d1_d(i);
			for (int j = 0; j < ds; j++)
			{
				intrinsic_effect = theta_l(MNS_map(tb + j)) + indl(id);
				lt = shark.growth_fn(lt, c + lsex(sx) + naoc_d(tb + j) + PAR5 * (intrinsic_effect), k + ksex(sx) + naok_d(tb + j) + (1 - PAR5) * intrinsic_effect, (Type)1 / 30.0, PAR1, exp(lognu));
			}
		}
		else
		{
			tb = t1(i);
			td = t2(i) - t1(i);
			if (td == 0)
			{
				// Released and recaptured in the same month
				intrinsic_effect = indl(id) + theta_l(ms);
				lt = shark.growth_fn(lt, c + lsex(sx) + naoc(tb) + PAR5 * intrinsic_effect, k + ksex(sx) + naok(tb) + (1 - PAR5) * intrinsic_effect, d1(i) - d2(i), PAR1, exp(lognu));
			}
			else
			{
				// Released and recaptured in different months
				intrinsic_effect = indl(id) + theta_l(ms);
				lt = shark.growth_fn(lt, c + lsex(sx) + naoc(tb) + PAR5 * intrinsic_effect, k + ksex(sx) + naok(tb) + (1 - PAR5) * intrinsic_effect, d1(i), PAR1, exp(lognu));
				for (int j = 1; j < td; j++)
				{
					intrinsic_effect = indl(id) + theta_l((ms + j) % 12);
					lt = shark.growth_fn(lt, c + lsex(sx) + naoc(tb + j) + PAR5 * intrinsic_effect, k + ksex(sx) + naok(tb + j) + (1 - PAR5) * intrinsic_effect, days(tb + j), PAR1, exp(lognu));
				}
				// Second partial month
				intrinsic_effect = indl(id) + theta_l((ms + td) % 12);
				lt = shark.growth_fn(lt, c + lsex(sx) + naoc(t2(i)) + PAR5 * intrinsic_effect, k + ksex(sx) + naok(t2(i)) + (1 - PAR5) * intrinsic_effect, 1 - d2(i), PAR1, exp(lognu));
			}
		}
		nll -= dnorm(l2(i), lt, exp(logSig), true);
	}
#ifdef _GRAPH
	// Prediction
	vector<Type> gam_k = prediction_design_matrix * gk;
	vector<Type> gam_c = prediction_design_matrix * gc;

	// Reference female with 150 cm length and one month interval
	int ref_sex = REFSEX;
	Type ref_len = REFLEN;
	Type ref_tim = REFTIM;
	int ref_mon = REFMON;

	vector<Type> gam_l(Pdim);
	intrinsic_effect = theta_l(ref_mon % 12);
	for (int i = 0; i < Pdim; i++)
	{
		gam_l(i) = shark.growth_fn(ref_len, c + lsex(ref_sex) + gam_c(i) + PAR5 * intrinsic_effect, k + ksex(ref_sex) + gam_k(i) + (1 - PAR5) * intrinsic_effect, ref_tim, PAR1, exp(lognu));
	}

	ADREPORT(gam_l);
	ADREPORT(gam_k);
	ADREPORT(gam_c);
#endif

#ifdef _PREDICT
	vector<Type> gam_kp = p_dm * gk;
	vector<Type> gam_cp = p_dm * gc;

	vector<Type> ref_p(days_p.size());
	Type tmp;
	for (int i = 0; i < days_p.size(); i++)
	{
		tmp = shark.growth_fn(ref_len, c + lsex(ref_sex) + gam_cp(i), k + ksex(ref_sex) + gam_kp(i), lag_ / 30., PAR1, exp(lognu));
		ref_p(i) = shark.growth_fn(tmp, c + lsex(ref_sex) + gam_cp(i + 1), k + ksex(ref_sex) + gam_kp(i + 1), days_p(i) - lag_ / 30., PAR1, exp(lognu));
	}

	ADREPORT(ref_p);
#endif

	return nll;
}
