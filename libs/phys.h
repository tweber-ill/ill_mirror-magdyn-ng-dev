/**
 * physics library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * @note Forked on 7-Nov-2018 from my "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Additional functions were forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Further functions and updates forked on 1-Feb-2021 and 19-Apr-2021 from my privately developed "geo" and "misc" projects (https://github.com/t-weber/geo and https://github.com/t-weber/misc).
 *
 * @note for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * TAS-Paths
 * Copyright (C) 2021       Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_PHYS__
#define __TLIBS2_PHYS__

#include "units.h"
#include "maths.h"
#include "string.h"

#include <stdexcept>
#include <optional>
#include <boost/units/pow.hpp>


namespace tl2 {


// --------------------------------------------------------------------------------
// constants
// --------------------------------------------------------------------------------
// import scipy.constants as co
// E2KSQ = 2.*co.neutron_mass/(co.Planck/co.elementary_charge*1000./2./co.pi)**2. / co.elementary_charge*1000. * 1e-20
template<class T = double> constexpr T KSQ2E = T(0.5) * hbar<T>/angstrom<T>/m_n<T> * hbar<T>/angstrom<T>/meV<T>;
template<class T = double> constexpr T E2KSQ = T(1) / KSQ2E<T>;
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// tas calculations
// @see M. D. Lumsden, J. L. Robertson, and M. Yethiraj, J. Appl. Crystallogr. 38(3), pp. 405–411 (2005), doi: 10.1107/S0021889805004875.
// @see (Shirane 2002), Ch. 1.3
// --------------------------------------------------------------------------------

/**
 * angle between ki and kf in the scattering triangle (a4)
 * @returns nullopt if the angle can't be reached
 *
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * 2*<ki|kf> = ki^2 + kf^2 - Q^2
 * cos phi = (ki^2 + kf^2 - Q^2) / (2 ki kf)
 */
template<typename t_real = double>
std::optional<t_real> calc_tas_angle_ki_kf(
	t_real ki, t_real kf, t_real Q, t_real sense = 1)
{
	t_real c = (ki*ki + kf*kf - Q*Q) / (t_real(2)*ki*kf);
	if(std::abs(c) > t_real(1))
		return std::nullopt;
	return sense*std::acos(c);
}


/**
 * angle between ki and kf in the scattering triangle (a4)
 * (version with units)
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_angle_ki_kf(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf, const t_wavenumber<Sys, T>& Q,
	bool bPosSense = true)
{
	t_dimensionless<Sys, T> ttCos = (ki*ki + kf*kf - Q*Q)/(T(2.)*ki*kf);
	if(units::abs(ttCos) > T(1.))
		throw std::runtime_error("Scattering triangle not closed.");

	t_angle<Sys, T> tt = units::acos(ttCos);

	if(!bPosSense)
		tt = -tt;
	return tt;
}


/**
 * angle between ki and Q in the scattering triangle
 * @returns nullopt if the angle can't be reached
 *
 * |Q> = |ki> - |kf>
 * |kf> = |ki> + |Q>
 * kf^2 = ki^2 + Q^2 - 2*<ki|Q>
 * 2*<ki|Q> = ki^2 + Q^2 - kf^2
 * cos phi = (ki^2 + Q^2 - kf^2) / (2 ki*Q)
 */
template<typename t_real = double>
std::optional<t_real> calc_tas_angle_ki_Q(
	t_real ki, t_real kf, t_real Q, t_real sense = 1)
{
	t_real c = (ki*ki + Q*Q - kf*kf) / (t_real(2)*ki*Q);
	if(std::abs(c) > t_real(1))
		return std::nullopt;
	return sense*std::acos(c);
}


/**
 * angle between ki and Q in the scattering triangle
 * (version with units)
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_angle_ki_Q(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf,
	const t_wavenumber<Sys, T>& Q,
	bool bPosSense = true,
	bool bAngleOutsideTriag = false)
{
	t_angle<Sys, T> angle;

	if(Q*angstrom<T> == T(0.))
	{
		angle = pi<T>/T(2) * radians<T>;
	}
	else
	{
		auto c = (ki*ki - kf*kf + Q*Q) / (T(2.)*ki*Q);
		if(units::abs(c) > T(1.))
			throw std::runtime_error("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(bAngleOutsideTriag)
		angle = pi<T>*radians<T> - angle;
	if(!bPosSense)
		angle = -angle;

	return angle;
}


/**
 * angle between kf and Q in the scattering triangle
 * (version with units)
 *
 * Q_vec = ki_vec - kf_vec
 * ki_vec = Q_vec + kf_vec
 * ki^2 = Q^2 + kf^2 + 2Q kf cos th
 * cos th = (ki^2 - Q^2 - kf^2) / (2Q kf)
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_angle_kf_Q(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf,
	const t_wavenumber<Sys, T>& Q,
	bool bPosSense = true,
	bool bAngleOutsideTriag = true)
{
	t_angle<Sys, T> angle;

	if(Q*angstrom<T> == T(0.))
		angle = pi<T>/T(2) * radians<T>;
	else
	{
		auto c = (ki*ki - kf*kf - Q*Q) / (T(2.)*kf*Q);
		if(units::abs(c) > T(1.))
			throw std::runtime_error("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(!bAngleOutsideTriag)
		angle = pi<T>*radians<T> - angle;
	if(!bPosSense)
		angle = -angle;

	return angle;
}


/**
 * get length of Q
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * Q^2 = ki^2 + kf^2 - 2*ki*kf*cos(a4)
 */
template<typename t_real = double>
t_real calc_tas_Q_len(t_real ki, t_real kf, t_real a4)
{
	t_real Qsq = ki*ki + kf*kf - t_real(2)*ki*kf*std::cos(a4);
	return std::sqrt(Qsq);
}


/**
 * get length of Q (version with units)
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T>
calc_tas_Q_len(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf, const t_angle<Sys, T>& tt)
{
	t_dimensionless<Sys, T> ctt = units::cos(tt);
	decltype(ki*ki) Qsq = ki*ki + kf*kf - T(2.)*ki*kf*ctt;

	if(T(Qsq*angstrom<T>*angstrom<T>) < T(0.))
	{
		// TODO
		Qsq = -Qsq;
	}

	t_wavenumber<Sys, T> Q = my_units_sqrt<t_wavenumber<Sys, T>>(Qsq);
	return Q;
}


/**
 * get tas a3 and a4 angles
 * @return [a3, a4, distance of Q to the scattering plane]
 * @see M. D. Lumsden, et al., doi: 10.1107/S0021889805004875.
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::tuple<bool, t_real, t_real, t_real> calc_tas_a3a4(
	const t_mat& B, t_real ki_lab, t_real kf_lab,
	const t_vec& Q_rlu, const t_vec& orient_rlu, const t_vec& orient_up_rlu,
	t_real sample_sense = 1, t_real a3_offs = pi<t_real>)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	// metric from crystal B matrix
	t_mat G = metric<t_mat>(B);

	// length of Q vector
	t_real Q_len_lab = norm<t_mat, t_vec>(G, Q_rlu);

	// angle xi between Q and orientation reflex
	t_real xi = angle<t_mat, t_vec>(G, Q_rlu, orient_rlu);

	// sign/direction of xi
	t_vec xivec = cross<t_mat, t_vec>(G, orient_rlu, Q_rlu);
	t_real xidir = inner<t_mat, t_vec>(G, xivec, orient_up_rlu);
	if(xidir < t_real(0))
		xi = -xi;

	// angle psi between ki and Q
	std::optional<t_real> psi =
		calc_tas_angle_ki_Q<t_real>(ki_lab, kf_lab, Q_len_lab, sample_sense);
	if(!psi)
		return std::make_tuple(false, 0, 0, 0);

	// crystal and scattering angle
	t_real a3 = - *psi - xi + a3_offs;
	std::optional<t_real> a4 =
		calc_tas_angle_ki_kf<t_real>(ki_lab, kf_lab, Q_len_lab);
	if(!a4)
		return std::make_tuple(false, a3, 0, 0);
	*a4 *= sample_sense;

	// distance of Q to the scattering plane
	t_real dist_Q_plane = inner<t_mat, t_vec>(G, Q_rlu, orient_up_rlu);
	dist_Q_plane /= norm<t_mat, t_vec>(G, orient_up_rlu);

	return std::make_tuple(true, a3, *a4, dist_Q_plane);
}


/**
 * get hkl position of a tas
 * @return Q_rlu
 * @see M. D. Lumsden, et al., doi: 10.1107/S0021889805004875.
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::optional<t_vec> calc_tas_hkl(
	const t_mat& B, t_real ki_lab, t_real kf_lab, t_real Q_len_lab, t_real a3,
	const t_vec& orient_rlu, const t_vec& orient_up_rlu,
	t_real sample_sense = 1, t_real a3_offs = pi<t_real>)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	auto [Binv, ok] = inv<t_mat>(B);
	if(!ok)
		return std::nullopt;

	// angle psi between ki and Q
	std::optional<t_real> psi =
		calc_tas_angle_ki_Q<t_real>(ki_lab, kf_lab, Q_len_lab, sample_sense);
	if(!psi)
		return std::nullopt;

	// angle xi between Q and orientation reflex
	t_real xi = a3_offs - a3 - *psi;

	t_vec rotaxis_lab = B * orient_up_rlu;
	t_mat rotmat = rotation<t_mat, t_vec>(rotaxis_lab, xi, false);

	t_vec orient_lab = B * orient_rlu;
	t_vec Q_lab = rotmat * orient_lab;
	Q_lab /= norm<t_vec>(Q_lab);
	Q_lab *= Q_len_lab;

	t_vec Q_rlu = Binv * Q_lab;
	return Q_rlu;
}


/**
 * get a1 or a5 angle
 * @returns nullopt of the angle can't be reached
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 *
 * Bragg: n lam = 2d sin(theta)
 * n 2pi / k = 2d sin(theta)
 * n pi / k = d sin(theta)
 * theta = asin(n pi / (k d))
 */
template<class t_real = double>
std::optional<t_real> calc_tas_a1(t_real k, t_real d)
{
	t_real sintheta = pi<t_real> / (k*d);
	if(std::abs(sintheta) > t_real(1))
		return std::nullopt;
	return std::asin(sintheta);
}


/**
 * get a2 or a6 angle
 * (version with units)
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_a1(const t_wavenumber<Sys, T>& k,
	const t_length<Sys, T>& d, bool bPosSense = true)
{
	const T order = T(1.);
	t_length<Sys, T> lam = T(2.)*pi<T> / k;
	auto dS = order*lam/(T(2.)*d);
	if(std::abs(T(dS)) > T(1))
		throw std::runtime_error("Invalid twotheta angle.");

	t_angle<Sys, T> theta = units::asin(dS);
	if(!bPosSense)
		theta = -theta;
	return theta;
}


/**
 * get k from crystal angle
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 *
 * k = n pi / (d sin(theta))
 */
template<class t_real = double>
t_real calc_tas_k(t_real theta, t_real d)
{
	t_real sintheta = std::abs(std::sin(theta));
	return pi<t_real> / (d * sintheta);
}


/**
 * get k from crystal angle
 * (version with units)
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T> calc_tas_k(const t_angle<Sys, T>& _theta,
	const t_length<Sys, T>& d, bool bPosSense = true)
{
	t_angle<Sys, T> theta = _theta;
	if(!bPosSense)
		theta = -theta;

	const T order = T(1.);

	// https://en.wikipedia.org/wiki/Bragg%27s_law
	t_length<Sys, T> lam = T(2.)*d/order * units::sin(theta);
	t_wavenumber<Sys, T> k = T(2.)*pi<T> / lam;

	return k;
}


/**
 * get ki from kf and energy transfer
 */
template<class t_real = double>
t_real calc_tas_ki(t_real kf, t_real E)
{
	return std::sqrt(kf*kf + E2KSQ<t_real>*E);
}


/**
 * get kf from ki and energy transfer
 */
template<class t_real = double>
t_real calc_tas_kf(t_real ki, t_real E)
{
	return std::sqrt(ki*ki - E2KSQ<t_real>*E);
}


/**
 * get energy transfer from ki and kf
 */
template<class t_real>
t_real calc_tas_E(t_real ki, t_real kf)
{
	return (ki*ki - kf*kf) / E2KSQ<t_real>;
}


template<class Sys, class T = double>
t_energy<Sys, T> k2E(const t_wavenumber<Sys, T>& k)
{
	T dk = k*angstrom<T>;
	T dE = KSQ2E<T> * dk*dk;
	return dE * meV<T>;
}


template<class Sys, class T = double>
t_wavenumber<Sys, T> E2k(const t_energy<Sys, T>& _E, bool &bImag)
{
	bImag = (_E < T(0.)*meV<T>);
	t_energy<Sys, T> E = bImag ? -_E : _E;
	const T dE = E / meV<T>;
	const T dk = std::sqrt(E2KSQ<T> * dE);
	return dk / angstrom<T>;
}


/**
 * get energy transfer from ki and kf
 * (version with units)
 */
template<class Sys, class T = double>
t_energy<Sys, T> get_energy_transfer(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf)
{
	return k2E<Sys, T>(ki) - k2E<Sys, T>(kf);
}


/**
 * (hbar*ki)^2 / (2*mn)  -  (hbar*kf)^2 / (2mn)  =  E
 * 1) ki^2  =  +E * 2*mn / hbar^2  +  kf^2
 * 2) kf^2  =  -E * 2*mn / hbar^2  +  ki^2
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T> get_other_k(const t_energy<Sys, T>& E,
	const t_wavenumber<Sys, T>& kfix, bool bFixedKi)
{
	auto kE_sq = E*T(2.)*(m_n<T>/hbar<T>)/hbar<T>;
	if(bFixedKi)
		kE_sq = -kE_sq;

	auto k_sq = kE_sq + kfix*kfix;
	if(k_sq*angstrom<T>*angstrom<T> < T(0.))
		throw std::runtime_error("Scattering triangle not closed.");

	return my_units_sqrt<t_wavenumber<Sys, T>>(k_sq);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
/**
 * Bose distribution (occupation number including detailed balance)
 * @see (Shirane 2002), p. 28
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 *
 * bose(+E, T) / bose(-E, T) = [ 1/(exp(E/kT) - 1) + 1 ] / [ 1/(exp(E/kT) - 1) ]
 *                           = 1 + 1/1/(exp(E/kT) - 1)
 *                           = exp(E/kT)
 * which is the detailed balance, S(+Q, +E) / S(-Q, -E), see (Shirane 2002), p. 26.
 */
template<class t_real = double>
t_real bose(t_real E, t_real T)
{
	const t_real _kB = kB<t_real> * kelvin<t_real>/meV<t_real>;

	t_real n = t_real(1)/(std::exp(std::abs(E)/(_kB*T)) - t_real(1));
	if(E >= t_real(0))
		n += t_real(1);

	return n;
}


/**
 * Bose factor with a lower cutoff energy
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class t_real = double>
t_real bose_cutoff(t_real E, t_real T, t_real E_cutoff=t_real(0.02))
{
	t_real dB;

	E_cutoff = std::abs(E_cutoff);
	if(std::abs(E) < E_cutoff)
		dB = bose<t_real>(sign(E)*E_cutoff, T);
	else
		dB = bose<t_real>(E, T);

	return dB;
}


/**
 * Bose factor
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class Sys, class T = double>
T bose(const t_energy<Sys, T>& E, const t_temperature<Sys, T>& temp,
	t_energy<Sys, T> E_cutoff = -meV<T>)
{
	if(E_cutoff < T(0)*meV<T>)
		return bose<T>(T(E/meV<T>), T(temp/kelvin<T>));
	else
		return bose_cutoff<T>(T(E/meV<T>), T(temp/kelvin<T>),
			T(E_cutoff/meV<T>));
}


// ----------------------------------------------------------------------------
// polarisation
// ----------------------------------------------------------------------------

/**
 * polarisation density matrix
 *   (based on a proof from a lecture by P. J. Brown, 2006)
 *
 * eigenvector expansion of a state: |psi> = a_i |xi_i>
 * mean value of operator with mixed states:
 * <A> = p_i * <a_i|A|a_i>
 * <A> = tr( A * p_i * |a_i><a_i| )
 * <A> = tr( A * rho )
 * polarisation density matrix: rho = 0.5 * (1 + <P|sigma>)
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9
 * @see (Desktop Bronstein 2008), Ch. 21 (Zusatzkapitel.pdf), pp. 11-12 and p. 24
 */
template<class t_vec, class t_mat>
t_mat pol_density_mat(const t_vec& P, typename t_vec::value_type c=0.5)
requires is_vec<t_vec> && is_mat<t_mat>
{
	return (unit<t_mat>(2,2) + proj_su2<t_vec, t_mat>(P, true)) * c;
}


/**
 * Blume-Maleev equation
 * calculate equation indirectly with density matrix
 *   (based on a proof from a lecture by P. J. Brown, 2006)
 *
 * V   = N*1 + <Mperp|sigma>
 * I   = tr( <V|V> rho )
 * P_f = tr( <V|sigma|V> rho ) / I
 *
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_mat, class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec>
blume_maleev_indir(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_mat<t_mat> && is_vec<t_vec>
{
	// spin-1/2
	constexpr t_cplx c = 0.5;

	// vector of pauli matrices
	const auto sigma = su2_matrices<std::vector<t_mat>>(false);

	// density matrix
	const auto density = pol_density_mat<t_vec, t_mat>(P_i, c);

	// potential
	const auto V_mag = proj_su2<t_vec, t_mat>(Mperp, true);
	const auto V_nuc = N * unit<t_mat>(2);
	const auto V = V_nuc + V_mag;
	const auto VConj = herm(V);

	// scattering intensity
	t_cplx I = trace(VConj*V * density);

	// ------------------------------------------------------------------------
	// scattered polarisation vector
	const auto m0 = (VConj * sigma[0]) * V * density;
	const auto m1 = (VConj * sigma[1]) * V * density;
	const auto m2 = (VConj * sigma[2]) * V * density;

	t_vec P_f = create<t_vec>({ trace(m0), trace(m1), trace(m2) });
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * Blume-Maleev equation
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec>
blume_maleev(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_vec<t_vec>
{
	const t_vec MperpConj = conj(Mperp);
	const t_cplx NConj = std::conj(N);
	constexpr t_cplx imag(0, 1);

	t_cplx N2 = N * NConj;
	t_cplx M2 = inner<t_vec>(Mperp, Mperp);
	t_vec Mx2 = cross<t_vec>({ MperpConj, Mperp });

	// ------------------------------------------------------------------------
	// intensity
	// nuclear and magnetic non-chiral
	t_cplx I = N2 + M2;

	// magnetic chiral
	I += imag * inner<t_vec>(P_i, Mx2);

	// nuclear-magnetic
	t_cplx I_nm = N * inner<t_vec>(Mperp, P_i);
	I += I_nm + std::conj(I_nm);
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// polarisation vector
	// nuclear
	t_vec P_f = N2 * P_i;                           // rotates P

	// magnetic non-chiral
	t_vec rot_ch = Mperp * inner<t_vec>(Mperp, P_i);
	P_f += rot_ch + tl2::conj(rot_ch);              // rotates P
	P_f -= M2 * P_i;                                // rotates P

	// magnetic chiral
	P_f -= imag * Mx2;                              // creates P

	// nuclear-magnetic
	t_vec rot_nm = imag * NConj * cross<t_vec>({ Mperp, P_i });
	t_vec create_nm = NConj * Mperp;
	P_f += rot_nm + tl2::conj(rot_nm);              // rotates P
	P_f += create_nm + tl2::conj(create_nm);        // creates P
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * Blume-Maleev in tensor form
 *   (based on a lecture by P. J. Brown, 2006, 2009)
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_mat, class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_mat, t_vec, t_vec>
blume_maleev_tensor(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using t_real = typename t_cplx::value_type;

	const t_vec MperpConj = conj(Mperp);
	const t_cplx NConj = std::conj(N);

	t_cplx N2 = N * NConj;
	t_cplx M2 = inner<t_vec>(Mperp, Mperp);
	t_mat Mo = t_real(2) * outer<t_mat, t_vec>(Mperp, Mperp);
	t_vec NM = t_real(2) * N * MperpConj;

	auto [ Mor, Moi ] = split_cplx<t_mat, t_mat>(Mo);
	auto [ NMr, NMi ] = split_cplx<t_vec, t_vec>(NM);

	// cross product vector of imaginary component of Mo
	t_vec Moivec = create<t_vec>({ Moi(2, 1), Moi(0, 2), Moi(1, 0) });

	// ------------------------------------------------------------------------
	// intensity
	// nuclear-magnetic and magnetic chiral intensity
	t_cplx I = inner<t_vec>(P_i, NMr + Moivec);

	// nuclear and magnetic non-chiral intensity
	I += N2 + M2;
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// rotates polarisation
	// nuclear and magnetic non-chiral components
	t_mat Prot = diag<t_mat>({ N2 - M2, N2 - M2, N2 - M2 });

	// magnetic non-chiral component
	// Mor * P_i corresponds to rot_ch + tl2::conj(rot_ch) in the case above
	Prot += Mor;

	// nuclear-magnetic component
	Prot += skewsymmetric<t_mat>(NMi);

	Prot /= I;
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// creates polarisation (nuclear-magnetic and magnetic chiral components)
	// Moivec correspond to imag * Mx2 in the case above
	t_vec Pcreate = NMr - Moivec;
	Pcreate /= I;
	// ------------------------------------------------------------------------

	t_vec P_f = Prot * P_i + Pcreate;
	return std::make_tuple(I, Prot, Pcreate, P_f);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// term symbols
// ----------------------------------------------------------------------------
/**
 * convert [S, L, J] into a term symbol, (2S+1)^L_J
 * @see https://en.wikipedia.org/wiki/Term_symbol
 */
template<class t_real = double, class t_str = std::string>
t_str to_termsymbol(t_real S, t_real L, t_real J)
{
	static const std::vector<t_str> vecL = {
		"S", "P", "D", "F",
		"G", "H", "I", "K",
		"L", "M", "N", "O" };

	t_str strS = var_to_str<t_real, t_str>(t_real(2)*S + 1);
	t_str strL = vecL[std::size_t(L)];
	t_str strJ = var_to_str<t_real, t_str>(J);

	return strS + strL + strJ;
}


/**
 * convert a term symbol, (2S+1)^L_J, into [S, L, J]
 * @see https://en.wikipedia.org/wiki/Term_symbol
 */
template<class t_real = double, class t_str = std::string>
std::tuple<t_real, t_real, t_real> from_termsymbol(const t_str& term)
{
	using t_ch = typename t_str::value_type;
	static const std::unordered_map<t_ch, t_real> mapSubOrbitals =
	{
		{'S', 0.}, {'P', 1.}, {'D',  2.}, {'F',  3.},
		{'G', 4.}, {'H', 5.}, {'I',  6.}, {'K',  7.},
		{'L', 8.}, {'M', 9.}, {'N', 10.}, {'O', 11.},
	};

	t_real S0 = 1., L = 0., J = 0., div = 1.;
	t_ch cSub = 'S', ch = 0;

	std::istringstream istr(term);
	istr >> S0 >> cSub >> J >> ch;
	if(ch == '/')
	{
		istr >> div;
		J /= div;
	}

	auto iter = mapSubOrbitals.find(cSub);
	if(iter == mapSubOrbitals.end())
		throw std::runtime_error("Invalid orbital.");
	L = iter->second;

	return std::make_tuple((S0 - 1.)/2., L, J);
}


/**
 * transforms the electron configuration, e.g. 1s2 -> [1, 0, 2]
 * @return [n, l, #electrons]
 */
template<class t_str = std::string>
std::tuple<uint16_t, uint16_t, uint16_t>
from_electron_config(const t_str& ecfg)
{
	using t_ch = typename t_str::value_type;
	static const std::unordered_map<t_ch, uint16_t> mapSubOrbitals =
	{
		{'s', 0}, {'p', 1}, {'d',  2}, {'f',  3},
		{'g', 4}, {'h', 5}, {'i',  6}, {'k',  7},
		{'l', 8}, {'m', 9}, {'n', 10}, {'o', 11},
	};

	uint16_t n = 0, l = 0, num_E = 0;
	t_ch cSub = 's';

	std::istringstream istr(ecfg);
	istr >> n >> cSub >> num_E;
	auto iter = mapSubOrbitals.find(cSub);
	if(iter == mapSubOrbitals.end())
		throw std::runtime_error("Invalid orbital.");
	l = iter->second;

	return std::make_tuple(n, l, num_E);
}


/**
 * gets term symbol from the electron configuration
 * @return [S, L, J]
 */
template<class t_real = double, class t_str = std::string>
std::tuple<t_real, t_real, t_real>
hund(const t_str& ecfgs)
{
	std::vector<t_str> vecOrbitals;
	tl2::get_tokens<t_str, t_str>(ecfgs, " ,;", vecOrbitals);

	t_real S = 0, L = 0, J = 0;

	// all orbitals
	for(const t_str& ecfg : vecOrbitals)
	{
		auto [ n, l, e ] = from_electron_config(ecfg);
		auto [ _S, _L, _J ] = hund(l, e);

		S += _S;
		L += _L;
		J += _J;
	}

	return std::make_tuple(S, L, J);
}


/**
 * effective g factor
 * @see (Khomskii 2014), equ. (2.13)
 * @see https://en.wikipedia.org/wiki/Land%C3%A9_g-factor
 */
template<class T = double>
T eff_gJ(T S, T L, T J, T gL = T(1), T gS = T(2))
{
	T Jfact = T(2)*J*(J + T(1));
	if(equals_0(Jfact))
		return T(0);

	T g = T(0.5) * (gL + gS) -
		(S*(S + T(1)) - L*(L + T(1))) * (gL - gS) / Jfact;

	return g;
}


/**
 * effective magneton number in units of muB
 * @see (Khomskii 2014), p. 33
 */
template<class T = double>
T eff_magnetons(T gJ, T J)
{
	return gJ * std::sqrt(J * (J + T(1)));
}
// ----------------------------------------------------------------------------


}
#endif
