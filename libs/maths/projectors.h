/**
 * maths library -- projection operators
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
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

#ifndef __TLIBS2_MATHS_PROJ_H__
#define __TLIBS2_MATHS_PROJ_H__

#include <cmath>
#include <tuple>
#include <vector>
#include <initializer_list>

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// cartesion projection operators
// ----------------------------------------------------------------------------

/**
 * matrix to project onto vector: P = |v><v|
 * from: |x'> = <v|x> * |v> = |v><v|x> = |v><v| * |x>
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_mat, class t_vec>
t_mat projector(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	if(is_normalised)
	{
		return tl2::outer<t_mat, t_vec>(vec, vec);
	}
	else
	{
		const auto len = tl2::norm<t_vec>(vec, false);
		return tl2::outer<t_mat, t_vec>(vec, vec) / len;
	}
}


/**
 * project vec1 onto vec2
 *
 * proj_op = |vec2><vec2¦/ len(vec2)^2,  len(vec2) = sqrt(<vec2|vec2>)
 * proj = proj_op * vec1 = |vec2> * <vec2|vec1> / <vec2|vec2>
 *
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_vec>
t_vec project(const t_vec& vec, const t_vec& vecProj, bool is_normalised = true)
requires is_vec<t_vec>
{
	if(is_normalised)
	{
		return tl2::inner<t_vec>(vecProj, vec) * vecProj;
	}
	else
	{
		const auto len = tl2::norm<t_vec>(vecProj, false);
		return tl2::inner<t_vec>(vecProj, vec) * vecProj / len;
	}
}


/**
 * project vector vec onto another vector vecProj
 * don't multiply with direction vector
 *
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real project_scalar(const t_vec& vec, const t_vec& vecProj, bool is_normalised = true)
requires is_vec<t_vec>
{
	if(is_normalised)
	{
		return tl2::inner<t_vec>(vecProj, vec);
	}
	else
	{
		const auto len = tl2::norm<t_vec>(vecProj);
		return tl2::inner<t_vec>(vecProj, vec) / len;
	}
}


/**
 * matrix to project onto orthogonal complement (plane perpendicular to vector): P = 1-|v><v|
 * from completeness relation: 1 = sum_i |v_i><v_i| = |x><x| + |y><y| + |z><z|
 *
 * @see (Arens 2015), p. 814 for the projection tensor
 */
template<class t_mat, class t_vec>
t_mat ortho_projector(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	const std::size_t size = vec.size();
	return unit<t_mat>(size) -
		tl2::projector<t_mat, t_vec>(vec, is_normalised);
}


/**
 * matrix to mirror on plane perpendicular to vector: P = 1 - 2*|v><v|
 * subtracts twice its projection onto the plane normal from the vector
 *
 * @see (Arens 2015), p. 710
 */
template<class t_mat, class t_vec>
t_mat ortho_mirror_op(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using T = typename t_vec::value_type;
	const std::size_t size = vec.size();

	return unit<t_mat>(size) -
		T(2)*tl2::projector<t_mat, t_vec>(vec, is_normalised);
}


/**
 * matrix to mirror [a, b, c, ...] into, e.g.,  [a, b', 0, 0]
 * @see (Scarpino 2011), p. 268
 */
template<class t_mat, class t_vec>
t_mat ortho_mirror_zero_op(const t_vec& vec, std::size_t row)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using T = typename t_vec::value_type;
	const std::size_t N = vec.size();

	t_vec vecSub = zero<t_vec>(N);
	for(std::size_t i=0; i<row; ++i)
		vecSub[i] = vec[i];

	// norm of rest vector
	T n = T(0);
	for(std::size_t i=row; i<N; ++i)
		n += vec[i]*vec[i];
	vecSub[row] = std::sqrt(n);

	const t_vec vecOp = vec - vecSub;

	// nothing to do -> return unit matrix
	if(equals_0<t_vec>(vecOp))
		return unit<t_mat>(vecOp.size(), vecOp.size());

	return ortho_mirror_op<t_mat, t_vec>(vecOp, false);
}


/**
 * project vector vec onto plane through the origin and perpendicular to vector vecNorm
 * (e.g. used to calculate magnetic interaction vector M_perp)
 *
 * @see (Stoecker 1999), Chapter "Analytische Geometrie"
 */
template<class t_vec>
t_vec ortho_project(const t_vec& vec, const t_vec& vecNorm, bool is_normalised = true)
requires is_vec<t_vec>
{
	return vec - tl2::project<t_vec>(vec, vecNorm, is_normalised);
}


/**
 * find orthonormal substitute basis for vector space (Gram-Schmidt algo)
 * remove orthogonal projections to all other basis vectors: |i'> = (1 - sum_{j<i} |j><j|) |i>
 *
 * @see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 * @see (Arens 2015), p. 744
 */
template<class t_vec,
	template<class...> class t_cont_in = std::initializer_list,
	template<class...> class t_cont_out = std::vector>
t_cont_out<t_vec> orthonorm_sys(const t_cont_in<t_vec>& sys)
requires is_vec<t_vec>
{
	t_cont_out<t_vec> newsys;
	newsys.reserve(sys.size());

	for(const t_vec& vecSys : sys)
	{
		t_vec vecOrthoProj = vecSys;

		// subtract projections to other basis vectors
		for(const t_vec& vecNewSys : newsys)
			vecOrthoProj -= project<t_vec>(vecSys, vecNewSys, true);

		// normalise
		vecOrthoProj /= norm<t_vec>(vecOrthoProj);
		newsys.emplace_back(std::move(vecOrthoProj));
	}

	return newsys;
}


/**
 * find orthonormal substitute basis for vector space (Gram-Schmidt algo)
 * remove orthogonal projections to all other basis vectors: |i'> = (1 - sum_{j<i} |j><j|) |i>
 *
 * @see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 * @see (Arens 2015), p. 744
 */
template<class t_mat, class t_vec>
t_mat orthonorm_sys(const t_mat& sys, bool use_col = true)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using size_t = decltype(sys.size1());
	t_mat newsys = sys;

	for(size_t colidx = 0; colidx < sys.size2(); ++colidx)
	{
		t_vec vecSys;
		if(use_col)
			vecSys = col<t_mat, t_vec>(sys, colidx);
		else
			vecSys = row<t_mat, t_vec>(sys, colidx);
		t_vec vecOrthoProj = vecSys;

		// subtract projections to other basis vectors
		for(size_t colidx2 = 0; colidx2 < colidx; ++colidx2)
		{
			t_vec vecNewSys;
			if(use_col)
				vecNewSys = col<t_mat, t_vec>(newsys, colidx2);
			else
				vecNewSys = row<t_mat, t_vec>(newsys, colidx2);
			vecOrthoProj -= project<t_vec>(vecSys, vecNewSys, true);
		}

		// normalise
		vecOrthoProj /= norm<t_vec>(vecOrthoProj);

		if(use_col)
			set_col<t_mat, t_vec>(newsys, vecOrthoProj, colidx);
		else
			set_row<t_mat, t_vec>(newsys, vecOrthoProj, colidx);
	}

	return newsys;
}
// ----------------------------------------------------------------------------

}

#endif
