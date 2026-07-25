/**
 * maths library -- constants
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

#ifndef __TLIBS2_MATHS_CONSTS_H__
#define __TLIBS2_MATHS_CONSTS_H__

#include <cmath>
#include <numbers>

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------
template<typename T = double> constexpr T pi{std::numbers::pi_v<T>};
template<typename T = double> constexpr T golden{std::numbers::phi_v<T>};

template<typename INT = int> bool is_even(INT i) { return (i%2 == 0); }
template<typename INT = int> bool is_odd(INT i) { return !is_even<INT>(i); }

template<class T = double> constexpr T r2d(T rad) { return rad/pi<T>*T(180); }     // rad -> deg
template<class T = double> constexpr T d2r(T deg) { return deg/T(180)*pi<T>; }     // deg -> rad
}

#endif
