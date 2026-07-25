/**
 * maths library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * @note the tlibs2 math library is based on code from my following projects:
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

#ifndef __TLIBS2_MATHS_H__
#define __TLIBS2_MATHS_H__


//#define USE_LAPACK 1
#define __TLIBS2_QR_METHOD__ 0


#include "maths/decls.h"
#include "maths/constants.h"
#include "maths/operators.h"
#include "maths/containers.h"
#include "maths/scalar.h"
#include "maths/ndim.h"
#include "maths/threedim.h"
#include "maths/projectors.h"
#include "maths/tensor.h"
#include "maths/solids.h"
#include "maths/polygon.h"
#include "maths/complex.h"
#include "maths/quaternion.h"
#include "maths/coordtrafos.h"
#include "maths/lapack.h"
#include "maths/qhull.h"
#include "maths/rtree.h"


#endif
