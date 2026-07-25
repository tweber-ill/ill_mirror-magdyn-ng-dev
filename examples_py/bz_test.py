#
# brillouin zone scripting test
# @author Tobias Weber <tweber@ill.fr>
# @date April-2023
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# magpie
# Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

import sys
import os

sys.path.append(os.getcwd())


import magpy

bz = magpy.BZCalcD()
bz.SetEps(1e-6)

bz.SetCrystal(5, 5, 5, 90, 90, 90)

num_ops = bz.SetSymOpsFromSpaceGroup("F d -3 m")
print("Using %d centring symops." % num_ops)

num_peaks = bz.CalcPeaks(4, True)
print("Using %d reflections." % num_peaks)

calc_ok = bz.CalcBZ()
if calc_ok:
	print("Brillouin zone calculation successful.")
else:
	print("Brillouin zone calculation failed.")

if calc_ok:
	print("\nJSON Output:")
	json = bz.PrintJSON(6)
	print(json)
