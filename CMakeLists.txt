#
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2022
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# magpie & mag-core
# Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

cmake_minimum_required(VERSION 3.10)
project(magpie)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra -Weffc++)

#add_compile_options(-coverage)
#link_libraries(--coverage)


# -----------------------------------------------------------------------------
#options
# -----------------------------------------------------------------------------
option(USE_QT      "use GUI"           TRUE)
option(USE_HDF     "use HDF"           TRUE)
option(USE_MINUIT  "use Minuit"        TRUE)
# -----------------------------------------------------------------------------


list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/setup/cmake"
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_compile_options(-g -ggdb)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)

	message("Debug building for ${CMAKE_SYSTEM_NAME} systems.")

elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
	add_compile_options(-g -ggdb3 -pg -O0)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)

	message("Profile building for ${CMAKE_SYSTEM_NAME} systems.")

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_compile_options("-DNDEBUG")
	add_compile_options("-Wno-#pragma-messages")

	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		# see: https://gcc.gnu.org/onlinedocs/gcc/Developer-Options.html
		add_compile_options(-fcompare-debug-second)
	endif()

	message("Release building for ${CMAKE_SYSTEM_NAME} systems.")
endif()


set(MINGW_EXTRALIBS)
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	# pass linker --subsystem option
	add_compile_options(-Wl,--subsystem,windows)
	set(MINGW_EXTRALIBS "ws2_32")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_options(-mmacosx-version-min=11.0)
endif()



# -----------------------------------------------------------------------------
# packages
# -----------------------------------------------------------------------------
find_package(Threads REQUIRED)

find_package(Boost REQUIRED COMPONENTS program_options)
set(Boost_USE_MULTITHREADED TRUE)
add_definitions(${Boost_CXX_FLAGS})

find_package(Gemmi REQUIRED)

find_package(GFortran)
find_package(Lapacke REQUIRED)
add_definitions(-DUSE_LAPACK=1)
include_directories("${Lapacke_INCLUDE_DIRS}")

find_package(Qhull REQUIRED)
if(Qhull_FOUND)
	add_definitions(-DUSE_QHULL)
	include_directories("${Qhull_INCLUDE_DIRS}")
endif()

if(USE_MINUIT)
	find_package(Minuit2)
endif()

if(Minuit2_FOUND)
	message("Using Minuit 2.")
	include_directories("${Minuit2_INCLUDE_DIRS}")
	add_definitions(-D__TLIBS2_MAGDYN_USE_MINUIT__)
else()
	message("Disabling Minuit support.")
	set(Minuit2_LIBRARIES "")
endif()

if(USE_HDF)
	find_package(HDF5 COMPONENTS CXX)
endif()
if(HDF5_FOUND)
	message("Using HDF5 version ${HDF5_VERSION}.")
	include_directories("${HDF5_CXX_INCLUDE_DIRS}")
	add_definitions(-DUSE_HDF5)
else()
	message("Disabling HDF5 file support.")
	set(HDF5_CXX_LIBRARIES "")
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# building
# -----------------------------------------------------------------------------
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/gui")

include_directories(
	SYSTEM "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.." "."
	"${PROJECT_SOURCE_DIR}/ext" "${PROJECT_SOURCE_DIR}/../../ext/qcp"
	gemmi::headers
)


if(USE_QT)
	message("Using Qt.")

	find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGLWidgets PrintSupport Svg)
	include_directories("${Qt6Core_INCLUDE_DIRS}/..")
	list(APPEND QtLibraries Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGLWidgets Qt6::PrintSupport Qt6::Svg)

	set(CMAKE_AUTOUIC TRUE)
	set(CMAKE_AUTOMOC TRUE)

	find_package(QCP)
	include_directories("${QCP_INCLUDE_DIRS}")

	add_executable(magpie
		# main program
		gui/startup/main.cpp gui/startup/main_gui.cpp
		gui/magdyn.cpp gui/magdyn.h
		gui/dialogs.cpp
		gui/structure.cpp gui/tables.cpp
		gui/config.cpp gui/file.cpp
		gui/export/sunny.cpp gui/export/sw.cpp
		gui/export/script.cpp gui/export/sqe.cpp
		gui/kernel.cpp gui/dispersion.cpp
		gui/defs.cpp gui/defs.h gui/gui_defs.h

		# widgets
		gui/widgets/graph.cpp gui/widgets/graph.h
		gui/widgets/siteitem.h

		# panels
		gui/panels/sites.cpp gui/panels/terms.cpp
		gui/panels/coordinates.cpp gui/panels/reciprocal.cpp
		gui/panels/sample.cpp gui/panels/environment.cpp
		gui/panels/dispersion.cpp gui/panels/hamilton.cpp
		gui/panels/variables.cpp gui/panels/export.cpp

		# cli tools
		gui/startup/healthcheck.cpp
		gui/startup/benchmark.cpp

		# bz calculator
		gui/dialogs/bz/bz_lib.h gui/dialogs/bz/plot_cut.h
		gui/dialogs/bz/plot.cpp gui/dialogs/bz/plot.h
		gui/dialogs/bz/bz.cpp gui/dialogs/bz/bz_file.cpp gui/dialogs/bz/bz_calc.cpp
		gui/dialogs/bz/bz_ops.cpp gui/dialogs/bz/bz.h
		gui/dialogs/bz/bz_conf.cpp gui/dialogs/bz/bz_conf.h
		gui/dialogs/bz/globals.cpp gui/dialogs/bz/globals.h
		gui/dialogs/bz/bz_formulas.cpp

		# dialogs
		gui/dialogs/struct_plot.cpp gui/dialogs/struct_plot.h
		gui/dialogs/ground_state.cpp gui/dialogs/ground_state.h
		gui/dialogs/topology.cpp gui/dialogs/topology.h
		gui/dialogs/diff.cpp gui/dialogs/diff.h
		gui/dialogs/powder.cpp gui/dialogs/powder.h
		gui/dialogs/ffact.cpp gui/dialogs/ffact.h
		gui/dialogs/dispersion3d.cpp gui/dialogs/dispersion3d_file.cpp
		gui/dialogs/dispersion3d_calc.cpp gui/dialogs/dispersion3d_plot.cpp
		gui/dialogs/dispersion3d.h
		gui/dialogs/table_import.cpp gui/dialogs/table_import.h
		gui/dialogs/trafos.cpp gui/dialogs/trafos.h
		gui/dialogs/pol.cpp gui/dialogs/pol.h
		gui/dialogs/assign.cpp gui/dialogs/assign.h
		gui/dialogs/matrixelems.cpp gui/dialogs/matrixelems.h
		gui/dialogs/notes.cpp gui/dialogs/notes.h
		gui/dialogs/infos.cpp gui/dialogs/infos.h
		gui/dialogs/glinfos.cpp gui/dialogs/glinfos.h
		gui/dialogs/settings.h

		# library
		libs/magdyn.h
		libs/qt/gl.cpp libs/qt/gl.h
		libs/qt/glplot.cpp libs/qt/glplot.h
		libs/qt/recent.cpp libs/qt/recent.h
		libs/symops.h

		# plotting tools
		gui/dialogs/plot2d/plot2d.cpp gui/dialogs/plot2d/plot2d.h
		gui/dialogs/plot3d/plot3d.cpp gui/dialogs/plot3d/plot3d_file.cpp
		gui/dialogs/plot3d/plot3d_calc.cpp gui/dialogs/plot3d/plot3d_render.cpp
		gui/dialogs/plot3d/plot3d_grid.cpp
		gui/dialogs/plot3d/plot3d.h

#		../ext/qcp/qcustomplot.cpp
	)

	target_link_libraries(magpie
		Threads::Threads
		${QCP_LIBRARIES}
		${QtLibraries}
		${Boost_LIBRARIES}
		${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES} ${GFortran_LIBRARIES}
		${HDF5_CXX_LIBRARIES}
		${Minuit2_LIBRARIES}
		${MINGW_EXTRALIBS}
		${Qhull_LIBRARIES}
		gemmi::gemmi_cpp

#		/usr/local/lib/liblapack.a /usr/local/lib/libblas.a libgfortran.so    # for static lapack(e)
	)
else()
	message("Compiling without Qt.")
	add_definitions(-DDONT_USE_QT)

	add_executable(magpie
		gui/startup/main.cpp
		gui/defs.cpp gui/defs.h

		# cli tools
		gui/startup/healthcheck.cpp
		gui/startup/benchmark.cpp

		# library
		libs/magdyn.h
	)

	target_link_libraries(magpie
		Threads::Threads
		${Boost_LIBRARIES}
		${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES} ${GFortran_LIBRARIES}
		${HDF5_CXX_LIBRARIES}
		${Minuit2_LIBRARIES}
		${MINGW_EXTRALIBS}
		gemmi::gemmi_cpp
	)
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# rpath
# -----------------------------------------------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_link_options(magpie
		# create an __info_plist section in the binary
		PRIVATE LINKER:-sectcreate,__TEXT,__info_plist,${PROJECT_SOURCE_DIR}/setup/osx/Info.plist
	)

	add_custom_command(TARGET magpie POST_BUILD
		COMMAND install_name_tool -add_rpath /usr/local/lib $<TARGET_FILE:magpie> || true
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# installation and packaging
# -----------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS magpie RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${PROJECT_SOURCE_DIR}/../res/magpie.svg DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/magpie/res)
install(FILES ${PROJECT_SOURCE_DIR}/../AUTHORS DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/magpie)
install(FILES ${PROJECT_SOURCE_DIR}/../LICENSE DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/magpie)
install(FILES ${PROJECT_SOURCE_DIR}/../LICENSES DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/magpie)
install(FILES ${PROJECT_SOURCE_DIR}/../LITERATURE DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/magpie)
install(FILES ${PROJECT_SOURCE_DIR}/../CITATION.cff DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/magpie)
install(FILES ${PROJECT_SOURCE_DIR}/../setup/deb/magpie.desktop DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/applications)


set(CPACK_PACKAGE_NAME "Magpie")
set(CPACK_PACKAGE_VENDOR "T. Weber, ILL, Grenoble")
set(CPACK_PACKAGE_VERSION "1.0")
set(CPACK_PACKAGE_DESCRIPTION "A Graphical Magnon Software")
set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/res/magpie.svg")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_STRIP_FILES True)

include(CPack)
# -----------------------------------------------------------------------------
