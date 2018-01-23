#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import sys
import shutil
sys.path.append("./which")
import which
from subprocess import call

abs_path_root = os.getcwd()


# For parallel make
try:
    import multiprocessing
    nb_core = multiprocessing.cpu_count()
except ImportError:
    nb_core = 1


try:
    make = which.which('make')
except which.WhichError:
    print "cannot find make"
    # make = '/usr/bin/make'

try:
    cmake = which.which('cmake')
except which.WhichError:
    print "cannot find cmake"
    cmake = '/home/frongere/mysofts/cmake-3.8.2-Linux-x86_64/bin/cmake'

try:
    git = which.which('git')
except which.WhichError:
    print "cannot find git"


def build_ceres_solver(build_type):
    print("\n\n==============================================================================")
    print("Building CERES-SOLVER library (shared)")
    print("==============================================================================")

    eigen_include_dir = os.path.join(os.getcwd(), "eigen")

    os.chdir("ceres-solver")

    try:
        os.mkdir('build')
    except OSError:
        pass

    os.chdir("build")

    # FIXME: We do not provide at that time support for sparse linear solvers into ceres. Other dependencies should be added
    # but we have to check for the licences

    call([cmake,
          "-DCMAKE_BUILD_TYPE=%s" % build_type,
          "-DEIGEN_INCLUDE_DIR=%s" % eigen_include_dir,
          "-DMINIGLOG=ON",
          "-DCXX11=ON",  # FIXME: from the ceres documetation, should not be used with MSVC
          "-DEIGENSPARSE=ON",  # TODO: voir si le jeu de licence LGPL dont parle la doc ceres est toujours en cours... si oui, mettre OFF
          "-DBUILD_SHARED_LIBS=ON",
          ".."])

    call([make, "-j%u" % nb_core])

    os.chdir("../..")


if __name__ == "__main__":

    print("==============================================================================")
    print("Building thirdparty libraries for MathUtils project")
    print("==============================================================================\n\n")

    if len(sys.argv) == 1:
        # By default, we build in Release mode
        build_type = "Debug"
    else:
        build_type = sys.argv[1]
        if build_type not in ("Debug", "Release", "RelWithDebInfo", "MinSizeRel"):
            raise NameError("Build type %s is not known. It has to be chosen among Debug, Release, RelWithDebInfo, MinSizeRel" % build_type)


    # Ensuring that thirdparty git submodules are up-to-date
    print("==============================================================================")
    print("Updating GIT SUBMODULES")
    print("==============================================================================")

    print("-- DONE")

    call([git, "submodule", "init"])
    call([git, "submodule", "update"])

    build_ceres_solver(build_type)