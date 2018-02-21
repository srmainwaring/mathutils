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


def build_eigen(build_type):
    print("\n\n==============================================================================")
    print("Building Eigen library (header only)")
    print("==============================================================================")

    print("""
     _____ _                  _____ 
    | ____(_) __ _  ___ _ __ |___ / 
    |  _| | |/ _` |/ _ \ '_ \  |_ \ 
    | |___| | (_| |  __/ | | |___) |
    |_____|_|\__, |\___|_| |_|____/ 
             |___/
    """)

    os.chdir('eigen')

    try:
        os.mkdir('build')
    except OSError:
        pass

    os.chdir('build')
    call([cmake,
          '..'])

    os.chdir('../..')


def build_ceres_solver(build_type):
    print("\n\n==============================================================================")
    print("Building CERES-SOLVER library (shared)")
    print("==============================================================================")

    print("""
      ____                        ____        _                
     / ___|___ _ __ ___  ___     / ___|  ___ | |_   _____ _ __ 
    | |   / _ \ '__/ _ \/ __|____\___ \ / _ \| \ \ / / _ \ '__|
    | |__|  __/ | |  __/\__ \_____|__) | (_) | |\ V /  __/ |   
     \____\___|_|  \___||___/    |____/ \___/|_| \_/ \___|_|   
    """)


    eigen_include_dir = os.path.join(os.getcwd(), "eigen")
    # eigen_dir = os.path.join(os.getcwd(), "eigen", "build")

    os.chdir("ceres-solver")

    try:
        os.mkdir('build')
    except OSError:
        pass

    os.chdir("build")

    # FIXME: We do not provide at that time support for sparse linear solvers into ceres. Other dependencies should be added
    # but we have to check for the licences

    # print('\nEigen include dir located at : %s\n' % eigen_include_dir)
    # print('\nEigen3_DIR = %s' % eigen_dir)

    options = [
        "-DCMAKE_BUILD_TYPE=%s" % build_type,
        "-DEIGEN_INCLUDE_DIR=%s" % eigen_include_dir,
        # "-DEigen3_DIR=%s" % eigen_dir,
        "-DEIGEN_PREFER_EXPORTED_EIGEN_CMAKE_CONFIGURATION=FALSE",
        "-DMINIGLOG=ON",
        "-DCXX11=ON",  # FIXME: from the ceres documetation, should not be used with MSVC
        "-DEIGENSPARSE=ON",  # TODO: voir si le jeu de licence LGPL dont parle la doc ceres est toujours en cours... si oui, mettre OFF
        "-DBUILD_SHARED_LIBS=ON",
        "-DEXPORT_BUILD_DIR=ON",
        # "-DBUILD_DOCUMENTATION=ON"
    ]

    print('CMAKE is going to be applied on Ceres-solver with the following options :')
    print('-------------------------------------------------------------------------\n')
    for opt in options:
        print('\t' + opt)

    print("\nRunning CMAKE...")
    print("------------------\n")

    # Building command line
    cmd = [cmake] + options + ['..']
    # print(cmd)

    call(cmd)

    # Calling make
    call([make, "-j", str(nb_core)])

    os.chdir("../..")


def build():

    print("==============================================================================")
    print("Building thirdparty libraries for MathUtils project")
    print("==============================================================================\n\n")

    print("""
     __  __       _   _     _   _ _   _ _       _____         _ 
    |  \/  | __ _| |_| |__ | | | | |_(_) |___  |___ / _ __ __| |
    | |\/| |/ _` | __| '_ \| | | | __| | / __|   |_ \| '__/ _` |
    | |  | | (_| | |_| | | | |_| | |_| | \__ \  ___) | | | (_| |
    |_|  |_|\__,_|\__|_| |_|\___/ \__|_|_|___/ |____/|_|  \__,_|
                                                                
                      _         
     _ __   __ _ _ __| |_ _   _ 
    | '_ \ / _` | '__| __| | | |
    | |_) | (_| | |  | |_| |_| |
    | .__/ \__,_|_|   \__|\__, |
    |_|                   |___/
    """)

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

    call([git, "submodule", "init"])
    call([git, "submodule", "update"])
    
    print("-- DONE")

    build_eigen(build_type)
    build_ceres_solver('Release')  # For Ceres, there is no real interest in building it in Debug mode


if __name__ == "__main__":
    build()
