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




def build_ffts(build_type):
    print("\n\n==============================================================================")
    print("Building FFTS as a shared library")
    print("==============================================================================")

    os.chdir('ffts')

    try:
        os.mkdir("build")
    except OSError:
        pass

    os.chdir('build')

    call([cmake, '..',
          '-DENABLE_SHARED=ON',
          '-DCMAKE_BUILD_TYPE=%s' % build_type
          ])

    call([make, '-j', str(nb_core)])


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

    call([git, "submodule", "init"])
    call([git, "submodule", "update"])


    build_ffts(build_type)