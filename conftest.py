# Copyright (c) 2022 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
import os
import pathlib
import sys
import sysconfig

# Check Python requirement
MAJOR = sys.version_info[0]
MINOR = sys.version_info[1]

# Working directory
WORKING_DIRECTORY = pathlib.Path(__file__).parent.absolute()


def build_dirname(extname=None):
    """Returns the name of the build directory."""
    extname = '' if extname is None else os.sep.join(extname.split('.')[:-1])
    path = pathlib.Path(
        WORKING_DIRECTORY, 'build',
        'lib.%s-%d.%d' % (sysconfig.get_platform(), MAJOR, MINOR), extname)
    if path.exists():
        return path
    return pathlib.Path(
        WORKING_DIRECTORY, 'build',
        'lib.{}-{}'.format(sysconfig.get_platform(),
                           sys.implementation.cache_tag), extname)


def push_front_syspath():
    """Add the build directory to the front of sys.path."""
    sys.path.insert(0, str(build_dirname().resolve()))


push_front_syspath()
