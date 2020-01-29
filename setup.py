#!/usr/bin/env python3
# Copyright (c) 2020 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""This script is the entry point for building, distributing and installing
this module using distutils/setuptools."""
import datetime
import distutils.command.build
import os
import pathlib
import platform
import re
import setuptools
import setuptools.command.build_ext
import setuptools.command.install
import subprocess
import sys
import sysconfig

# Check Python requirement
MAJOR = sys.version_info[0]
MINOR = sys.version_info[1]
if not (MAJOR >= 3 and MINOR >= 6):
    raise RuntimeError("Python %d.%d is not supported, "
                       "you need at least Python 3.6." % (MAJOR, MINOR))

# Working directory
WORKING_DIRECTORY = pathlib.Path(__file__).parent.absolute()


def build_dirname(extname=None):
    """Returns the name of the build directory"""
    extname = '' if extname is None else os.sep.join(extname.split(".")[:-1])
    return str(
        pathlib.Path(WORKING_DIRECTORY, "build",
                     "lib.%s-%d.%d" % (sysconfig.get_platform(), MAJOR, MINOR),
                     extname))


def execute(cmd):
    """Executes a command and returns the lines displayed on the standard
    output"""
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.stdout.read().decode()


def update_meta(path, version):
    """Updating the version number description in conda/meta.yaml."""
    with open(path, "r") as stream:
        lines = stream.readlines()
    pattern = re.compile(r'{% set version = ".*" %}')

    for idx, line in enumerate(lines):
        match = pattern.search(line)
        if match is not None:
            lines[idx] = '{%% set version = "%s" %%}\n' % version

    with open(path, "w") as stream:
        stream.write("".join(lines))


def revision():
    """Returns the software version"""
    os.chdir(WORKING_DIRECTORY)
    module = pathlib.Path(WORKING_DIRECTORY, 'src', 'pytide', 'version.py')
    stdout = execute("git describe --tags --dirty --long --always").strip()
    pattern = re.compile(r'([\w\d\.]+)-(\d+)-g([\w\d]+)(?:-(dirty))?')
    match = pattern.search(stdout)

    # If the information is unavailable (execution of this function outside the
    # development environment), file creation is not possible
    if not stdout:
        pattern = re.compile(r'return "(\d+\.\d+\.\d+)"')
        with open(module, "r") as stream:
            for line in stream:
                match = pattern.search(line)
                if match:
                    return match.group(1)
        raise AssertionError()

    assert match is not None, 'No GIT tag found'
    version = match.group(1)
    sha1 = match.group(3)

    stdout = execute("git log  %s -1 --format=\"%%H %%at\"" % sha1)
    stdout = stdout.strip().split()
    date = datetime.datetime.utcfromtimestamp(int(stdout[1]))

    # Conda configuration files are not present in the distribution, but only
    # in the GIT repository of the source code.
    meta = pathlib.Path(WORKING_DIRECTORY, 'conda', 'meta.yaml')
    if meta.exists():
        update_meta(meta, version)

    # Updating the version number description for sphinx
    conf = pathlib.Path(WORKING_DIRECTORY, 'docs', 'source', 'conf.py')
    with open(conf, "r") as stream:
        lines = stream.readlines()
    pattern = re.compile(r'(\w+)\s+=\s+(.*)')

    for idx, line in enumerate(lines):
        match = pattern.search(line)
        if match is not None:
            if match.group(1) == 'version':
                lines[idx] = "version = %r\n" % version
            elif match.group(1) == 'release':
                lines[idx] = "release = %r\n" % version
            elif match.group(1) == 'copyright':
                lines[idx] = "copyright = '(%s, CNES/CLS)'\n" % date.year

    with open(conf, "w") as stream:
        stream.write("".join(lines))

    # Finally, write the file containing the version number.
    with open(module, 'w') as handler:
        handler.write('''"""
# Copyright (c) 2020 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
Get software version information
================================
"""


def release() -> str:
    """Returns the software version number"""
    return "{version}"


def date() -> str:
    """Returns the creation date of this release"""
    return "{date}"
'''.format(version=version, date=date.strftime("%d %B %Y")))
    return version


class CMakeExtension(setuptools.Extension):
    """Python extension to build"""
    def __init__(self, name):
        super(CMakeExtension, self).__init__(name, sources=[])


class BuildExt(setuptools.command.build_ext.build_ext):
    """Build the Python extension using cmake"""

    #: Preferred C++ compiler
    CXX_COMPILER = None

    #: Preferred Eigen root
    EIGEN3_INCLUDE_DIR = None

    #: Preferred MKL root
    MKL_ROOT = None

    #: Run CMake to configure this project
    RECONFIGURE = None

    def run(self):
        """A command's raison d'etre: carry out the action"""
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    @staticmethod
    def eigen():
        """Get the default Eigen3 path in Anaconda's environnement."""
        eigen_include_dir = pathlib.Path(sys.prefix, "include", "eigen3")
        if eigen_include_dir.exists():
            return "-DEIGEN3_INCLUDE_DIR=" + str(eigen_include_dir)
        eigen_include_dir = pathlib.Path(sys.prefix, "Library", "include",
                                         "eigen3")
        if not eigen_include_dir.exists():
            eigen_include_dir = eigen_include_dir.parent
        if not eigen_include_dir.exists():
            raise RuntimeError(
                "Unable to find the Eigen3 library in the conda distribution "
                "used.")
        return "-DEIGEN3_INCLUDE_DIR=" + str(eigen_include_dir)

    @staticmethod
    def mkl():
        """Get the default MKL path in Anaconda's environnement."""
        mkl_header = pathlib.Path(sys.prefix, "include", "mkl.h")
        if mkl_header.exists():
            os.environ["MKLROOT"] = sys.prefix
            return
        mkl_header = pathlib.Path(sys.prefix, "Library", "include", "mkl.h")
        if mkl_header.exists():
            os.environ["MKLROOT"] = str(pathlib.Path(sys.prefix, "Library"))
            return
        raise RuntimeError(
            "Unable to find the MKL library in the conda distribution "
            "used.")

    @staticmethod
    def is_conda():
        """Detect if the Python interpreter is part of a conda distribution."""
        result = pathlib.Path(sys.prefix, 'conda-meta').exists()
        if not result:
            try:
                # pylint: disable=unused-import
                import conda
                # pylint: enable=unused-import
            except ImportError:
                result = False
            else:
                result = True
        return result

    def set_cmake_user_options(self):
        """Sets the options defined by the user."""
        is_conda = self.is_conda()
        result = []

        if self.CXX_COMPILER is not None:
            result.append("-DCMAKE_CXX_COMPILER=" + self.CXX_COMPILER)

        if self.EIGEN3_INCLUDE_DIR is not None:
            result.append("-DEIGEN3_INCLUDE_DIR=" + self.EIGEN3_INCLUDE_DIR)
        elif is_conda:
            result.append(self.eigen())

        if self.MKL_ROOT is not None:
            os.environ["MKLROOT"] = self.MKL_ROOT
        elif is_conda:
            self.mkl()

        return result

    def build_cmake(self, ext):
        """Execute cmake to build the Python extension"""
        # These dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(WORKING_DIRECTORY, self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = build_dirname(ext.name)

        cfg = 'Debug' if self.debug else 'Release'

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir),
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DCMAKE_PREFIX_PATH=" + sys.prefix
        ] + self.set_cmake_user_options()

        build_args = ['--config', cfg]

        if platform.system() != 'Windows':
            build_args += ['--', '-j%d' % os.cpu_count()]
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            if platform.system() == 'Darwin':
                cmake_args += ['-DCMAKE_OSX_DEPLOYMENT_TARGET=10.14']
        else:
            cmake_args += [
                '-G', 'Visual Studio 15 2017',
                '-DCMAKE_GENERATOR_PLATFORM=x64',
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                    cfg.upper(), extdir)
            ]
            build_args += ['--', '/m']
            if self.verbose:
                build_args += ['/verbosity:n']

        if self.verbose:
            build_args.insert(0, "--verbose")

        os.chdir(str(build_temp))

        # Has CMake ever been executed?
        if pathlib.Path(build_temp, "CMakeFiles",
                        "TargetDirectories.txt").exists():
            # The user must force the reconfiguration
            configure = self.RECONFIGURE is not None
        else:
            configure = True

        if configure:
            self.spawn(['cmake', str(WORKING_DIRECTORY)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.', '--target', 'core'] +
                       build_args)
        os.chdir(str(WORKING_DIRECTORY))


class Build(distutils.command.build.build):
    """Build everything needed to install"""
    user_options = distutils.command.build.build.user_options
    user_options += [
        ('eigen-root=', None, 'Preferred Eigen3 include directory'),
        ('cxx-compiler=', None, 'Preferred C++ compiler'),
        ('mkl-root=', None, 'Preferred MKL installation prefix'),
        ('reconfigure', None, 'Forces CMake to reconfigure this project')
    ]

    def initialize_options(self):
        """Set default values for all the options that this command supports"""
        super().initialize_options()
        self.cxx_compiler = None
        self.eigen_root = None
        self.mkl_root = None
        self.reconfigure = None

    def run(self):
        """A command's raison d'etre: carry out the action"""
        if self.cxx_compiler is not None:
            BuildExt.CXX_COMPILER = self.cxx_compiler
        if self.eigen_root is not None:
            BuildExt.EIGEN3_INCLUDE_DIR = self.eigen_root
        if self.reconfigure is not None:
            BuildExt.RECONFIGURE = True
        super().run()


def main():
    setuptools.setup(name='pytide',
                     version=revision(),
                     classifiers=[
                         "Development Status :: 3 - Alpha",
                         "Topic :: Scientific/Engineering :: Physics",
                         "License :: OSI Approved :: BSD License",
                         "Natural Language :: English",
                         "Operating System :: POSIX",
                         "Operating System :: MacOS",
                         "Operating System :: Microsoft :: Windows",
                         "Programming Language :: Python :: 3.6",
                         "Programming Language :: Python :: 3.7",
                         "Programming Language :: Python :: 3.8"
                     ],
                     description='Tidal constituents analysis in Python.',
                     url='https://github.com/CNES/pangeo-pytide',
                     author='CNES/CLS',
                     license="BSD License",
                     ext_modules=[CMakeExtension(name="pytide.core")],
                     setup_requires=[],
                     scripts=["src/scripts/mit_gcm_detiding.py"],
                     install_requires=["numpy"],
                     tests_require=["netCDF4", "numpy"],
                     package_dir={'': 'src'},
                     packages=setuptools.find_packages(where="src"),
                     cmdclass={
                         'build': Build,
                         'build_ext': BuildExt
                     },
                     zip_safe=False)


if __name__ == "__main__":
    main()
