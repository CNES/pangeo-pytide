#!/usr/bin/env python3
# Copyright (c) 2022 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""This script is the entry point for building, distributing and installing
this module using distutils/setuptools."""
from typing import ClassVar, List, Optional
import datetime
# The setuptools must be imported before distutils
import distutils.command.build
import os
import pathlib
import platform
import re
import subprocess
import sys
import sysconfig

import setuptools
import setuptools.command.build_ext
import setuptools.command.install

# Check Python requirement
MAJOR = sys.version_info[0]
MINOR = sys.version_info[1]
if not (MAJOR >= 3 and MINOR >= 6):
    raise RuntimeError("Python %d.%d is not supported, "
                       "you need at least Python 3.6." % (MAJOR, MINOR))

# Working directory
WORKING_DIRECTORY = pathlib.Path(__file__).parent.absolute()


def build_dirname(extname=None):
    """Returns the name of the build directory."""
    extname = '' if extname is None else os.sep.join(extname.split(".")[:-1])
    path = pathlib.Path(
        WORKING_DIRECTORY, "build",
        "lib.%s-%d.%d" % (sysconfig.get_platform(), MAJOR, MINOR), extname)
    if path.exists():
        return path
    return pathlib.Path(
        WORKING_DIRECTORY, "build",
        "lib.%s-%s" % (sysconfig.get_platform(), sys.implementation.cache_tag),
        extname)


def execute(cmd):
    """Executes a command and returns the lines displayed on the standard
    output."""
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stream = process.stdout
    assert stream is not None
    return stream.read().decode()


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
    """Returns the software version."""
    os.chdir(WORKING_DIRECTORY)
    module = pathlib.Path(WORKING_DIRECTORY, 'src', 'pytide', 'version.py')

    # If the ".git" directory exists, this function is executed in the
    # development environment, otherwise it's a release.
    if not pathlib.Path(WORKING_DIRECTORY, '.git').exists():
        pattern = re.compile(r'return "(\d+\.\d+\.\d+)"')
        with open(module, "r") as stream:
            for line in stream:
                match = pattern.search(line)
                if match:
                    return match.group(1)
        raise AssertionError()

    stdout = execute("git describe --tags --dirty --long --always").strip()
    pattern = re.compile(r'([\w\d\.]+)-(\d+)-g([\w\d]+)(?:-(dirty))?')
    match = pattern.search(stdout)
    if match is None:
        # No tag found, use the last commit
        pattern = re.compile(r'([\w\d]+)(?:-(dirty))?')
        match = pattern.search(stdout)
        assert match is not None, f"Unable to parse git output {stdout!r}"
        version = "0.0"
        sha1 = match.group(1)
    else:
        version = match.group(1)
        commits = int(match.group(2))
        sha1 = match.group(3)
        if commits != 0:
            version += f".dev{commits}"

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
        handler.write('''# Copyright (c) 2022 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
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
    """Python extension to build."""

    def __init__(self, name):
        super(CMakeExtension, self).__init__(name, sources=[])


class BuildExt(setuptools.command.build_ext.build_ext):
    """Build the Python extension using cmake."""

    #: Preferred C++ compiler
    CXX_COMPILER: ClassVar[Optional[str]] = None

    #: Preferred Eigen root
    EIGEN3_INCLUDE_DIR: ClassVar[Optional[str]] = None

    #: Selected CMAKE generator
    GENERATOR: ClassVar[Optional[str]] = None

    #: Preferred MKL root
    MKL_ROOT: ClassVar[Optional[str]] = None

    #: Run CMake to configure this project
    RECONFIGURE: ClassVar[Optional[bool]] = None

    def run(self):
        """Carry out the action."""
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    @staticmethod
    def eigen():
        """Get the default Eigen3 path in Anaconda's environment."""
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
        """Get the default MKL path in Anaconda's environment."""
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
                import conda  # noqa: F401

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
        elif is_conda and platform.system() != 'Darwin':
            self.mkl()

        return result

    def get_cmake_args(self, cfg: str, extdir: str) -> List[str]:
        """build cmake arguments.

        # Args:
        * `cfg`: config, one of {"debug", "release"}
        * `extdir`: output directory.
        """
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir),
            "-DCMAKE_PREFIX_PATH=" + sys.prefix,
            "-DPython_EXECUTABLE=" + sys.executable,
        ] + self.set_cmake_user_options()

        is_windows = platform.system() == "Windows"

        if self.GENERATOR is not None:
            cmake_args.append("-G" + self.GENERATOR)
        elif is_windows:
            cmake_args.append("-G" + "Visual Studio 16 2019")

        if is_windows:
            cmake_args += [
                "-DCMAKE_GENERATOR_PLATFORM=x64",
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(
                    cfg.upper(), extdir),
            ]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            if platform.system() == "Darwin":
                cmake_args += ["-DCMAKE_OSX_DEPLOYMENT_TARGET=10.14"]
        return cmake_args

    def get_build_args(self, cfg: str) -> List[str]:
        """make compiler build arguments.

        # Args:
        * `cfg`: config, one of {"debug", "release"}
        """
        build_args = ["--config", cfg]
        is_windows = platform.system() == "Windows"
        if is_windows:
            build_args += ['--', '/m']
            if self.verbose:  # type: ignore
                build_args += ["/verbosity:n"]
        else:
            build_args += ["--", "-j%d" % os.cpu_count()]
            if self.verbose:  # type: ignore
                build_args.insert(0, "--verbose")
        return build_args

    def build_cmake(self, ext):
        """execute cmake to build the python extension."""
        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(WORKING_DIRECTORY, self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = str(
            pathlib.Path(self.get_ext_fullpath(ext.name)).parent.resolve())

        cfg = 'Debug' if self.debug else 'Release'

        os.chdir(str(build_temp))

        # Has CMake ever been executed?
        if pathlib.Path(build_temp, "CMakeFiles",
                        "TargetDirectories.txt").exists():
            # The user must force the reconfiguration
            configure = self.RECONFIGURE is not None
        else:
            configure = True

        if configure:
            cmake_args = self.get_cmake_args(cfg, extdir)
            self.spawn(["cmake", str(WORKING_DIRECTORY)] + cmake_args)
        if not self.dry_run:  # type: ignore
            build_args = self.get_build_args(cfg)
            self.spawn(["cmake", "--build", ".", "--target", "core"] +
                       build_args)
        os.chdir(str(WORKING_DIRECTORY))


class Build(distutils.command.build.build):
    """Build everything needed to install."""
    user_options = distutils.command.build.build.user_options
    user_options += [
        ('cxx-compiler=', None, 'Preferred C++ compiler'),
        ('eigen-root=', None, 'Preferred Eigen3 include directory'),
        ('generator=', None, 'Selected CMake generator'),
        ('mkl-root=', None, 'Preferred MKL installation prefix'),
        ('reconfigure', None, 'Forces CMake to reconfigure this project')
    ]

    def initialize_options(self):
        """Set default values for all the options that this command
        supports."""
        super().initialize_options()
        self.cxx_compiler = None
        self.eigen_root = None
        self.generator = None
        self.mkl_root = None
        self.reconfigure = None

    def run(self):
        """Carry out the action."""
        if self.cxx_compiler is not None:
            BuildExt.CXX_COMPILER = self.cxx_compiler
        if self.mkl_root is not None:
            BuildExt.MKL_ROOT = self.mkl_root
        if self.eigen_root is not None:
            BuildExt.EIGEN3_INCLUDE_DIR = self.eigen_root
        if self.generator is not None:
            BuildExt.GENERATOR = self.generator
        if self.reconfigure is not None:
            BuildExt.RECONFIGURE = True
        super().run()


def long_description():
    """Reads the README file."""
    with open(pathlib.Path(WORKING_DIRECTORY, "README.md")) as stream:
        return stream.read()


def main():
    setuptools.setup(
        name='pytide',
        version=revision(),
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Topic :: Scientific/Engineering :: Physics",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: POSIX",
            "Operating System :: MacOS",
            "Operating System :: Microsoft :: Windows",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
        ],
        description='Tidal constituents analysis in Python.',
        url='https://github.com/CNES/pangeo-pytide',
        author='CNES/CLS',
        license="BSD License",
        ext_modules=[CMakeExtension(name="pytide.core")],
        long_description=long_description(),
        setup_requires=[],
        scripts=["src/scripts/mit_gcm_detiding.py"],
        install_requires=["numpy"],
        tests_require=["netCDF4", "numpy"],
        package_dir={'': 'src'},
        packages=setuptools.find_packages(where="src"),
        cmdclass={
            'build': Build,
            'build_ext': BuildExt
        },  # type: ignore
        zip_safe=False)


if __name__ == "__main__":
    main()
