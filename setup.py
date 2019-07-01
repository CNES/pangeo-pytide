#!/usr/bin/env python3
# Copyright (c) 2019 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
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


def execute(cmd):
    """Executes a command and returns the lines displayed on the standard
    output"""
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.stdout.read().decode()


def revision():
    """Returns the software version"""
    cwd = pathlib.Path().absolute()
    module = os.path.join(cwd, 'src', 'pytide', 'version.py')
    stdout = execute("git describe --tags --dirty --long --always").strip()
    pattern = re.compile(r'([\w\d\.]+)-(\d+)-g([\w\d]+)(?:-(dirty))?')
    match = pattern.search(stdout)

    # If the information is unavailable (execution of this function outside the
    # development environment), file creation is not possible
    if not stdout:
        pattern = re.compile(r'\s+result = "(.*)"')
        with open(module, "r") as stream:
            for line in stream:
                match = pattern.search(line)
                if match:
                    return match.group(1)
        raise AssertionError()

    # TODO: Delete this test which will be useless after the creation of the
    # first tag
    if match is None:
        pattern = re.compile(r'([\w\d]+)(?:-(dirty))?')
        match = pattern.search(stdout)
        version = "0.1"
        sha1 = match.group(1)
    else:
        version = match.group(1)
        sha1 = match.group(3)

    stdout = execute("git log  %s -1 --format=\"%%H %%at\"" % sha1)
    stdout = stdout.strip().split()
    date = datetime.datetime.fromtimestamp(int(stdout[1]))
    sha1 = stdout[0]

    # Updating the version number description in "meta.yaml"
    meta = os.path.join(cwd, 'conda', 'meta.yaml')
    with open(meta, "r") as stream:
        lines = stream.readlines()
    pattern = re.compile(r'\s+version:\s+(.*)')
    for idx, line in enumerate(lines):
        match = pattern.search(line)
        if match is not None:
            lines[idx] = '  version: %s\n' % version
    with open(meta, "w") as stream:
        stream.write("".join(lines))

    # Finally, write the file containing the version number.
    with open(module, 'w') as handler:
        handler.write('''"""
Get software version information
================================
"""


def release(full: bool = False) -> str:
    """Returns the software version number"""
    # {sha1}
    result = "{version}"
    if full:
        result += " ({date})"
    return result
'''.format(sha1=sha1, version=version, date=date.strftime("%d %B %Y")))
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

    def run(self):
        """A command's raison d'etre: carry out the action"""
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    @staticmethod
    def eigen():
        """Get the default Eigen3 path in Anaconda's environnement."""
        eigen_include_dir = os.path.join(sys.prefix, "include", "eigen3")
        if os.path.exists(eigen_include_dir):
            return "-DEIGEN3_INCLUDE_DIR=" + eigen_include_dir
        eigen_include_dir = os.path.join(sys.prefix, "Library", "include")
        if not os.path.exists(eigen_include_dir):
            raise RuntimeError(
                "Unable to find the Eigen3 library in the conda distribution "
                "used.")
        return "-DEIGEN3_INCLUDE_DIR=" + eigen_include_dir

    @staticmethod
    def is_conda():
        """Detect if the Python interpreter is part of a conda distribution."""
        result = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))
        if not result:
            return 'conda-bld' in sysconfig.get_config_var("abs_srcdir")
        return result

    def build_cmake(self, ext):
        """Execute cmake to build the Python extension"""
        cwd = pathlib.Path().absolute()

        # These dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(
            ext.name)).absolute().parent

        cfg = 'Debug' if self.debug else 'Release'

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir),
            "-DPYTHON_EXECUTABLE=" + sys.executable
        ]

        if self.CXX_COMPILER is not None:
            cmake_args.append("-DCMAKE_CXX_COMPILER=" + self.CXX_COMPILER)

        is_conda = self.is_conda()

        if self.EIGEN3_INCLUDE_DIR is not None:
            cmake_args.append("-DEIGEN3_INCLUDE_DIR=" +
                              self.EIGEN3_INCLUDE_DIR)
        elif is_conda:
            cmake_args.append(self.eigen())

        build_args = ['--config', cfg]

        if platform.system() != 'Windows':
            build_args += ['--', '-j%d' % os.cpu_count()]
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            if platform.system() == 'Darwin':
                cmake_args += ['-DCMAKE_OSX_DEPLOYMENT_TARGET=10.14']
        else:
            cmake_args += [
                '-DCMAKE_GENERATOR_PLATFORM=x64',
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                    cfg.upper(), extdir)
            ]
            build_args += ['--', '/m']

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.', '--target', 'core'] +
                       build_args)
        os.chdir(str(cwd))


class Build(distutils.command.build.build):
    """Build everything needed to install"""
    user_options = distutils.command.build.build.user_options
    user_options += [('eigen-root=', 'e',
                      'Preferred Eigen3 include directory'),
                     ('cxx-compiler=', 'x', 'Preferred C++ compiler')]

    def initialize_options(self):
        """Set default values for all the options that this command supports"""
        super().initialize_options()
        self.cxx_compiler = None
        self.eigen_root = None

    def run(self):
        """A command's raison d'etre: carry out the action"""
        if self.cxx_compiler is not None:
            BuildExt.CXX_COMPILER = self.cxx_compiler
        if self.eigen_root is not None:
            BuildExt.EIGEN3_INCLUDE_DIR = self.eigen_root
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
                         "Programming Language :: Python :: 3.7"
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
