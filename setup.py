import distutils.command.build
import pathlib
import platform
import setuptools
import setuptools.command.build_ext
import setuptools.command.install
import os
import sys
import sysconfig


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
                     version='0.1',
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
                     install_requires=["numpy"],
                     test_requires=["netCDF4", "numpy"],
                     package_dir={'': 'src'},
                     packages=setuptools.find_packages(where="src"),
                     cmdclass={
                         'build': Build,
                         'build_ext': BuildExt
                     },
                     zip_safe=False)


if __name__ == "__main__":
    main()
