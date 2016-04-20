"""
jGEM ("j"unction and coverage based Gene and exon Extractor and Merger)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()


# Use build_ext from Cython if found
# Use build_ext from Cython
# command_classes = {}
try:
    # import Cython.Distutils
    # command_classes['build_ext'] = Cython.Distutils.build_ext
    from Cython.Build import cythonize
    have_cython = True
except:
    have_cython = False

try:
    import numpy
    have_numpy = True
except:
    have_numpy = False

def get_extension_modules():
    extensions = []        
    if have_numpy and have_cython:
        # jgem
        extensions.append( Extension( "jgem.cy.bw", [ "jgem/cy/bw.pyx" ], include_dirs=[numpy.get_include()]  ) )
        # Reading UCSC "big binary index" files
        extensions.append( Extension( "jgem.bxbbi.bpt_file", [ "jgem/bxbbi/bpt_file.pyx" ] ) )
        extensions.append( Extension( "jgem.bxbbi.cirtree_file", [ "jgem/bxbbi/cirtree_file.pyx" ] ) )
        extensions.append( Extension( "jgem.bxbbi.bbi_file", [ "jgem/bxbbi/bbi_file.pyx" ], include_dirs=[numpy.get_include()] ) )
        extensions.append( Extension( "jgem.bxbbi.bigwig_file", [ "jgem/bxbbi/bigwig_file.pyx" ], include_dirs=[numpy.get_include()] ) )
        extensions = cythonize(extensions)
    return extensions 

setup(
    name='jgem',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.9.0',

    description='jGEM ("j"unction and coverage based Gene and exon Extractor and Merger)',
    long_description=long_description,

    # The project's main homepage.
    url='',

    # Author details
    author='Ken Sugino',
    author_email='ken.sugino@gmail.com',

    # Choose your license
    license='BSD',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        #'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        #'Programming Language :: Python :: 3.2',
        #'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.5',
    ],

    # What does your project relate to?
    keywords='RNASeq bioinformatics',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy','scipy','pandas','matplotlib','cython','xlrd'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        #'dev': ['check-manifest'],
        #'test': ['coverage'],
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        #'sample': ['example_data/*'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        #'console_scripts': [
        #    'sample=sample:main',
        #],
    },
    ext_modules = get_extension_modules(),
    #cmdclass=command_classes,
)



