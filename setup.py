#!/usr/bin/env python

import os,sys, shlex, subprocess
from os.path import join
import platform
import distutils
from distutils import sysconfig
from setuptools import find_packages
from distutils.core import setup, Extension

cfgDict = distutils.sysconfig.get_config_vars()


progcompiler = 'g++ --std=c++0x'
sources='src/compute_similarity.cpp'

cmd = '%s %s -o bin/compute_similarity -I%s %s -L%s -lpython%s %s %s' % \
    (progcompiler, 
    cfgDict['LINKFORSHARED'],
    cfgDict['INCLUDEPY'],
    sources,
    cfgDict['LIBPL'],
    cfgDict['VERSION'], 
    cfgDict['LIBS'], 
    cfgDict['LIBM'])
    
# compile compute_similarity software
#print 'cmd = ', cmd 
cmd = shlex.split( cmd )
try :
    subprocess.check_call( cmd )
except subprocess.CalledProcessError as e :
    print >>sys.stderr, e.output
    print >>sys.stderr, e.cmd
    sys.exit( 1 )
except SystemExit as e :
    ext = e.code
    sys.exit( ext )
except :
    print >>sys.stderr, "Unexpected error:", sys.exc_info()[0]
    sys.exit(1)

try:
    import networkx
except ImportError:
    sys.stderr.write('networkx is not installed, you can find it at: '
                     'http://networkx.lanl.gov/\n')
    sys.exit(1)

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit(1)

try:
    import Pycluster
except ImportError:
    sys.stderr.write('pycluster is not installed you can find it at: '
                     'http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm#pycluster\n')
    sys.exit(1)

__version__ = ''
with open('lib/porthoDA/__init__.py') as inp:
  for line in inp:
      if line.startswith('__version__'):
          exec(line.strip())
          break
          
PACKAGES = [ 'porthoDA',]

PACKAGE_DATA = {
    'data': ['*.dat',]
}

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join('lib', *pkg.split('.'))

SCRIPTS = ['scripts/porthoDA']

if (platform.system() == 'Windows' or 
    len(sys.argv) > 1 and sys.argv[1] not in ('build', 'install')):
    for script in list(SCRIPTS):
        SCRIPTS.append(script + '.bat')

#extension = Extension('compute_similarity', 
#                      sources = ['src/compute_similarity.cpp'],
#                      extra_compile_args = ['--std=c++0x'] )

setup(
    name='porthoDA',
    version=__version__,
    author='Tristan Bitard-Feildel',
    author_email='t.bitard.feildel@uni-muenster.de',
    url='http://www.bornberglab.org/',
    description='A Python Package for proteinortho used with protein domain ',
    long_description="",
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    package_data=PACKAGE_DATA,
    license='GPL',
    keywords=('protein, domain, orthology, '
              'domain arrangement, domain similarity, '),
    classifiers=[
                 'Development Status :: 5 - Production/Stable',
                 'Intended Audience :: Education',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GPL License',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                ],
    scripts=SCRIPTS,
    requires= ['NumPy (>=1.6.2)',
               'Networkx (>=1.7)',
               'Matplotlib (>=1.1)',
               'Pycluster (>=1.49)',],
    provides=['porthoDA ({0:s})'.format(__version__)],
    #ext_modules=[extension],
    #entry_points = {
        #'console_scripts': ['watcher = porthoDA.apps.__init__::porthoDA_main'],
    #}
)




