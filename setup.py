#!/usr/bin/env python

from distutils.core import setup

setup(name='cruzdb',
      version='0.1',
      description='''Interface to UCSC genomic databases.
Also allows things like up/downstream queries''',
      long_description=open('README.rst').read(),
      author='Brent Pedersen',
      author_email='bpederse@gmail.com',
      url='https://github.com/brentp/cruzdb/',
      packages=['cruzdb', 'cruzdb.tests'],
      # can also use mysql-ctypes
      #requires=['mysql-python', 'sqlalchemy>=0.7'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'Intended Audience :: System Administrators',
          'License :: OSI Approved :: MIT License',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Topic :: Database',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],

 )
