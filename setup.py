#!/usr/bin/env python
import os
from distutils.core import setup


# stolen from mpld3
def get_version():
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(os.path.join("cruzdb", "__init__.py"), "r") as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")

setup(name='cruzdb',
      version=get_version(),
      description='''Interface to UCSC genomic databases.
Also allows things like up/downstream/k-nearest-neighbor queries and mirroring
of tables to local sqlite databases''',
      long_description=open('README.rst').read(),
      author='Brent Pedersen',
      author_email='bpederse@gmail.com',
      url='https://github.com/brentp/cruzdb/',
      packages=['cruzdb', 'cruzdb.tests'],
      # can also use mysql-ctypes
      #requires=['mysql-python', 'sqlalchemy>=0.7', 'six'],
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
