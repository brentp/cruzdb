#!/bin/sh
nosetests --with-doctest --doctest-extension=".rst" README.rst
PYTHONPATH=. python cruzdb/tests/__init__.py
# or python -c "import cruzdb; cruzdb.test()"
