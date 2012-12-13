#!/bin/sh
#nosetests --with-doctest --doctest-extension=".rst" README.rst
python -c "import sys; sys.path.insert(0, \".\"); import cruzdb; cruzdb.test()"
