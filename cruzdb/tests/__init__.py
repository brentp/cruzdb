import os

def load_tests(loader, standard_tests, pattern):
    # top level directory cached on loader instance
    this_dir = os.path.dirname(__file__)
    package_tests = loader.discover(start_dir=this_dir, pattern=pattern)
    standard_tests.addTests(package_tests)
    return standard_tests

def test():
    import unittest
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite = load_tests(loader, suite, "test*")
    unittest.TextTestRunner().run(suite)

if __name__ == "__main__":
    test()
