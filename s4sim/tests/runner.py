import unittest


def suite():
    """Returns unittest.TestSuite of tests.
    """
    from os.path import dirname

    py_dir = dirname(dirname(__file__))
    return unittest.defaultTestLoader.discover(py_dir, top_level_dir=dirname(py_dir))


def run():
    """Run all tests.
    """
    tests = suite()
    unittest.TextTestRunner(verbosity=2).run(tests)


if __name__ == "__main__":
    run()
