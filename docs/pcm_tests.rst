.. _PCM_Tests:

===========================
The pycoxmunk testing suite
===========================

Pycoxmunk comes with a series of tests that are designed to ensure correct functionality of the library is maintained
after changes to the codebase. When one or more of the tests fails, this may indicate a problem with changes to
pycoxmunk's code. Alternatively, failures may result from new functionality or improvements that result in better
output.

When contributing to `pycoxmunk`, please ensure that you run the tests before opening a pull request.

To run the tests, follow these steps:

#. Install `pytest`, preferably using conda/mamba: `mamba install pytest`
#. Navigate to the tests directory: `cd pycoxmunk/Tests/`
#. To run all tests, simply type: `pytest .`
#. To run only one test, type the name of that test: `pytest test_PixMask.py`

Each code file in the `pycoxmunk` library (beginning `CM_`) has an equivalent file (prefixed `test_`) in the `Tests`
directory. If you modify one of the code files, you can run the appropriate `test_` file to check that the code still
functions as intended. Depending on your modifications, you may also have to update the tests code, or add additional
tests to check the functionality that you have added. Coverage statistics are also collected for the tests, the aim of
which is to ensure that as much of `pycoxmunk`'s code as possible is covered by the tests and that there are no gaps
into which bugs could occur or problems arise.

Writing new tests
=================

If you need to write new tests, please follow the style of the existing test code. It is better to have numerous
functions that each test a separate part of your code than to have one large test function that tests the entire code at
once.
