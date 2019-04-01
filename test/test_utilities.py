""" :module TestUtilities: Hosting various utilities useful for testing. """

import os, shutil

def _remove_test_files(files):
    """ """
    """ Remove all files and directories listed.

    :param files: Files and directories to remove.
    :type files: list

    """

    # Loop over files
    for f in files:
        # Check if file.
        if os.path.isfile(f):
            os.remove(f)
        # Check if dir.
        elif os.path.isdir(f):
            shutil.rmtree(f)

def check_keys(test_class, expected_keys, dictionary):
    """ Check that all expected keys are contained in the passed items. """

    present_keys = dictionary.keys()

    for xk in expected_keys:
        test_class.assertIn(xk, present_keys)


