import unittest

class RunSmashImportTest(unittest.TestCase):

    def test_can_be_imported(self):
        """Test that all imports in run_antismash.py can be imported"""
        import run_antismash
        # Silence pylint about unused imports. We just want to test if importing works, after all.
        run_antismash
