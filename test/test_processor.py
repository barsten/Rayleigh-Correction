from unittest import TestCase

from raycorr.processor import main

class RayCorrTest(TestCase):
    def test_end_to_end(self):
        main(['../testdata/subset_0_of_MER_RR__1PTACR20050713_094325_000002592039_00022_17611_0000.dim'])
        # validate here
        # read output
        # read ref
        # compare