import unittest

import numpy as np

from scripts.analyze_inputs import describe_change
from scripts.analyze_interface_features import sequence_features
from scripts.compare_structures import kabsch


class DescribeChangeTests(unittest.TestCase):
    def test_replacement(self) -> None:
        result = describe_change("AAABBBCCC", "AAAXXXCCC")
        self.assertEqual(result["type"], "replacement")
        self.assertEqual((result["a_start"], result["a_end"]), (4, 6))

    def test_insertion(self) -> None:
        result = describe_change("AAACCC", "AAABBBCCC")
        self.assertEqual(result["type"], "insert")
        self.assertEqual(result["b_sequence"], "BBB")


class KabschTests(unittest.TestCase):
    def test_recovers_rigid_transform(self) -> None:
        reference = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        mobile = reference + np.array([5.0, -3.0, 2.0])
        rotation, translation = kabsch(mobile, reference)
        aligned = mobile @ rotation + translation
        np.testing.assert_allclose(aligned, reference, atol=1e-10)


class InterfaceFeatureTests(unittest.TestCase):
    def test_rac1b_insertion_charge(self) -> None:
        features = sequence_features("VGETYGKDITSRGKDKPIA")
        self.assertEqual(features["length"], 19)
        self.assertEqual(features["net_charge_proxy"], 1)
        self.assertEqual(features["positive_fraction"], 0.211)
        self.assertEqual(features["negative_fraction"], 0.158)


if __name__ == "__main__":
    unittest.main()
