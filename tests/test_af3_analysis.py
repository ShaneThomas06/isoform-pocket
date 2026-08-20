"""Unit tests for AlphaFold Server analysis helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.analyze_af3_multiseed import mean_matrix, model_number
from scripts.map_p120_hotspot import distance, model_identity


class AlphaFoldAnalysisHelperTests(unittest.TestCase):
    def test_model_number(self) -> None:
        self.assertEqual(model_number(Path("job_summary_confidences_4.json")), 4)

    def test_mean_matrix(self) -> None:
        matrices = [
            [[0.0, 0.2], [0.2, 1.0]],
            [[0.0, 0.4], [0.4, 1.0]],
        ]
        result = mean_matrix(matrices)
        for observed_row, expected_row in zip(result, [[0.0, 0.3], [0.3, 1.0]]):
            for observed, expected in zip(observed_row, expected_row):
                self.assertAlmostEqual(observed, expected)

    def test_distance_ignores_confidence_field(self) -> None:
        self.assertAlmostEqual(distance((0.0, 0.0, 0.0, 20.0), (3.0, 4.0, 0.0, 90.0)), 5.0)

    def test_model_identity(self) -> None:
        path = Path("RAC1B_p120ARM350-824_seed303/extracted/job_model_2.cif")
        self.assertEqual(model_identity(path), ("RAC1B", 303, 2))


if __name__ == "__main__":
    unittest.main()
