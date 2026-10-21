"""Unit checks for D0 classification semantics; no study data are read."""

import unittest

from analyze_d0 import classify


class ClassifyTests(unittest.TestCase):
    def classify(self, **changes: object) -> str:
        args: dict[str, object] = {
            "all_rc0": True,
            "max_ll_gap": 0.0,
            "tol": 2e-4,
            "n_boundary": 0,
            "n_flat_eig": 0,
            "n_floor_under_marker": 0,
        }
        args.update(changes)
        return classify(**args)  # type: ignore[arg-type]

    def test_fit_failure_has_precedence(self) -> None:
        self.assertEqual(
            self.classify(all_rc0=False, max_ll_gap=1.0, n_boundary=1),
            "FIT-FAILED",
        )

    def test_start_sensitive_above_tolerance(self) -> None:
        self.assertEqual(self.classify(max_ll_gap=2.00001e-4), "START-SENSITIVE")

    def test_tolerance_boundary_is_inclusive(self) -> None:
        self.assertEqual(self.classify(max_ll_gap=2e-4), "NO-CAUTION-OBSERVED")

    def test_each_geometry_trigger(self) -> None:
        for change in (
            {"n_boundary": 1},
            {"n_flat_eig": 1},
            {"n_floor_under_marker": 1},
        ):
            with self.subTest(change=change):
                self.assertEqual(self.classify(**change), "GEOMETRY-CAUTION")

    def test_no_caution_observed(self) -> None:
        self.assertEqual(self.classify(), "NO-CAUTION-OBSERVED")


if __name__ == "__main__":
    unittest.main()
