#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest


PATH = pathlib.Path(__file__).parents[1] / "spikes/scripts/spikein/stable_seed.py"
SPEC = importlib.util.spec_from_file_location("stable_seed", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class StableSeedTests(unittest.TestCase):
    def test_frozen_value(self):
        self.assertEqual(
            MODULE.stable_seed(13, "spike-independent-v1", ["DRR127476", "Fnuc", "0.0001"]),
            873630871,
        )

    def test_scheduler_and_batch_are_not_inputs(self):
        components = ["DRR127476", "CRCpanel", "0.0001", "Fnuc"]
        first = MODULE.stable_seed(13, "spike-community-v1", components)
        second = MODULE.stable_seed(13, "spike-community-v1", list(components))
        self.assertEqual(first, second)

    def test_biological_identifiers_change_seed(self):
        one = MODULE.stable_seed(13, "spike-independent-v1", ["DRR127476", "Fnuc", "0.0001"])
        two = MODULE.stable_seed(13, "spike-independent-v1", ["DRR127477", "Fnuc", "0.0001"])
        self.assertNotEqual(one, two)


if __name__ == "__main__":
    unittest.main()
