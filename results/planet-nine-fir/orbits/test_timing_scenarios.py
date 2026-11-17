import copy, json, tempfile, unittest
from pathlib import Path
from timing_scenarios import SCENARIOS, scenario_inputs

class TimingScenarios(unittest.TestCase):
    def test_roles_and_input_immutable(self):
        original = json.loads(Path(__file__).with_name('orbit_input.json').read_text())
        saved = copy.deepcopy(original)
        for name in SCENARIOS:
            derived = scenario_inputs(original, name)
            self.assertEqual(derived['akari_midpoint_utc'], SCENARIOS[name][0])
            self.assertEqual(original, saved)
        self.assertIn('counterfactual', SCENARIOS['jun1'][1])
        self.assertIn('contributor unverified', SCENARIOS['dec1'][1])
        self.assertEqual(scenario_inputs(original, 'historical'), original)

if __name__ == '__main__': unittest.main()
