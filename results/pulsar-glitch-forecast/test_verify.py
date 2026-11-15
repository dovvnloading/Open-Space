"""Regression checks for study population, split integrity and external measurement replay."""
import io
import json
import unittest
from contextlib import redirect_stdout

import verify


class PartitionTests(unittest.TestCase):
    def events(self, count):
        return [verify.Event(i, 100. * i, 1., float(i % 5 + 1)) for i in range(count)]

    def test_missing_size_can_close_but_cannot_start_an_interval(self):
        events = self.events(17)
        events[5] = verify.Event(5, 500., 1., None)
        partitions, count = verify.derive_partitions({'J': events}, {})
        self.assertEqual(count, 15)
        ordered = [p for half in partitions['J'] for p in half]
        self.assertIn((4, 5), [(a.row_id, b.row_id) for a, b, _ in ordered])
        self.assertNotIn(5, [a.row_id for a, _, _ in ordered])

    def test_intervening_epoch_and_overlapping_error_brackets(self):
        events = self.events(17)
        events[9] = verify.Event(9, 900., 99., 5.)
        splits, count = verify.derive_partitions({'J': events}, {'J': [250.]})
        ids = [p[0].row_id for half in splits['J'] for p in half]
        self.assertEqual(count, 13)
        self.assertEqual(set(range(16)) - set(ids), {2, 8, 9})
        self.assertEqual(len(splits['J'][0]), 8)
        self.assertEqual(len(splits['J'][1]), 5)

    def test_rxte_nicer_gap_belongs_to_neither_split(self):
        epochs = [54000 + i for i in range(10)] + [58000 + i for i in range(7)]
        events = [verify.Event(i, t, 0., 1.) for i, t in enumerate(epochs)]
        partitions, count = verify.derive_partitions({verify.CONTROL: events}, {})
        train, test = partitions[verify.CONTROL]
        self.assertEqual(count, 15)
        self.assertEqual([p[0].row_id for p in train], list(range(9)))
        self.assertEqual([p[0].row_id for p in test], list(range(10, 16)))

    def test_cannot_omit_eligible_source_or_duplicate_names(self):
        for rows in ([], [{'name': 'J'}, {'name': 'J'}]):
            with self.subTest(rows=rows), self.assertRaisesRegex(ValueError, 'population'):
                verify.verify_screen(rows, {}, {'J': None})

    def test_duplicate_train_id_and_train_test_reuse(self):
        events = self.events(17)
        parts, _ = verify.derive_partitions({'J': events}, {})
        train, test = parts['J']
        def ids(points):
            return ' '.join(str(p[0].row_id) for p in points)

        for bad_train, bad_test in [(train[:-1] + train[:1], test), (train, train[:1] + test[1:])]:
            row = dict(name='J', train_start_ids=ids(bad_train), test_start_ids=ids(bad_test))
            with self.subTest(row=row), self.assertRaisesRegex(ValueError, 'split'):
                verify.verify_screen([row], {'J': events}, parts)


class ExternalTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((verify.HERE / 'parkes-intervals.json').read_text())
        # The structural checks run before fitting errors are compared.
        self.training = [(1., 200.), (3., 300.)]

    def verify(self):
        with redirect_stdout(io.StringIO()):
            verify.verify_external(self.data, self.training, 58214.)

    def test_duplicate_and_spurious_rows_are_rejected(self):
        self.data['points'].append(self.data['points'][0].copy())
        with self.assertRaisesRegex(ValueError, 'external forecast rows'):
            self.verify()

    def test_external_must_use_all_catalogue_training_before_cutoff(self):
        self.data['jbo_training_cutoff_mjd'] = 100
        with self.assertRaisesRegex(ValueError, 'cutoff'):
            self.verify()

    def test_new_series_cannot_overlap_catalogue(self):
        self.data['observing_span_mjd'][0] = 58000
        with self.assertRaisesRegex(ValueError, 'after all catalogue events'):
            self.verify()

    def test_gain_bounds_include_turning_points(self):
        self.assertEqual(verify.gain_bounds(0., 4., 1., 2., 3.), (-2., 2.))


if __name__ == '__main__':
    unittest.main()
