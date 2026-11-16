"""The number reader: what counts as a stated quantity, and with what precision."""

import unittest

from open_space.verify import numbers


def read(text):
    return numbers.extract(text)


def values(text, kind="value"):
    return [number.value for number in read(text) if number.kind == kind]


class PrecisionTest(unittest.TestCase):
    def test_decimal_places_set_the_claimed_interval(self):
        (number,) = read("The advance is 12.99 seconds.")
        self.assertEqual(number.value, 12.99)
        self.assertAlmostEqual(number.ulp, 0.01)
        self.assertEqual(number.sig, 4)
        self.assertTrue(number.precise)
        low, high = number.interval()
        self.assertAlmostEqual(low, 12.985)
        self.assertAlmostEqual(high, 12.995)

    def test_trailing_zero_after_the_point_is_significant(self):
        (number,) = read("r = 0.280 in that band")
        self.assertEqual(number.sig, 3)
        self.assertAlmostEqual(number.ulp, 0.001)

    def test_leading_zeros_are_not_significant(self):
        (number,) = read("a tail fraction of 0.00093")
        self.assertEqual(number.sig, 2)
        self.assertFalse(number.precise)

    def test_thousands_separators_and_exact_integers(self):
        (number,) = read("It passes 49,751 km away.")
        self.assertEqual(number.value, 49751)
        self.assertEqual(number.ulp, 1)
        self.assertEqual(number.sig, 5)

    def test_round_numbers_keep_only_their_certain_digits(self):
        (number,) = read("in about 205,000 years")
        self.assertEqual(number.sig, 3)
        self.assertTrue(number.approximate)
        self.assertIn((205000.0, 1000.0, "rounded to its trailing zeros"), number.readings)

    def test_rounding_rule(self):
        self.assertTrue(numbers.round_matches(12.9937949, 12.99, 0.005))
        self.assertTrue(numbers.round_matches(12.985, 12.99, 0.005))
        self.assertFalse(numbers.round_matches(12.9961, 12.99, 0.005))


class NotationTest(unittest.TestCase):
    def test_signs_and_ranges(self):
        found = read("offsets +4.97 and −14.70, range 10.33–15.77 and 0.39-0.99")
        self.assertEqual([number.value for number in found], [4.97, -14.70, 10.33, 15.77, 0.39, 0.99])
        self.assertTrue(found[0].signed and found[1].signed)
        self.assertFalse(found[3].signed)

    def test_scientific_notation_in_three_spellings(self):
        self.assertEqual(values("rate −8.19×10⁻¹² s/s"), [-8.19e-12])
        self.assertEqual(values(r"$1.3271244\times10^{20}$ m"), [1.3271244e20])
        self.assertEqual(values("about 2.5e-3 of it"), [2.5e-3])

    def test_exponent_shared_by_a_list(self):
        first, second = read("stayed between 3.85 and 4.49 × 10⁻⁸ W")
        self.assertEqual(first.value, 3.85)
        self.assertAlmostEqual(first.readings[0][0], 3.85e-8)
        self.assertAlmostEqual(second.value, 4.49e-8)

    def test_exponent_applied_to_a_bracket(self):
        found = read(r"measured $(-9.8\pm2.8)\times10^{-12}$ s/s")
        self.assertEqual([number.value for number in found], [-9.8, 2.8])
        self.assertAlmostEqual(found[0].readings[0][0], -9.8e-12)
        self.assertAlmostEqual(found[1].readings[0][0], 2.8e-12)

    def test_percent_and_word_multipliers(self):
        (percent,) = read("it removes about 98% of the light")
        self.assertTrue(percent.percent)
        self.assertIn((0.98, 0.01, "percent as a fraction"), percent.readings)
        (million,) = read("spans 1.13 million parsecs")
        self.assertAlmostEqual(million.readings[0][0], 1.13e6)

    def test_table_header_scale_applies_to_its_column(self):
        table = "| Target | Rate (10⁻¹² s/s) | Share (%) |\n|---|---:|---:|\n| A | 8.19 | 42.5 |\n"
        rate, share = [number for number in read(table) if number.line == 3]
        self.assertAlmostEqual(rate.readings[0][0], 8.19e-12)
        self.assertAlmostEqual(share.readings[0][0], 0.425)


class NotQuantitiesTest(unittest.TestCase):
    def test_catalogue_designations_are_identifiers_or_skipped(self):
        text = ("TIC 384244835, KOI-7865.01, HD 34445, PSR J1813−2113, 4U 1538−522, GX 1+4, "
                "(41013) 1999 UH22, 2I/Borisov, CD−27 2812, ZTF J0720+6439, GRB 211211A")
        self.assertEqual(values(text), [])
        self.assertIn(384244835.0, values(text, "identifier"))
        self.assertIn(34445.0, values(text, "identifier"))

    def test_dates_years_and_references(self):
        text = "On 24 Dec 2026 and 8, 9 and 11 Oct 2027; in 2069–2071; see Table 3, Equation 4 and Section 5.2."
        self.assertEqual(values(text), [])
        self.assertEqual(values(text, "year"), [2069.0, 2071.0])
        self.assertEqual(len(values(text, "date")), 6)
        self.assertEqual(values(text, "reference"), [3.0, 4.0, 5.2])

    def test_code_links_digests_and_timestamps_are_ignored(self):
        text = ("Run `python x.py --n 12345` then read [the paper (2023)](https://x.org/2505.15580v1). "
                "SHA-256 f9d42c0bf385063181769d3614ce39cf7eb5890331ae6dba1197e82f8219b492 "
                "at 2021-12-11T13:09:59.651 and 13:10:30.351 UTC, with 87 detections.\n"
                "```sh\npython y.py --seed 999\n```\n")
        self.assertEqual(values(text), [87.0])

    def test_time_scales_before_a_number_do_not_make_it_a_name(self):
        self.assertEqual(values("after MJD 59090 on plate 6497"), [59090.0, 6497.0])

    def test_latex_commands_do_not_glue_to_numbers(self):
        self.assertEqual(values(r"$\sqrt{1.5}\,0.5=0.612$ and $3.6\pm0.1$"), [1.5, 0.5, 0.612, 3.6, 0.1])
        self.assertEqual(values(r"$(2\pi)^{5/3}=21.3944$"), [2.0, 21.3944])


class MachineSpellingTest(unittest.TestCase):
    def test_decimal_places_of_machine_numbers(self):
        self.assertEqual(numbers.decimals_of("12.990"), 3)
        self.assertEqual(numbers.decimals_of("-7"), 0)
        self.assertEqual(numbers.decimals_of("1.5e-3"), 4)
        self.assertIsNone(numbers.decimals_of("abc"))


if __name__ == "__main__":
    unittest.main()
