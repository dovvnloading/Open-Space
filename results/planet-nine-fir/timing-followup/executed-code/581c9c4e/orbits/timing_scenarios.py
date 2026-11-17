"""Metadata-only timing scenarios; does not read image pixels or recompute masks."""
import argparse, copy, hashlib, json
from pathlib import Path
import numpy as np
import propagate as p

SCENARIOS = {
    'historical': ('2006-12-26T21:03:07.5', 'unverified historical assumption'),
    'dec1': ('2006-12-26T15:59:18.072', 'paper-supported December scan label; catalog contributor unverified'),
    'dec2': ('2006-12-26T17:38:34.553', 'paper-supported December scan label; catalog contributor unverified'),
    'jun1': ('2006-06-26T12:28:36.556', 'counterfactual seasonal sensitivity'),
    'jun2': ('2006-06-26T14:07:41.262', 'counterfactual seasonal sensitivity'),
    'jun3': ('2006-06-26T15:46:46.012', 'counterfactual seasonal sensitivity'),
}

def scenario_inputs(inputs, name):
    result = copy.deepcopy(inputs)
    result['akari_midpoint_utc'] = SCENARIOS[name][0]
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario', required=True, choices=SCENARIOS)
    parser.add_argument('--out', required=True, type=Path, help='New exclusive directory')
    args = parser.parse_args()
    path = Path(__file__).with_name('orbit_input.json')
    raw = path.read_bytes(); original = json.loads(raw)
    if original['akari_midpoint_utc'] != SCENARIOS['historical'][0]:
        raise ValueError('historical input identity changed')
    args.out.mkdir(parents=True, exist_ok=False)
    states = p.family(scenario_inputs(original, args.scenario))
    result = {'scenario': args.scenario, 'scenario_time_utc': SCENARIOS[args.scenario][0],
              'role': SCENARIOS[args.scenario][1], 'input_sha256': hashlib.sha256(raw).hexdigest(),
              'assumptions': 'Geocentric; metadata-only; no mask, spacecraft-position or image test',
              'orbit_summary': states[-1], 'epochs': {}}
    for day in ['2016-12-30T03:58:11.008277', '2016-12-31T01:26:08.177374', '2024-12-27T07:30:00']:
        a,b = p.sky(states,p.epoch(day))
        result['epochs'][day] = {'ra_deg': [float(a.min()),float(a.max())], 'dec_deg': [float(b.min()),float(b.max())]}
    (args.out/'scenario.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__ == '__main__': main()
