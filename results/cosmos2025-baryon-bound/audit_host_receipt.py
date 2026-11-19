#!/usr/bin/env python3
"""Adversarial local review of COSMOS2025 source-route receipt, not a science review."""
import hashlib
import json
import pathlib
import sys

receipt_path = pathlib.Path(sys.argv[1])
code_path = pathlib.Path(sys.argv[2])
receipt = json.loads(receipt_path.read_text())
checks = receipt['checks']
assert len(checks) == 9
assert len({row['hostname'] for row in checks}) == 7
assert all(r['dns']['exception'] == 'gaierror' and r['dns']['errno'] == -3 for r in checks)
assert all(r['https']['status'] == 'not_attempted_dns_failed' for r in checks)
assert receipt['catalogue_bytes_received'] == 0
assert receipt['catalogue_sha256'] is None
assert receipt['catalogue_download_attempted'] is False
assert receipt['access_result'] == 'HOST_CHECKS_NO_HTTP_ROUTE'
assert hashlib.sha256(code_path.read_bytes()).hexdigest() == receipt['script_sha256']
report = {
    'review': 'LOCAL_SELF_AUDIT_ONLY_NOT_EXTERNAL_REFEREE',
    'passed_invariants': 9,
    'url_checks': len(checks),
    'unique_hosts': len({row['hostname'] for row in checks}),
    'network_access_in_environment': 'not_observed',
    'catalogue_analysis': 'NOT_RUN',
    'candidate_count': None,
    'exceedance_probability': None,
    'astrophysical_disposition': 'INCONCLUSIVE_DATA_ACCESS_FAILURE_NOT_NULL',
    'referee_disposition': 'BLOCK_DO_NOT_MERGE',
    'source_receipt_sha256': hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
    'probe_script_sha256': receipt['script_sha256'],
}
output = pathlib.Path(sys.argv[3])
output.write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')
print(json.dumps(report, sort_keys=True))
