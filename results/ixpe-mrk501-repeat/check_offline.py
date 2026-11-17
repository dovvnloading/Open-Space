"""Check frozen input bytes and deposited Stokes comparison with the stdlib."""
import hashlib
import ast
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def close(x,y):
    if not math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-11):raise AssertionError(f'{x} != {y}')

def verify():
    frozen=ROOT/'frozen';packet=read(frozen/'conditional_packet.json')
    diffs=ROOT/'git-blob-differences'
    expected_names={'baseline_inputs.json','centroid_rule.py','forward_fold.py','rmf_sensitivity.py','summary.py','synthetic_stokes.py','test_development.py','test_forward_fold.py','test_prospective.py'}
    if {x.name for x in diffs.iterdir()}!=expected_names:raise AssertionError('Git blob difference set changed')
    for blob in diffs.iterdir():
        old=blob.read_bytes();executed=(frozen/blob.name).read_bytes()
        if old==executed or old.replace(b'\r\n',b'\n')!=executed.replace(b'\r\n',b'\n'):raise AssertionError(f'not a newline-only difference: {blob.name}')
        if blob.suffix=='.py' and ast.dump(ast.parse(old),include_attributes=False)!=ast.dump(ast.parse(executed),include_attributes=False):raise AssertionError('Python AST changed: '+blob.name)
        if blob.suffix=='.json' and json.loads(old)!=json.loads(executed):raise AssertionError('JSON content changed: '+blob.name)
    for name,sha in packet['repository_sha256'].items():
        actual=hashlib.sha256((frozen/name).read_bytes()).hexdigest()
        if actual!=sha:raise AssertionError(f'frozen file changed: {name}')
    receipt=ROOT/'evidence'/'receipt.json';comparison=ROOT/'evidence'/'comparison.json'
    expected={'receipt.json':'c88d7913437e5ba7809574b4ced4fdc280960179c65f92fe82ef256cfa067ea3','comparison.json':'c5d87cd4c0a49f4b712e1022104cff20e58bca7ae9de32bb4e8f499774362046'}
    for p in (receipt,comparison):
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected[p.name]:raise AssertionError(f'evidence changed: {p.name}')
    base=read(frozen/'frozen_baseline.json');result=read(comparison)
    if result['packet_sha256']!=hashlib.sha256((frozen/'conditional_packet.json').read_bytes()).hexdigest():raise AssertionError('packet digest mismatch')
    if result['receipt']!=read(receipt):raise AssertionError('receipt mismatch')
    for label in ('nominal','xy_centroid'):
        b=base[label];h=result['summary'][label];reported=result['comparison']['results'][label]
        for item in (b,h):
            sums=[0.,0.,0.]
            for du in item['du_sufficient']:
                for i in range(3):sums[i]+=du['source_sum'][i]-du['scale']*du['background_sum'][i]
            close(item['q'],sums[1]/sums[0]);close(item['u'],sums[2]/sums[0])
            c=[[sum(du['source_outer'][i][j]+du['scale']**2*du['background_outer'][i][j] for du in item['du_sufficient']) for j in range(3)] for i in range(3)]
            intensity,q,u=sums
            jac=[[-q/intensity**2,1/intensity,0],[-u/intensity**2,0,1/intensity]]
            for i in range(2):
                for j in range(2):
                    reconstructed=sum(jac[i][k]*c[k][l]*jac[j][l] for k in range(3) for l in range(3))
                    close(item['cov_qu'][i][j],reconstructed)
        diff=[h['q']-b['q'],h['u']-b['u']]
        cov=[[b['cov_qu'][i][j]+h['cov_qu'][i][j] for j in range(2)] for i in range(2)]
        for i in range(2):
            close(diff[i],reported['difference'][i])
            for j in range(2):close(cov[i][j],reported['covariance'][i][j])
        det=cov[0][0]*cov[1][1]-cov[0][1]*cov[1][0]
        if det<=0:raise AssertionError('non-positive covariance')
        chi=(cov[1][1]*diff[0]**2-(cov[0][1]+cov[1][0])*diff[0]*diff[1]+cov[0][0]*diff[1]**2)/det
        close(chi,reported['chi2']);close(math.exp(-chi/2),reported['conditional_p'])
        if reported['screening_crossing'] != (math.exp(-chi/2)<5.7e-7):raise AssertionError('threshold mismatch')
        print(f'{label}: delta=({diff[0]:.8f}, {diff[1]:.8f}), chi2={chi:.8f}, p={math.exp(-chi/2):.8f}, crossing={reported["screening_crossing"]}')
    print('Frozen hashes, nine newline-only/parsed-content comparisons, and deposited arithmetic passed; no raw event data were inspected.')

if __name__=='__main__':verify()
