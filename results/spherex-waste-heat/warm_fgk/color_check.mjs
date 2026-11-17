/** Response-integrated WISE color temperature for the two warm, red published stars. */
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';

function compareColors(objects, passbands) {
  const responseZeros = [0,8.180e-15,2.415e-15,6.515e-17,5.090e-18];
  const bandIntegrals = new Map();
  const integrate = (band, f) => {
    const rows = passbands[band];
    let area = 0;
    let previous = f(rows[0]);
    for (let i = 1; i < rows.length; i++) {
      const current = f(rows[i]);
      area += (rows[i][0] - rows[i-1][0]) * (previous + current) / 2;
      previous = current;
    }
    return area;
  };
  const zero = responseZeros.map((z,b) => b ? z * integrate(b, ([wave,r])=>wave*r) : 0);
  function bandFlux(b, temperature) {
    const key = b + '/' + temperature;
    if (!bandIntegrals.has(key)) {
      const counts = integrate(b, ([wave,r]) => r / wave**4 / Math.expm1(14387.76877/(wave*temperature)));
      bandIntegrals.set(key, counts/zero[b]);
    }
    return bandIntegrals.get(key);
  }
  function temperature(blue, red, measured) {
    let lo=100, hi=1800;
    const ratio = t => bandFlux(blue,t)/bandFlux(red,t);
    if (!(measured >= ratio(lo) && measured <= ratio(hi))) throw new Error('Color outside model grid');
    for (let iteration=0; iteration<42; iteration++) {
      const mid=(lo+hi)/2;
      if (ratio(mid)<measured) lo=mid;
      else hi=mid;
    }
    return (lo+hi)/2;
  }
  const flux = mag => 10**(-0.4*mag);
  const results=[];
  for (const star of objects) {
    for (const catalog of ['allwise','catwise']) {
      for (const stellar_K of [4500,5500,6500]) {
        for (const dW1 of [-0.1,-0.05,0,0.05,0.1]) {
          for (const dW2 of [-0.1,-0.05,0,0.05,0.1]) {
            const star1=flux(Number(star.pred_w1_native)-0.004+dW1);
            const star2=flux(Number(star.pred_w2_native)-0.032+dW2);
            const star3=star1*bandFlux(3,stellar_K)/bandFlux(1,stellar_K);
            const star4=star1*bandFlux(4,stellar_K)/bandFlux(1,stellar_K);
            const exc1=flux(star[catalog+'_w1'])-star1;
            const exc2=flux(star[catalog+'_w2'])-star2;
            const obs3=flux(star.allwise_w3_or_limit), obs4=flux(star.allwise_w4_or_limit);
            const T12=temperature(1,2,exc1/exc2);
            const T23=temperature(2,3,exc2/(obs3-star3));
            const predict=(band,photosphere,T) => photosphere+exc2*bandFlux(band,T)/bandFlux(2,T);
            const published=Number(star.bb_temperature_k);
            results.push({source_id:star.source_id,catalog,stellar_K,dW1,dW2,T12_K:T12.toFixed(2),T23_K:T23.toFixed(2),
              W3_observed_over_T12: (obs3/predict(3,star3,T12)).toFixed(4),
              W4_observed_over_T12: (obs4/predict(4,star4,T12)).toFixed(4),
              W3_observed_over_published_T: (obs3/predict(3,star3,published)).toFixed(4),
              W4_observed_over_published_T: (obs4/predict(4,star4,published)).toFixed(4)});
          }
        }
      }
    }
  }
  return results;
}

const root = new URL('.',import.meta.url);
const inputs = {};
function read(name) {
  const data = readFileSync(new URL(name,root));
  inputs[name] = createHash('sha256').update(data).digest('hex');
  return data.toString('utf8');
}
function readTable(name) {
  const [header,...rows]=read(name).trim().split(/\r?\n/);
  const keys=header.split(',');
  return rows.map(line=>{
    const cells=line.split(',');
    if (cells.length!==keys.length) throw new Error('Unexpected columns: '+name);
    return Object.fromEntries(keys.map((k,i)=>[k,cells[i]]));
  });
}
const references=['../literature/fgk_baseline.csv','../literature/fgk_catwise_comparison.csv','../literature/contardo_hogg.csv'].map(readTable);
const red = readTable('../fgk_spectra/all_53.csv').filter(s=>s.continuum_check==='redder_than_both_in_both_masks' && +s.published_bb_k>=100 && +s.published_bb_k<=600);
const objects=red.map(s=>Object.assign({},...references.map(ref=>{
  const row=ref.find(r=>r.source_id===s.source_id);
  if (!row) throw new Error('Missing target '+s.source_id);
  return row;
})));
const passbands={};
for (const b of [1,2,3,4]) {
  const file=b===2?'RSR-W2-numeric.txt':'../blind/RSR-W'+b+'.txt';
  passbands[b]=read(file).trim().split(/\r?\n/).filter(l=>!l.startsWith('#')).map(l=>l.trim().split(/\s+/).map(Number));
}
const rows=compareColors(objects,passbands);
const keys=Object.keys(rows[0]);
writeFileSync(new URL('color_grid.csv',root),[keys.join(','),...rows.map(r=>keys.map(k=>r[k]).join(','))].join('\n')+'\n');
writeFileSync(new URL('color_input_hashes.json',root),JSON.stringify(inputs,null,2)+'\n');
