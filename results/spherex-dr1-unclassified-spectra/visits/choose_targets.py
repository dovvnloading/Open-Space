"""Rebuild the native-exposure followup list from the fixed cube rankings."""
import csv
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
for field,folder,ranksdir in [('nep',ROOT,'nep_two_band'),('wide',ROOT/'wide','full')]:
    def rows(path):
        with path.open() as f:
            return list(csv.DictReader(f))
    spectra={(r['source_id'],int(r['channel'])):r for r in rows(folder/'spectra.csv')}
    ranks=rows(folder/ranksdir/'ranks.csv')
    adj=rows(folder/ranksdir/'peer_channel_check.csv')
    picks={}
    def add(name,ch,reason):
        ch=int(ch)
        picks.setdefault((name,ch),[]).append(reason)
    for row in ranks:
        if (field=='wide' and int(row['template_rank'])<=5) or (field=='nep' and float(row['score'])>=5):
            add(row['source_id'],row['peak_channel'],'template')
    for row in ranks:
        if int(row['peer_rank'])==1:
            add(row['source_id'],row['peak_channel'],'shape')
    for row in sorted(adj,key=lambda r:int(r['adjusted_rank'])):
        if abs(float(row['adjusted_adjusted']))>=5:
            add(row['source_id'],row['channel_adjusted'],'field_adjusted')
    with (HERE/field/'targets.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['source_id','detector','wave_um','scale_um','channel','selection','label'])
        used=set()
        for (name,ch),reasons in picks.items():
            r=spectra[name,ch]
            width=round(float(r['wave_hi_um'])-float(r['wave_lo_um']),6)
            writer.writerow([name,r['detector'],r['wave_um'],width,ch,'+'.join(reasons),f'ch{ch}' if name in used else ''])
            used.add(name)
