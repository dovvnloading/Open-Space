import base64, xml.etree.ElementTree as ET, struct, datetime, pathlib, json
p=pathlib.Path(__file__).parent
root=ET.parse(p/'selected.xml').getroot(); ns={'v':'http://www.ivoa.net/xml/VOTable/v1.3'}
data=base64.b64decode(root.find('.//v:STREAM',ns).text); offset=0
def string():
 global offset
 n=struct.unpack_from('>i',data,offset)[0]; offset+=4
 v=data[offset:offset+n].decode(); offset+=n; return v
def number(fmt):
 global offset
 v=struct.unpack_from('>'+fmt,data,offset)[0]; offset+=struct.calcsize(fmt); return v
obj={'obsid':string(),'name':string(),'start_mjd':number('d'),'end_mjd':number('d'),'public_mjd':number('i'),'exposure_a_s':number('d'),'exposure_b_s':number('d'),'status':string()}
for k in ('start_mjd','end_mjd','public_mjd'):
 obj[k+'_utc']=(datetime.datetime(1858,11,17)+datetime.timedelta(days=obj[k])).isoformat()
print(json.dumps(obj,indent=2))
(p/'obs_metadata.json').write_text(json.dumps(obj,indent=2)+'\n')
