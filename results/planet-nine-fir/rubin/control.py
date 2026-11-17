from pathlib import Path
import sys
from urllib.parse import urlencode
from query import TAP,run
q="SELECT TOP 1 o.oid,o.n_det,o.meanra,o.meandec,o.firstmjd,o.lastmjd FROM alerce_tap.object AS o WHERE o.sid=1 AND o.oid=170019696403809402"
run('positive_control',TAP+urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','MAXREC':5,'QUERY':q}),Path(sys.argv[1]))
