SELECT a.fiberID AS anchor, a.bestObjID AS bestobjid, a.ra AS anchor_ra, a.dec AS anchor_dec,
 b.plate,b.mjd,b.fiberID AS fiber,b.ra,b.dec,b.run2d,b.class_noqso,b.z_noqso,b.zwarning_noqso,b.survey
FROM SpecObjAll AS a
JOIN SpecObjAll AS b ON a.bestObjID = b.bestObjID
WHERE a.plate=6497 AND a.mjd=56329
 AND a.run2d='v5_13_2' AND a.fiberID IN (88,96,100,139,192,210,406,424,676,748,858,963)
 AND a.bestObjID > 0
ORDER BY a.fiberID,b.mjd,b.plate,b.fiberID
