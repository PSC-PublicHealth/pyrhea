import pyrheautils
import facilitybase


onn = pyrheautils.outputNotesName

sarahTest = int(onn.split('sarahltac_')[1].split('_')[0])
bedMult, bedMultFrac, rest = onn.split('_capped_')[1].split('_',2)
bedMult = int(bedMult) + float(bedMultFrac) / 10**(len(bedMultFrac))

facilitybase.HackBedMultiplier = bedMult

print "*******************"
print "*******************"
print "*******************"
print "*******************"
print "sarah test %s"%sarahTest
print "bedMult: %s"%bedMult
print "*******************"
print "*******************"
print "*******************"
print "*******************"

ratioTxt = """ADVO_3435_L 1 1
PRES_100_L 0.53498661 0.53498661
RML_5601_L 0.42933198 0.42933198
THC_225_L 0.97081868 0.97081868
THC_2544_L 0.33560435 0.42933198
THC_365_L 0.37286791 0.42933198
THC_4058_L 0.42933198 0.42933198
THC_6130_L 0.61589405 0.74264496
VIBR_9509_L 0.42933198 0.93920844"""

test = 1
ratios = {l.split(' ')[0]: float(l.split(' ')[sarahTest]) for l in ratioTxt.split('\n')}



facilitiesReplacementData = {}

for abbr, ratio in ratios.items():
    facilitiesReplacementData[abbr] = [
        [['scaleLengthOfStay', 'value'], ratio],
        [['scaleLengthOfStay', 'prov'], "test"]
    ]
