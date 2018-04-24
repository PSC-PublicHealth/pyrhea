import yaml
import os.path
import phacsl.utils.formats.yaml_tools as yaml_tools
import sys
cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../../sim"))
import pyrheautils


provenance = "Per email chains of April 4-6 with Bruce Lee, Sheryl Siegmund, Sarah Bartsch, Jay DePasse, Joel Welling, and Jim Leonard"


ratioTxt = """ADVO_3435_L 1 1
PRES_100_L 0.53498661 0.53498661
RML_5601_L 0.42933198 0.42933198
THC_225_L 0.97081868 0.97081868
THC_2544_L 0.33560435 0.42933198
THC_365_L 0.37286791 0.42933198
THC_4058_L 0.42933198 0.42933198
THC_6130_L 0.61589405 0.74264496
VIBR_9509_L 0.42933198 0.93920844"""

sarahTest = 1
ratios = {l.split(' ')[0]: float(l.split(' ')[sarahTest]) for l in ratioTxt.split('\n')}

facilitiesReplacementData = {}

for abbr, ratio in ratios.items():
    facilitiesReplacementData[abbr] = [
        [['scaleLengthOfStay', 'value'], ratio],
        [['scaleLengthOfStay', 'prov'], provenance]
    ]

facilityfactsDir = "../../../models/ChicagoLand/facilityfacts"
facilityfactsUpdateDir = "../../../models/ChicagoLand/facilityfacts_new"

for fac,changes in facilitiesReplacementData.items():
    siteFName = "%s.yaml"%fac
    fname = os.path.join(facilityfactsDir, siteFName)
    updateFname = os.path.join(facilityfactsUpdateDir, siteFName)

    with open(fname) as f:
        yData = yaml.load(f)

    pyrheautils.replaceYamlData(changes, yData)

    yaml_tools.save_one(updateFname, yData)
