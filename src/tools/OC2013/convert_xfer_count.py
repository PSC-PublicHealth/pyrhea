import yaml
import csv
import os.path

rheaDanteDir = "/home/jim/projects/pyrhea/rhea-dante"
pyrheaDir = "../../.."

xferMatrixCSVFile = "data/OC_2013/Transferring_matrix_abbrev_9-18-2014_with-update-10-2-2014_copy_RHEA_Direct_fix_HSOU+silos+SCLE.csv"
xferMatrixCSVFile = os.path.join(rheaDanteDir, xferMatrixCSVFile)
outputFileName = "transfer_counts.yaml"
outputFile = outputFileName
# don't overwrite the current output file except with intent
#outputFile = os.path.join(pyrheaDir, 'models/OrangeCounty2013', outputFileName)

def processRow(row, xferData):
    "process a row of transfer data.  Mostly just removing the zeros and inserting it into xferData"
    abbrev = row[""]
    del(row[""])

    newRow = {}
    for k,v in row.items():
        v = int(v)
        if v > 0:
            newRow[k] = v
    xferData[abbrev] = newRow

def main():
    xferData = {}
    with open(xferMatrixCSVFile) as f:
        dr = csv.DictReader(f)
        for row in dr:
            processRow(row, xferData)

    with open(outputFile, "w") as f:
        f.write(yaml.dump(xferData))


main()
        

    
