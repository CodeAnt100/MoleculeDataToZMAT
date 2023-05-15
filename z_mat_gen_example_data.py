# BH3
atomArray = ["B", "H1", "H2", "H3"]

BHLength = 1.1900
bondLengthDict = {("B", "H1"): BHLength,
                  ("B", "H2"): BHLength,
                  ("B", "H3"): BHLength}

BHAngle = 120
bondAngleDict = {("H1", "B", "H2"): BHAngle,
                 ("H2", "B", "H3"): BHAngle,
                 ("H1", "B", "H3"): BHAngle}

# --------------------------------------------------------------------------------------------------------------------
# C2H6
atomArray = ["C1", "C2", "H1", "H2", "H3", "H4", "H5", "H6"]

CHLength = 1.0940
CCLength = 1.5351
CCHAngle = 111.17

bondLengthDict = {("C1", "C2"): CCLength}
for i in range(6):
    bondLengthDict[("C" + str(i // 3 + 1), "H" + str(i + 1))] = CHLength

bondAngleDict = {}
for i in range(3):
    bondAngleDict[("C2", "C1", "H" + str(i + 1))] = CCHAngle
    bondAngleDict[("C1", "C2", "H" + str(i + 4))] = CCHAngle

print(bondLengthDict)
print(bondAngleDict)

# --------------------------------------------------------------------------------------------------------------------
# AsH2O
atomArray = ["As", "O", "H1", "H2"]

AsHLength = 1.5134
AsOLength = 1.672011

HAsHAngle = 101.84
HAsOAngle = 106.65

bondLengthDict = {
    ("As", "O"): AsOLength,
    ("As", "H1"): AsHLength,
    ("As", "H2"): AsHLength,
}

bondAngleDict = {
    ("H1", "As", "O"): HAsOAngle,
    ("H2", "As", "O"): HAsOAngle,
    ("H1", "As", "H2"): HAsHAngle,
}

# --------------------------------------------------------------------------------------------------------------------
# NH3
atomArray = ["N", "H1", "H2", "H3"]

NHLength = 1
HNHAngle = 106.7

bondLengthDict = {}
for i in range(3):
    bondLengthDict[("N", "H" + str(i + 1))] = NHLength

bondAngleDict = {}
for i in range(3):
    bondAngleDict[("H" + str(i + 1), "N", "H" + str(i % 3 + 1))] = HNHAngle


# --------------------------------------------------------------------------------------------------------------------
# Cl2OS
atomArray = ["S", "O", "Cl1", "Cl2"]

SClLength = 2.072
SOLength = 1.434

ClSClAngle = 97.08
ClSOAngle = 108.00

bondLengthDict = {
    ("S", "O"): SOLength,
    ("S", "Cl1"): SClLength,
    ("S", "Cl2"): SClLength,
}

bondAngleDict = {
    ("Cl1", "S", "O"): ClSOAngle,
    ("Cl2", "S", "O"): ClSOAngle,
    ("Cl1", "S", "Cl2"): ClSClAngle,
}





