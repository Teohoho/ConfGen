import numpy as np

import mdtraj as md
import ConfGen
import sys

# Load the systems (monteverdi)
#Sys1Top =   "../getConformations/ZAR1_individual_helics/h2/h2l.prmtop"
#Sys1Coord = "../getConformations/ZAR1_individual_helics/h2/h2l_min.inpcrd"
#Sys2Top =   "../getConformations/ZAR1_individual_helics/h3/h3s.prmtop"
#Sys2Coord = "../getConformations/ZAR1_individual_helics/h3/h3s_min.inpcrd"
#Sys3Top =   "../getConformations/ZAR1_individual_helics/h3/h3s.prmtop"
#Sys3Coord = "../getConformations/ZAR1_individual_helics/h3/h3s_min.inpcrd"
#Sys4Top =   "../getConformations/ZAR1_individual_helics/h3/h3s.prmtop"
#Sys4Coord = "../getConformations/ZAR1_individual_helics/h3/h3s_min.inpcrd"

# Load the systems (Home)
Sys1Top =   "../ZAR1_individual_helics/h2/h2l.prmtop"
Sys1Coord = "../ZAR1_individual_helics/h2/h2l_min.inpcrd"
Sys2Top =   "../ZAR1_individual_helics/h3/h3s.prmtop"
Sys2Coord = "../ZAR1_individual_helics/h3/h3s_min.inpcrd"
Sys3Top =   "../ZAR1_individual_helics/h3/h3s.prmtop"
Sys3Coord = "../ZAR1_individual_helics/h3/h3s_min.inpcrd"
Sys4Top =   "../ZAR1_individual_helics/h3/h3s.prmtop"
Sys4Coord = "../ZAR1_individual_helics/h3/h3s_min.inpcrd"

#Sys3 = ConfGen.TrajLoader.TrajectoryLoader(Sys3Coord, Sys3Coord)
Sys1 = md.load(Sys1Coord, top=Sys1Top)
Sys2 = md.load(Sys2Coord, top=Sys2Top)
Sys3 = md.load(Sys3Coord, top=Sys3Top)
Sys4 = md.load(Sys4Coord, top=Sys4Top)

outputFN = "./Output/"

# Align the systems to the Y axis
print("Sys1")
alignedPos1 = ConfGen.AlignSystem.alignToAxis(Sys1, axis="y")
print("Sys2")
alignedPos2 = ConfGen.AlignSystem.alignToAxis(Sys2, axis="y")
alignedPos3 = ConfGen.AlignSystem.alignToAxis(Sys3, axis="y")
                          #  CenterAxisSele=["resid 0 to 3","resid 21 to 24"])

# Save the aligned systems and visualize them in VMD, to check
# that everything went smoothly
Aligned1FN = outputFN + "Sys1.rst7"
Aligned2FN = outputFN + "Sys2.rst7"
Aligned3FN = outputFN + "Sys3.rst7"
Aligned4FN = outputFN + "Sys4.rst7"


#print (Aligned1FN.split(".")[-1])

ConfGen.TrajWriter.writeTraj(alignedPos1[0], out_FN=Aligned1FN)
ConfGen.TrajWriter.writeTraj(alignedPos2[0], out_FN=Aligned2FN)
ConfGen.TrajWriter.writeTraj(alignedPos3[0], out_FN=Aligned3FN)

## First search conformations of Helix2 relative to Helix 1
# Set minimum distance between helices
longAx = alignedPos1[1]
shortAx = alignedPos1[2]
axesRatio = longAx/shortAx
print (longAx)
print (shortAx)
print()
print (axesRatio)
offset = longAx

# Search
conformations = ConfGen.Search.Search(alignedPos2[0], offset=offset,
                                      axesRatio=axesRatio, ndelta=0, phi=15, theta=15)

# Write
Found2FN = outputFN + "sys2_moved.dcd"
ConfGen.TrajWriter.writeTraj(conformations, out_FN=Found2FN)

# Score
Sys1 = md.load(Aligned1FN, top=Sys1Top)
Sys2 = md.load(Found2FN, top=Sys2Top)

orderedFrames = ConfGen.Scoring.evaluateScore(Sys1, Sys2,[
                                        #[["ILE"], ["ILE"], ["ILE"]]
                                        #[["LEU","ALA"], ["ASN","GLY"], ["PHE"]],
                                        [["ILE"], ["LEU"], ["VAL"], ["PHE"], ["ALA"], ["MET"], ["TYR"], ["TRP"]],
                                        [["ASP", "GLU"], ["LYS", "ARG"]],
                                        [["ASP"], ["GLU"]],
                                        [["LYS"], ["ARG"]]
                                        ],
                                        score=[10, 10, -10, -10], cutoff=[1, 1, 1, 1])

# Write the absolute score and corresponding frame to an ASCII file
with open("{}/scores.txt".format(outputFN), "w") as f:
    f.write("RANK\t|SCORE\n")
    for frameIx in range(len(orderedFrames)):
        f.write ("{}\t|{}\n".format(frameIx, orderedFrames[frameIx][1]))
        ConfGen.TrajWriter.writeTraj(orderedFrames[frameIx][0].xyz,
                                     out_FN="{}/rank{}.pdb".format(outputFN,frameIx),
                                     topology=orderedFrames[frameIx][0].topology,
                                     verbosity=False
                                     )
        #ConfGen.TrajWriter.PDBHelper("{}/rank{}.pdb".format(outputFN,frameIx))

# We generate a .pml file to color the above configurations in PyMOL, for easier visualization
# Assign a color value to each frame, based on ranking
NFrames = np.arange(0,len(orderedFrames), dtype="int")
NFrames = np.interp(NFrames, (NFrames.min(), NFrames.max()), (-255,255)).astype("int")

with open(outputFN + "color.pml", "w") as f:
    for frameIx in range(len(orderedFrames)):
        f.write("load {}/rank{}.pdb\n".format(outputFN, frameIx))
        if (NFrames[frameIx] <= 0):
            if ((NFrames[frameIx] + 255) < 16):
                f.write("color 0xff0{0}0{0}, m. rank{1}\n".format(hex((NFrames[frameIx]) + 255)[2:], frameIx))
            else:
                f.write("color 0xff{0}{0}, m. rank{1}\n".format(hex((NFrames[frameIx]) + 255)[2:], frameIx))
        elif (NFrames[frameIx] > 0):
            if ((255 - NFrames[frameIx]) < 16):
                f.write("color 0x0{0}0{0}ff, m. rank{1}\n".format(hex((255 - NFrames[frameIx]))[2:], frameIx))
            else:
                f.write("color 0x{0}{0}ff, m. rank{1}\n".format(hex(255 - NFrames[frameIx])[2:], frameIx))

print ("PyMOL input file written!")

## Then search conformations of Helix3 relative to the highest scoring conformation
## of Helix 1/2
Sys12 = md.load("./Output/rank0.pdb")
#print("Sys3")
#alignedPos3 = ConfGen.AlignSystem.alignToAxis(Sys3, axis="y")
print("Sys12")
alignedPos12 = ConfGen.AlignSystem.alignToAxis(Sys12, axis="y",
                            CenterAxisSele=["(resid 0 to 3) and chainid 0",
                                            "(resid 21 to 24) and chainid 0"])

# Write the newly aligned Sys12 and Sys3 to disk
Aligned12FN = outputFN + "Sys12.pdb"
ConfGen.TrajWriter.writeTraj(alignedPos12[0], out_FN=Aligned12FN,
                             topology=Sys12.topology)
Sys12 = md.load(Aligned12FN)

longAx = alignedPos12[1]
shortAx = alignedPos12[2]
axesRatio = longAx/shortAx
print (longAx)
print (shortAx)
print()
print (axesRatio)
offset = longAx

# Search
conformations = ConfGen.Search.Search(alignedPos3[0], offset=offset, axesRatio=axesRatio,
                                      ndelta=0, phi=15, theta=15)

# Write
Found3FN = outputFN + "sys3_moved.dcd"
ConfGen.TrajWriter.writeTraj(conformations, out_FN=Found3FN)
Sys3 = md.load(Found3FN, top=Sys3Top)

orderedFrames = ConfGen.Scoring.evaluateScore(Sys12, Sys3,[
                                        #[["ILE"], ["ILE"], ["ILE"]]
                                        #[["LEU","ALA"], ["ASN","GLY"], ["PHE"]],
                                        [["ILE"], ["LEU"], ["VAL"], ["PHE"], ["ALA"], ["MET"], ["TYR"], ["TRP"]],
                                        [["ASP", "GLU"], ["LYS", "ARG"]],
                                        [["ASP"], ["GLU"]],
                                        [["LYS"], ["ARG"]]
                                        ],
                                        score=[1, 10, 0, 0], cutoff=[1, 1, 1, 1])

# Write the absolute score and corresponding frame to an ASCII file
with open("{}/scores_2ndRun.txt".format(outputFN), "w") as f:
    f.write("RANK\t|SCORE\n")
    for frameIx in range(len(orderedFrames)):
        f.write ("{}\t|{}\n".format(frameIx, orderedFrames[frameIx][1]))
        ConfGen.TrajWriter.writeTraj(orderedFrames[frameIx][0].xyz,
                                     out_FN="{}/rank{}_2ndRun.pdb".format(outputFN,frameIx),
                                     topology=orderedFrames[frameIx][0].topology,
                                     verbosity=False
                                     )
        #ConfGen.TrajWriter.PDBHelper("{}/rank{}_2ndRun.pdb".format(outputFN,frameIx))
## Then search conformations of Helix4 relative to the highest scoring conformation
## of Helix 1/2/3
Sys123 = md.load("./Output/rank0_2ndRun.pdb")
print ("Sys4")
alignedPos4 = ConfGen.AlignSystem.alignToAxis(Sys4, axis="y")
print ("Sys123")
alignedPos123 = ConfGen.AlignSystem.alignToAxis(Sys123, axis="y",
                            CenterAxisSele=["(resid 0 to 3) and chainid 0",
                                            "(resid 21 to 24) and chainid 0"])

# Write the newly aligned Sys12 and Sys3 to disk
Aligned123FN = outputFN + "Sys123.pdb"
ConfGen.TrajWriter.writeTraj(alignedPos123[0], out_FN=Aligned123FN,
                             topology=Sys123.topology)
Sys123 = md.load(Aligned123FN)

longAx = alignedPos123[1]
shortAx = alignedPos123[2]
axesRatio = longAx/shortAx
print (longAx)
print (shortAx)
print()
print (axesRatio)
offset = longAx

# Search
conformations = ConfGen.Search.Search(alignedPos4[0], offset=offset, axesRatio=axesRatio,
                                      ndelta=0, phi=15, theta=15)

# Write
Found4FN = outputFN + "sys4_moved.dcd"
ConfGen.TrajWriter.writeTraj(conformations, out_FN=Found4FN)
Sys4 = md.load(Found4FN, top=Sys4Top)

orderedFrames = ConfGen.Scoring.evaluateScore(Sys123, Sys4,[
                                        #[["ILE"], ["ILE"], ["ILE"]]
                                        #[["LEU","ALA"], ["ASN","GLY"], ["PHE"]],
                                        [["ILE"], ["LEU"], ["VAL"], ["PHE"], ["ALA"], ["MET"], ["TYR"], ["TRP"]],
                                        [["ASP", "GLU"], ["LYS", "ARG"]],
                                        [["ASP"], ["GLU"]],
                                        [["LYS"], ["ARG"]]
                                        ],
                                        score=[1, 10, 0, 0], cutoff=[1, 1, 1, 1])

# Write the absolute score and corresponding frame to an ASCII file
with open("{}/scores_3rdRun.txt".format(outputFN), "w") as f:
    f.write("RANK\t|SCORE\n")
    for frameIx in range(len(orderedFrames)):
        f.write ("{}\t|{}\n".format(frameIx, orderedFrames[frameIx][1]))
        ConfGen.TrajWriter.writeTraj(orderedFrames[frameIx][0].xyz,
                                     out_FN="{}/rank{}_3rdRun.pdb".format(outputFN,frameIx),
                                     topology=orderedFrames[frameIx][0].topology,
                                     verbosity=False
                                     )
        #ConfGen.TrajWriter.PDBHelper("{}/rank{}_3rdRun.pdb".format(outputFN,frameIx))
