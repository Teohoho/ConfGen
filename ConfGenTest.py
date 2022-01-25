import numpy as np

import ConfGen
import sys
# Load the systems (monteverdi)
#Sys1Top = "../getConformations/ZAR1_individual_helics/h2/h2l.prmtop"
#Sys1Coord = "../getConformations/ZAR1_individual_helics/h2/h2l_min.inpcrd"
#Sys2Top = "../getConformations/ZAR1_individual_helics/h3/h3s.prmtop"
#Sys2Coord = "../getConformations/ZAR1_individual_helics/h3/h3s_min.inpcrd"

# Load the systems (Home)
Sys1Top =   "../ZAR1_individual_helics/h2/h2l.prmtop"
Sys1Coord = "../ZAR1_individual_helics/h2/h2l_min.inpcrd"
Sys2Top =   "../ZAR1_individual_helics/h3/h3s.prmtop"
Sys2Coord = "../ZAR1_individual_helics/h3/h3s_min.inpcrd"

#Sys3Coord = "output/rank0.pdb"


#Sys3 = ConfGen.TrajLoader.TrajectoryLoader(Sys3Coord, Sys3Coord)
Sys1 = ConfGen.TrajLoader.TrajectoryLoader(Sys1Top, Sys1Coord)
Sys2 = ConfGen.TrajLoader.TrajectoryLoader(Sys2Top, Sys2Coord)

outputFN = "./output/"

# Align the systems to the Y axis
alignedPos1 = ConfGen.AlignSystem.alignToAxis(Sys1, axis="y")
alignedPos2 = ConfGen.AlignSystem.alignToAxis(Sys2, axis="y")
#alignedPos3 = ConfGen.AlignSystem.alignToAxis(Sys3, axis="y",
#                            CenterAxisSele=["resid 0 to 3","resid 21 to 24"])

# Save the aligned systems and visualize them in VMD, to check
# that everything went smoothly
Aligned1FN = outputFN + "Sys1.rst7"
Aligned2FN = outputFN + "Sys2.rst7"
#Aligned3FN = outputFN + "Sys3.rst7"

#print (Aligned1FN.split(".")[-1])

ConfGen.TrajWriter.writeTraj(alignedPos1[0], out_FN=Aligned1FN)
ConfGen.TrajWriter.writeTraj(alignedPos2[0], out_FN=Aligned2FN)
#ConfGen.TrajWriter.writeTraj(alignedPos3, out_FN=Aligned3FN)

# Set minimum distance between helices
longAx =  alignedPos1[1] + 0.1
shortAx = alignedPos1[2] + 0.1

print(longAx, shortAx)

# Search
conformations = ConfGen.Search.Search(alignedPos2[0], (longAx,shortAx), ndelta=0, phi=360, theta=1)

# Write
Found2FN = outputFN + "sys2_moved.dcd"
ConfGen.TrajWriter.writeTraj(conformations, out_FN=Found2FN)

sys.exit()

# Score
Sys1 = ConfGen.TrajLoader.TrajectoryLoader(Sys1Top, Aligned1FN)[0]
Sys2 = ConfGen.TrajLoader.TrajectoryLoader(Sys2Top, Found2FN)[0]

orderedFrames = ConfGen.Scoring.evaluateScore(Sys1, Sys2,[
                                        #[["ILE"], ["ILE"], ["ILE"]]
                                        #[["LEU","ALA"], ["ASN","GLY"], ["PHE"]],
                                        [["ILE"], ["LEU"], ["VAL"], ["PHE"], ["ALA"], ["MET"], ["TYR"], ["TRP"]],
                                        [["ASP", "GLU"], ["LYS", "ARG"]],
                                        [["ASP"], ["GLU"]],
                                        [["LYS"], ["ARG"]]
                                        ],
                                        score=[1, 10, 0, 0], cutoff=[1, 1, 1, 1])

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