import ConfGen
import sys
# Load the systems
Sys1Top = "../ZAR1_individual_helics/h2/h2l.prmtop"
Sys1Coord = "../ZAR1_individual_helics/h2/h2l_min.inpcrd"
Sys2Top = "../ZAR1_individual_helics/h3/h3s.prmtop"
Sys2Coord = "../ZAR1_individual_helics/h3/h3s_min.inpcrd"

Sys1 = ConfGen.TrajLoader.TrajectoryLoader(Sys1Top, Sys1Coord)
Sys2 = ConfGen.TrajLoader.TrajectoryLoader(Sys2Top, Sys2Coord)

# Align the systems to the Y axis
#alignedPos1 = ConfGen.AlignSystem.alignToAxis(Sys1[0], axis="y")
#alignedPos2 = ConfGen.AlignSystem.alignToAxis(Sys2[0], axis="y")

# Save the aligned systems and visualize them in VMD, to check
# that everything went smoothly
Aligned1FN = "../Sys1.rst7"
Aligned2FN = "../Sys2.rst7"
#ConfGen.TrajWriter.writeTraj(alignedPos1, out_FN = Aligned1FN)
#ConfGen.TrajWriter.writeTraj(alignedPos2, out_FN = Aligned2FN)

# Set minimum distance between helices
distMin = Sys1[1] + Sys2[1] + 0.1
print (distMin)
# Search
#conformations = ConfGen.Search.Search(alignedPos1, alignedPos2, "temp", distMin, ndelta=7)
# Write
Found2FN = "../sys2_moved.dcd"
#ConfGen.TrajWriter.writeTraj(conformations, out_FN=Found2FN)

# Score
Sys1 = ConfGen.TrajLoader.TrajectoryLoader(Sys1Top, Aligned1FN)[0]
Sys2 = ConfGen.TrajLoader.TrajectoryLoader(Sys2Top, Found2FN)[0]

orderedFrames = ConfGen.Scoring.evaluateScore(Sys1, Sys2,[
                                        #[["ILE"], ["ILE"], ["ILE"]]
                                        #[["LEU","ALA"], ["ASN","GLY"], ["PHE"]],
                                        [["ILE"], ["LEU"], ["VAL"], ["PHE"], ["ALA"], ["MET"], ["TYR"], ["TRP"]],
                                        [["ASP", "GLU"], ["LYS", "ARG"]]
                                        ],
                                        score=[1,10], cutoff=[1,1])

#print (orderedFrames)
outDir = "./output"
for frameIx in range(len(orderedFrames)):
    #print(frameIx)

    ConfGen.TrajWriter.writeTraj(orderedFrames[frameIx][0].xyz,
                                 out_FN="{}/rank{}.pdb".format(outDir,frameIx),
                                 topology=orderedFrames[frameIx][0].topology,
                                 b_factors=orderedFrames[frameIx][1])
