import ConfGen
import sys

# Load the systems
Sys1 = ConfGen.TrajLoader.TrajectoryLoader("ZAR1_individual_helics/h2/h2l.prmtop", "ZAR1_individual_helics/h2/h2l_min.inpcrd")
Sys2 = ConfGen.TrajLoader.TrajectoryLoader("ZAR1_individual_helics/h3/h3s.prmtop", "ZAR1_individual_helics/h3/h3s_min.inpcrd")

# Align the systems to the Y axis
alignedPos1 = ConfGen.AlignSystem.alignToAxis(Sys1[0], axis="y")
alignedPos2 = ConfGen.AlignSystem.alignToAxis(Sys2[0], axis="y")

# Save the aligned systems and visualize them in VMD, to check
# that everything went smoothly
#ConfGen.RSTWriter.writeTraj(alignedPos1, out_FN="Sys1.rst7")
#ConfGen.RSTWriter.writeTraj(alignedPos2, out_FN="Sys2.rst7")

distMin = Sys1[1] + Sys2[1] + 1
print (distMin)

ConfGen.Search.Search(alignedPos1, alignedPos2, outRoot, distMin)
#ConfGen.InitializeSystem.InitializeSys(Sys1, Sys2)

