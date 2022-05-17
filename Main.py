import ConfGen
import numpy as np
import mdtraj as md

def evalScore1(Sys1,Sys2):
    return (0)

class Main:
    def __init__(self):
        pass
    
    def EvaluateScore(self, Func, Sys1, Sys2):
        return(Func(Sys1,Sys2))

    #TODO: Generate a helical structure, from a user-provided sequence.
    #TODO: Add secondary structure prediction (?)
    #TODO: Add *something* related to topology. An easier way to provide top. info.

    def Input(self, systems):
        self.systems = []
        for SystemIx in range(len(systems)):
            system = md.load(systems[SystemIx][1],
                             top=systems[SystemIx][0])
            self.systems.append(system)

    def generateTopology(self, outputFN="./"):
        wholeSys = self.systems[0]
        for SystemIx in range(1, len(self.systems)):
            currentSys = self.systems[SystemIx]
            currentSys.xyz = currentSys.xyz + 1*SystemIx
            wholeSys = wholeSys.stack(self.systems[SystemIx])
        ConfGen.TrajWriter.writeTraj(wholeSys.xyz[0], out_FN=outputFN,
                                     topology=wholeSys.topology)

    def Search(self):
        print ("len(systems): {}".format(len(self.systems)))
        alignedSystems = []
        counter = 0
        orientation = {0:"negative", 1:"positive"}
        for SystemIx in range(len(self.systems)):
            alignedSys = ConfGen.AlignSystem.alignToAxis(self.systems[SystemIx],
                                                         NTermTop=orientation[counter % 2])
            alignedSystems.append(alignedSys)

        print (len(alignedSystems))
        fixedSys = alignedSystems[0]
        print(fixedSys.shape)
        for mobileIx in range(1,len(alignedSystems)):

            currentSearch = ConfGen.Search.SearchSurface(fixedSystem = fixedSys,
                                                     mobileSystem = alignedSystems[mobileIx],
                                                     vdWoffset = 0.5,
                                                     ndelta=0, phi=90, theta=90)

            outputFN = "/home/teo/2022/ConfGen/test_site/Output"
            #ConfGen.TrajWriter.writeTraj(currentSearch, out_FN="{}/fullSystemTest.dcd".format(outputFN))

            Scores = []

            for FrameIx in range(currentSearch.shape[0]):
                Scores.append(self.EvaluateScore(evalScore1, fixedSys, currentSearch))
            bestSysIx = np.argmin(Scores)
            print ("{}: {}".format(np.min(Scores), bestSysIx))

            fixedSys = np.vstack((fixedSys,
                                  currentSearch[bestSysIx]))
            fixedSys = fixedSys - np.average(fixedSys, axis=0)

            print(fixedSys.shape)
        return (fixedSys)

## Example, substitute your own files ##
A_Obj = Main()
A_Obj.Input([["/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb", "/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb"],
             ["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             ["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"]
            ])

outputFN = "/home/teo/2022/ConfGen/test_site/Output"
A_Obj.generateTopology(outputFN="{}/fullSystemTest.pdb".format(outputFN))

finalSystem = A_Obj.Search()
ConfGen.TrajWriter.writeTraj(finalSystem, out_FN="{}/fullSystemTest.dcd".format(outputFN))