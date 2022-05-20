import ConfGen, sys, itertools
import numpy as np
import mdtraj as md

def knobsIntoHolesScore(Sys1,Sys2):
    pass
    ## TO BE FINISHED
#    fullSystem = Sys1.stack(Sys2)
#    contacts = md.compute_contacts(fullSystem, scheme="closest-ca")
#
#    for aaIx in sys: #KNOB
#        contactsIX = sort(list(contacts[aaIx]))[0:5] ##Get first 4
#        if (contactsIX[0]) < threshold:
#            considera knob-into-hole
#        else:
#            nu e knob-into-hole

def hydrophobicScore(Sys1,Sys2):

    score = 0
    hydrophobicIndex = {"A":0.31, "C":1.54, "D": -0.77, "E": -0.64, "F":1.79,
                        "G":0, "H":0.13, "I":1.8, "K": -0.99, "L":1.7, "M":1.23,
                        "N": -0.6, "P":0.72, "Q": -0.22, "R": -1.01, "S": -0.04,
                        "T":0.26, "V":1.22, "W":2.25, "Y":0.96}

    fullSystem = Sys1.stack(Sys2)
    fixedSysChains = " "
    for i in range(fullSystem.n_chains-1):
        fixedSysChains = "{}{} ".format(fixedSysChains, i)
    helix1CAs = fullSystem.topology.select("name CA and chainid {}".format(fixedSysChains))
    helix2CAs = fullSystem.topology.select("name CA and chainid {}".format(fullSystem.n_chains-1))
    pairs = list(itertools.product(helix1CAs,helix2CAs))
    distances = md.compute_distances(fullSystem,pairs,periodic=False)[0]
    print (distances[0])

    for pairIx in range(len(pairs)):
        res1 = fullSystem.topology.atom(pairs[pairIx][0]).residue.code
        res2 = fullSystem.topology.atom(pairs[pairIx][1]).residue.code
        pairScore = (hydrophobicIndex[res1] + hydrophobicIndex[res2])/distances[pairIx]**2
        score = score + pairScore

    return (score)

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

    def Search(self, scoreFunc):
        print ("len(systems): {}".format(len(self.systems)))
        counter = 0
        orientation = {0:"negative", 1:"positive"}
        for SystemIx in range(len(self.systems)):
            alignedSys = ConfGen.AlignSystem.alignToAxis(self.systems[SystemIx],
                                                         NTermTop=orientation[counter % 2])
            self.systems[SystemIx].xyz[0] = alignedSys

        #print (len(alignedSystems))
        fixedSys = self.systems[0]
        for mobileIx in range(1,len(self.systems)):
            print("fixedSys: {}".format(fixedSys))
            print("mobileSys: {}".format(self.systems[mobileIx]))
            currentSearch = ConfGen.Search.SearchSurface(fixedSystem = fixedSys.xyz[0],
                                                     mobileSystem = self.systems[mobileIx].xyz[0],
                                                     vdWoffset = 0.5,
                                                     ndelta=0, phi=90, theta=90)

            outputFN = "/home/teo/2022/ConfGen/test_site/Output"
            ConfGen.TrajWriter.writeTraj(currentSearch, out_FN="{}/fullSystemTest.dcd".format(outputFN))

            Scores = []

            for FrameIx in range(currentSearch.shape[0]):
                self.systems[mobileIx].xyz[0] = currentSearch[FrameIx]
                Scores.append(self.EvaluateScore(scoreFunc, fixedSys, self.systems[mobileIx]))
            print(Scores)
            bestSysIx = np.argmax(Scores)
            self.systems[mobileIx].xyz[0] = currentSearch[bestSysIx]
            print ("{}: {}".format(np.min(Scores), bestSysIx))

            fixedSys = fixedSys.stack(self.systems[mobileIx])
            fixedSys.xyz[0] = fixedSys.xyz[0] - np.average(fixedSys.xyz[0], axis=0)

            #print(fixedSys.shape)
        return (fixedSys)

## Example, substitute your own files ##
A_Obj = Main()
A_Obj.Input([["/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb", "/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb"],
             ["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             ["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"]
            ])

outputFN = "/home/teo/2022/ConfGen/test_site/Output"
A_Obj.generateTopology(outputFN="{}/fullSystemTest.pdb".format(outputFN))

finalSystem = A_Obj.Search(hydrophobicScore)
finalSystem.save("{}/fullSystemTest.pdb".format(outputFN))
#ConfGen.TrajWriter.writeTraj(finalSystem, out_FN="{}/fullSystemTest.dcd".format(outputFN))