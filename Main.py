import ConfGen, sys, itertools
import numpy as np
import mdtraj as md
import subprocess

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
    #maxASA as determined by Tien et al. 2013 (doi:10.1371/journal.pone.0080635)
    maxASA = {"A": 1.29, "C": 1.67, "D": 1.93, "E": 2.23, "F":2.40,
              "G": 1.04, "H": 2.24, "I": 1.97, "K": 2.36, "L":2.01, "M":2.24,
              "N": 1.95, "P": 1.59, "Q": 2.25, "R": 2.74, "S":1.55,
              "T": 1.72, "V": 1.74, "W": 2.85, "Y": 2.63}
    #solvFreeEn as calculated by Kraml et al. 2019 (https://doi.org/10.1021/acs.jctc.9b00742)
    solvFreeEn = {"A": -2.6, "C": -3.4, "D": -45, "E": -40.9, "F": -2.1,
              "G": -4.2, "H": -8.2, "I": -1.2, "K": -36.5, "L": -1, "M": -3.1,
              "N": -8.9, "P": -1.9, "Q": -9.1, "R": -32.5, "S": -7.5,
              "T": -5, "V": -1.8, "W": -3.5, "Y": -4.6}

    #solvFreeEn = {"A": -2.6, "C": -3.4, "D": 0, "E": 0, "F": -2.1,
    #              "G": -4.2, "H": -8.2, "I": -1.2, "K": 0, "L": -1, "M": -3.1,
    #             "N":0, "P": -1.9, "Q": 0, "R": 0, "S": 0,
    #              "T": 0, "V": -1.8, "W": -3.5, "Y": -4.6}

    fullSystem = Sys1.stack(Sys2)
    sasa = md.shrake_rupley(fullSystem, mode="residue", probe_radius=0.75)[0]

    for resIx in range(fullSystem.n_residues):
        resName = fullSystem.topology.residue(resIx).code
        resAccessibility = (sasa[resIx]/maxASA[resName]) * solvFreeEn[resName]

        #print ("ResName: {}; sasa: {}; maxASA: {}; solvFreeEn: {}; ratio: {}; resAccessibility: {}".format(resName, sasa[resIx],
        #        maxASA[resName], solvFreeEn[resName], sasa[resIx]/maxASA[resName], resAccessibility))
        score = score + resAccessibility
    return (score)

#def hydrophobicScore(Sys1, Sys2):
#    fullSystem = Sys1.stack(Sys2)



class Main:
    def __init__(self):
        self.ShortToLong = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                            "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET",
                            "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER",
                            "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
        pass
    
    def EvaluateScore(self, Func, Sys1, Sys2):
        return(Func(Sys1,Sys2))

    #TODO: Add secondary structure prediction (?)
    #TODO: Add *something* related to topology. An easier way to provide top. info.

    def Gen3D(self, seq, TMResIDs, tempDir="./temp"):
        helices = []
        for TM in range(len(TMResIDs)):
            helices.append([self.ShortToLong[x] for x in seq[TMResIDs[TM][0]:TMResIDs[TM][1]+1]])
        print (helices)
        leapFile = "/home/teo/2022/ConfGen/test_site/realSys/tempLEaP.in"
        #leapFile = "D:/Munca/2022/ConfGen/test_site/tempLEaP.in"
        with open(leapFile, 'w') as f:
            f.write("source leaprc.protein.ff19SB\n")
            for TM in range(len(helices)):
                helixSeq = ""
                for i in range(len(helices[TM])):
                    helixSeq = helixSeq + " " + helices[TM][i]
                f.write("TM{} = sequence {{{}}}\n".format(TM,helixSeq))
                f.write("impose TM{} {{ {{ 1 999 }} }} {{ {{ $N $CA  $C $N -47 }} {{ $C $N  $CA $C -57.8 }} }}\n".format(TM))
                f.write("savepdb TM{0} {1}/TM{0}.pdb\n".format(TM, tempDir))
                helices[TM] = "{}/TM{}.pdb".format(tempDir,TM)
            f.write("quit\n")
        print (helices)
        tleapCommand = ['tleap', '-f', leapFile]
        subprocess.Popen(tleapCommand).communicate()

        return (helices)

    def Input(self, systems):
        self.systems = []
        for SystemIx in range(len(systems)):
            if (type(systems[SystemIx]) == str):
                system = md.load(systems[SystemIx])
            else:
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
        outputFN = "/home/teo/2022/ConfGen/test_site/realSys/Output"
        print ("len(systems): {}".format(len(self.systems)))
        counter = 0
        orientation = {0:"negative", 1:"positive"}
        for SystemIx in range(len(self.systems)):
            alignedSys = ConfGen.AlignSystem.alignToAxis(self.systems[SystemIx],
                                                         NTermTop=orientation[counter % 2])
            counter = counter + 1
            self.systems[SystemIx].xyz[0] = alignedSys
            self.systems[SystemIx].save("{}/TM{}_init.pdb".format(outputFN, SystemIx))

        #print (len(alignedSystems))
        fixedSys = self.systems[0]
        rounds = 3
        offset = 1.5
        #for mobileIx in range(1,len(self.systems)):
        #    print("fixedSys: {}".format(fixedSys))
        #    print("mobileSys: {}".format(self.systems[mobileIx]))
        #    currentSearch = ConfGen.Search.SearchSurface(fixedSystem=fixedSys,
        #                                                 mobileSystem=self.systems[mobileIx], offset=offset,
        #                                                 theta=90, phi=90, ndelta=0)
        #
        #    #outputFN = "D:/Munca/2022/ConfGen/newOutput"
        #    #ConfGen.TrajWriter.writeTraj(currentSearch, out_FN="{}/fullSearch_{}.dcd".format(outputFN,mobileIx))
        #    currentSearch.save_dcd("{}/fullSearch_{}.dcd".format(outputFN,mobileIx))
        #    print ("File saved successfully: {}/fullSearch_{}.dcd".format(outputFN,mobileIx))
        #
        #    Scores = []
        #
        #    for FrameIx in range(currentSearch.n_frames):
        #        #self.systems[mobileIx].xyz[0] = currentSearch[FrameIx]
        #        Scores.append(self.EvaluateScore(scoreFunc, fixedSys, currentSearch[FrameIx]))
        #    print(Scores)
        #    bestSysIx = np.argmin(Scores)
        #    self.systems[mobileIx] = currentSearch[bestSysIx]
        #    print ("{}: {}".format(np.min(Scores), bestSysIx))
        #
        #    fixedSys = fixedSys.stack(self.systems[mobileIx])
        #    fixedSys.save("{}/FixedSys_{}.pdb".format(outputFN, mobileIx))
        #    centerChain = fixedSys.topology.select("chainid {}".format(fixedSys.n_chains-1))
        #    fixedSys.xyz[0] = fixedSys.xyz[0] - np.average(fixedSys.xyz[0][centerChain[0]+1:centerChain[-1]], axis=0)
        #    fixedSys.save("{}/FixedSys_Centered_{}.pdb".format(outputFN, mobileIx))
        #    print ("File saved successfully: {}/FixedSys_Centered_{}.pdb".format(outputFN,mobileIx))

        tempList = []
        for mobileIx in range(1, len(self.systems)):
            mobileSys = self.systems[mobileIx]
        
            tempList2 = []
            #print("fixedSys: {}".format(fixedSys))
            # print("mobileSys: {}".format(self.systems[mobileIx]))
            if (len(tempList) > 0):
                print("len(tempList): {}".format(len(tempList[-1])))
                for fixedIx in range(len(tempList[-1])):
                    fixedSys = tempList[-1][fixedIx]
                    currentSearch = ConfGen.Search.SearchSurface(fixedSystem=fixedSys,
                                                                 mobileSystem=mobileSys, offset=offset,
                                                                 theta=30, phi=30, ndelta=0)
        
                    # outputFN = "D:/Munca/2022/ConfGen/newOutput"
                    # ConfGen.TrajWriter.writeTraj(currentSearch, out_FN="{}/fullSearch_{}.dcd".format(outputFN,mobileIx))
                    currentSearch.save("{}/fullSearch_{}_{}.dcd".format(outputFN, mobileIx, fixedIx))
                    print ("File saved successfully: {}/fullSearch_{}_{}.dcd".format(outputFN,mobileIx, fixedIx))
        
                    Scores = []
                    for FrameIx in range(currentSearch.n_frames):
                        # self.systems[mobileIx].xyz[0] = currentSearch[FrameIx]
                        Scores.append(self.EvaluateScore(scoreFunc, fixedSys, currentSearch[FrameIx]))
                    print(Scores)
                    ## Get top 3 (or N) scores
                    for i in range(rounds):
                        score = np.argsort(Scores)[i]
                        print("Score #{}: {}".format(i, score))
                        tempList2.append(fixedSys.stack(currentSearch[score]))
                        centerChain = tempList2[i].topology.select("chainid {}".format(tempList2[i].n_chains-1))
                        print(centerChain)
                        tempList2[i].save("{}/FixedSys_{}_{}_{}.pdb".format(outputFN, mobileIx, fixedIx, i))
                        print ("File saved successfully: {}/FixedSys_{}_{}_{}.pdb".
                               format(outputFN, mobileIx, fixedIx, i))
                        tempList2[i].xyz[0] = tempList2[i].xyz[0] - np.average(
                            tempList2[i].xyz[0][centerChain[0]:centerChain[-1]], axis=0)
                        tempList2[i].save(
                            "{}/FixedSys_Centered_{}_{}_{}.pdb".format(outputFN, mobileIx, fixedIx, i))
                        print("File saved successfully: {}/FixedSys_Centered_{}_{}_{}.pdb".
                              format(outputFN, mobileIx, fixedIx, i))
                tempList.append(tempList2)
                print("len(tempList): {}".format(len(tempList[-1])))
        
            else:
                currentSearch = ConfGen.Search.SearchSurface(fixedSystem=fixedSys,
                                                             mobileSystem=mobileSys, offset=offset,
                                                             theta=30, phi=30, ndelta=0)
                currentSearch.save("{}/fullSearch_{}.dcd".format(outputFN, mobileIx))
                Scores = []
                for FrameIx in range(currentSearch.n_frames):
                    # self.systems[mobileIx].xyz[0] = currentSearch[FrameIx]
                    Scores.append(self.EvaluateScore(scoreFunc, fixedSys, currentSearch[FrameIx]))
                print(Scores)
                for i in range(rounds):
                    print("Score #{}: {}".format(i, np.argsort(Scores)[i]))
                    tempList2.append(fixedSys.stack(currentSearch[np.argsort(Scores)[i]]))
                    centerChain = tempList2[i].topology.select("chainid {}".format(tempList2[i].n_chains-1))
                    print (centerChain)
                    tempList2[i].save("{}/FixedSys_{}_{}.pdb".format(outputFN, mobileIx, i))
                    print("File saved successfully: {}/FixedSys_{}_{}.pdb".
                          format(outputFN, mobileIx, i))
                    tempList2[i].xyz[0] = tempList2[i].xyz[0] - np.average(
                        tempList2[i].xyz[0][centerChain[0]:centerChain[-1]], axis=0)
                    print (np.average(tempList2[i].xyz[0][centerChain[0]:centerChain[-1]], axis=0))
                    tempList2[i].save("{}/FixedSys_Centered_{}_{}.pdb".format(outputFN, mobileIx, i))
                    print("File saved successfully: {}/FixedSys_Centered_{}_{}.pdb".
                          format(outputFN, mobileIx, i))
                tempList.append(tempList2)

                # bestSysIx = np.argmin(Scores)
                # self.systems[mobileIx] = currentSearch[bestSysIx]
                # print ("{}: {}".format(np.min(Scores), bestSysIx))

            # fixedSys = fixedSys.stack(self.systems[mobileIx])
            # fixedSys.save("{}/FixedSys_{}.pdb".format(outputFN, mobileIx))
            # centerChain = fixedSys.topology.select("chainid {}".format(fixedSys.n_chains-1))
            # fixedSys.xyz[0] = fixedSys.xyz[0] - np.average(fixedSys.xyz[0][centerChain[0]+1:centerChain[-1]], axis=0)
            # fixedSys.save("{}/FixedSys_Centered_{}.pdb".format(outputFN, mobileIx))

            #print(fixedSys.shape)
        return (fixedSys)

## Example, substitute your own files ##
A_Obj = Main()
#systemsIn = A_Obj.Gen3D(seq="AAAAAAAAKKAAAAAAAAAGGGGAAAAAAAAAKAAAAAAAAAGGGGAAAAAAAAAKAAAAAAAAA", TMResIDs=[[0,18],[23,41],[46,64] ], tempDir="/home/teo/2022/ConfGen/test_site")
#systemsIn = A_Obj.Gen3D(seq="AAAKKAIAAGGGGAAAKKAIAAGGGGAAAKKAIAAGGGGAAAKKAIAAGGGGAAAKKAIAAGGGGAAAKKAIAA", TMResIDs=[[0,8],[13,21],[26,34],[39,47],[52,60]], tempDir="/home/teo/2022/ConfGen/test_site")
ZAR4 = A_Obj.Gen3D(seq="MVDAVVTVFLEKTLNILEEKGRTVSDYRKQLEDLQSELKYMQSFLKDAERQKRTNETLRTLVADLRELVYEAEDILVDCQLADGDDGNEQRSSNAWLSRLHPARVPLQYKKSK",
                   TMResIDs=[[0,17],[27,49],[55,80],[88,110]], tempDir="/home/teo/2022/ConfGen/test_site/realSys/")

#print (systemsIn)

A_Obj.Input(
            ZAR4
             #[["/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb", "/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"]]
            #[["D:/Munca/2022/ConfGen/test_site/TM0.pdb", "D:/Munca/2022/ConfGen/test_site/TM0.pdb"],
            #["D:/Munca/2022/ConfGen/test_site/TM0.pdb", "D:/Munca/2022/ConfGen/test_site/TM0.pdb"],
            #["D:/Munca/2022/ConfGen/test_site/TM0.pdb", "D:/Munca/2022/ConfGen/test_site/TM0.pdb"],
            #["D:/Munca/2022/ConfGen/test_site/TM0.pdb", "D:/Munca/2022/ConfGen/test_site/TM0.pdb"],
            #["D:/Munca/2022/ConfGen/test_site/TM0.pdb", "D:/Munca/2022/ConfGen/test_site/TM0.pdb"]
            )

#outputFN = "/home/teo/2022/ConfGen/test_site/Output"
#A_Obj.generateTopology(outputFN="{}/fullSystemTest.pdb".format(outputFN))

finalSystem = A_Obj.Search(hydrophobicScore)
#finalSystem.save("{}/fullSystemTest.pdb".format(outputFN))
#ConfGen.TrajWriter.writeTraj(finalSystem, out_FN="{}/fullSystemTest.dcd".format(outputFN))