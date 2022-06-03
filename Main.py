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
    maxASA = {"A":129, "C":167, "D": 193, "E": 223, "F":240,
              "G":104, "H":224, "I":197, "K": 236, "L":201, "M":224,
              "N": 195, "P":159, "Q": 225, "R": 274, "S": 155,
              "T":172, "V":174, "W":285, "Y":263}
    #solvFreeEn as calculated by Kraml et al. 2019 (https://doi.org/10.1021/acs.jctc.9b00742)
    solvFreeEn = {"A": -2.6, "C": -3.4, "D": -45, "E": -40.9, "F": -2.1,
              "G": -4.2, "H": -8.2, "I": -1.2, "K": -36.5, "L": -1, "M": -3.1,
              "N": -8.9, "P": -1.9, "Q": -9.1, "R": -32.5, "S": -7.5,
              "T": -5, "V": -1.8, "W": -3.5, "Y": -4.6}

    fullSystem = Sys1.stack(Sys2)
    sasa = md.shrake_rupley(fullSystem, mode="residue")[0]

    for resIx in range(fullSystem.n_residues):
        resName = fullSystem.topology.residue(resIx).code
        resAccessibility = (sasa[resIx]/maxASA[resName]) * solvFreeEn[resName]
        #print ("ResName: {}; sasa: {}; solvFreeEn: {}; resAccessibility: {}".format(resName, sasa[resIx], solvFreeEn[resName], resAccessibility))
        score = score + resAccessibility
    return (score)

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
        leapFile = "/home/teo/2022/ConfGen/test_site/tempLEaP.in"
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
            currentSearch = ConfGen.Search.SearchSurface(fixedSystem=fixedSys,
                                                         mobileSystem=self.systems[mobileIx], vdWoffset=0.5,
                                                         theta=180, phi=180, ndelta=0)

            outputFN = "/home/teo/2022/ConfGen/test_site/Output"
            ConfGen.TrajWriter.writeTraj(currentSearch, out_FN="{}/fullSearch_{}.dcd".format(outputFN,mobileIx))

            Scores = []

            for FrameIx in range(currentSearch.shape[0]):
                self.systems[mobileIx].xyz[0] = currentSearch[FrameIx]
                Scores.append(self.EvaluateScore(scoreFunc, fixedSys, self.systems[mobileIx]))
            print(Scores)
            bestSysIx = np.argmin(Scores)
            self.systems[mobileIx].xyz[0] = currentSearch[bestSysIx]
            print ("{}: {}".format(np.min(Scores), bestSysIx))

            fixedSys = fixedSys.stack(self.systems[mobileIx])
            fixedSys.save("{}/FixedSys_{}.pdb".format(outputFN, mobileIx))
            fixedSys.xyz[0] = fixedSys.xyz[0] - np.average(fixedSys.xyz[0], axis=0)
            fixedSys.save("{}/FixedSys_Centered_{}.pdb".format(outputFN, mobileIx))


            #print(fixedSys.shape)
        return (fixedSys)

## Example, substitute your own files ##
A_Obj = Main()
systemsIn = A_Obj.Gen3D(seq="KKKKKKKKKIKKKKKKKKKGGGGKKKKKKKKKIKKKKKKKKKGGGGKKKKKKKKKIKKKKKKKKK", TMResIDs=[[0,18],[23,41],[46,64] ], tempDir="/home/teo/2022/ConfGen/test_site")

print (systemsIn)

A_Obj.Input(
            systemsIn[0:2]
             #[["/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb", "/home/teo/2022/ConfGen/test_site/testSystem/2MTS.pdb"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"],
             #["/home/teo/2022/ConfGen/test_site/testSystem/h3s.prmtop", "/home/teo/2022/ConfGen/test_site/testSystem/h3s_min.inpcrd"]]
            )

#outputFN = "/home/teo/2022/ConfGen/test_site/Output"
#A_Obj.generateTopology(outputFN="{}/fullSystemTest.pdb".format(outputFN))

finalSystem = A_Obj.Search(hydrophobicScore)
#finalSystem.save("{}/fullSystemTest.pdb".format(outputFN))
#ConfGen.TrajWriter.writeTraj(finalSystem, out_FN="{}/fullSystemTest.dcd".format(outputFN))