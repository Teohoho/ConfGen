import ConfGen
import numpy as np
import mdtraj as md
import subprocess

def hydrophobicScore(Sys1,Sys2, probe_radius=0.14):

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

    fullSystem = Sys1.stack(Sys2)
    sasa = md.shrake_rupley(fullSystem, mode="residue", probe_radius=probe_radius)[0]

    for resIx in range(fullSystem.n_residues):
        resName = fullSystem.topology.residue(resIx).code
        resAccessibility = (sasa[resIx]/maxASA[resName]) * solvFreeEn[resName]

        score = score + resAccessibility
    return (score)

class Main:
    def __init__(self, outputDir):
        self.ShortToLong = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                            "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET",
                            "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER",
                            "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
        self.outputDir = outputDir
        pass
    
    def EvaluateScore(self, Func, Sys1, Sys2, **kwargs):
        return(Func(Sys1,Sys2, **kwargs))

    def Gen3D(self, seq, TMResIDs):
        helices = []
        for TM in range(len(TMResIDs)):
            helices.append([self.ShortToLong[x] for x in seq[TMResIDs[TM][0]:TMResIDs[TM][1]+1]])
        print ("The following helices will be generated:")
        for i in range(len(helices)):
            print (helices[i])
        leapFile = "{}/tempLEaP.in".format(self.outputDir)
        with open(leapFile, 'w') as f:
            f.write("source leaprc.protein.ff19SB\n")
            for TM in range(len(helices)):
                helixSeq = ""
                for i in range(len(helices[TM])):
                    helixSeq = helixSeq + " " + helices[TM][i]
                f.write("TM{} = sequence {{{}}}\n".format(TM,helixSeq))
                f.write("impose TM{} {{ {{ 1 999 }} }} {{ {{ $N $CA  $C $N -47 }} {{ $C $N  $CA $C -57.8 }} }}\n".format(TM))
                f.write("savepdb TM{0} {1}/TM{0}.pdb\n".format(TM, self.outputDir))
                helices[TM] = "{}/TM{}.pdb".format(self.outputDir, TM)
            f.write("quit\n")
        print("The generated helices have been saved at:")
        for i in range(len(helices)):
            print (helices[i])
        tleapCommand = ['tleap', '-f', leapFile]
        subprocess.Popen(tleapCommand, stdout=subprocess.DEVNULL).communicate()

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

    def Search(self, scoreFunc, rounds=1, probe_radius=0.14, **kwargs):
        outputFN = "/home/teo/2022/ConfGen/test_site/realSys/Output"
        counter = 0
        orientation = {0:"negative", 1:"positive"}
        for SystemIx in range(len(self.systems)):
            alignedSys = ConfGen.AlignSystem.alignToAxis(self.systems[SystemIx],
                                                         NTermTop=orientation[counter % 2])
            counter = counter + 1
            self.systems[SystemIx].xyz[0] = alignedSys
            self.systems[SystemIx].save("{}/TM{}_init.pdb".format(self.outputDir, SystemIx))

        fixedSys = self.systems[0]

        tempList = [fixedSys]
        for mobileIx in range(1, len(self.systems)):
            mobileSys = self.systems[mobileIx]

            tempList2 = []
            for fixedIx in range(len(tempList[-1])):
                fixedSys = tempList[-1][fixedIx]
                fixedSys.save("/home/teo/2022/ConfGen/tempOut/fixedSys_{}_{}.pdb".format(mobileIx,fixedIx))
                currentSearch = ConfGen.Search.SearchSurface(fixedSystem=fixedSys,
                                                             mobileSystem=mobileSys, **kwargs)

                currentSearch.save("{}/fullSearch_{}_{}.dcd".format(self.outputDir, mobileIx, fixedIx))
                print ("File saved successfully: {}/fullSearch_{}_{}.dcd".format(self.outputDir, mobileIx, fixedIx))

                Scores = []
                for FrameIx in range(currentSearch.n_frames):
                    Scores.append(self.EvaluateScore(scoreFunc, fixedSys, currentSearch[FrameIx], probe_radius=probe_radius))
                print(Scores)

                ## Get top 'rounds' scores
                for i in range(rounds):
                    ix = i + rounds*fixedIx
                    score = np.argsort(Scores)[i]
                    print("Score #{}: {}".format(i+1, score))
                    tempList2.append(fixedSys.stack(currentSearch[score]))
                    centerChain = tempList2[ix].topology.select("chainid {}".format(tempList2[ix].n_chains-1))
                    tempList2[ix].save("{}/FixedSys_{}_{}_{}.pdb".format(self.outputDir, mobileIx, fixedIx, i))
                    print ("File saved successfully: {}/FixedSys_{}_{}_{}.pdb".
                           format(self.outputDir, mobileIx, fixedIx, i))
                    tempList2[ix].xyz[0] = tempList2[ix].xyz[0] - np.average(
                        tempList2[ix].xyz[0][centerChain[0]:centerChain[-1]], axis=0)
                    tempList2[ix].save(
                        "{}/FixedSys_Centered_{}_{}_{}.pdb".format(self.outputDir, mobileIx, fixedIx, i))
                    print("File saved successfully: {}/FixedSys_Centered_{}_{}_{}.pdb".
                          format(self.outputDir, mobileIx, fixedIx, i))
            tempList.append(tempList2)
            print("len(tempList): {}".format(len(tempList[-1])))


## Example, substitute your own files ##
A_Obj = Main(outputDir="/home/teo/2022/ConfGen/test_site/realSys/Output")
ZAR4 = A_Obj.Gen3D(seq="MVDAVVTVFLEKTLNILEEKGRTVSDYRKQLEDLQSELKYMQSFLKDAERQKRTNETLRTLVADLRELVYEAEDILVDCQLADGDDGNEQRSSNAWLSRLHPARVPLQYKKSK",
                   TMResIDs=[[0,17],[27,49],[55,80],[88,110]])

A_Obj.Input(ZAR4)
A_Obj.Search(hydrophobicScore, probe_radius=0.75, rounds=3, theta=90, phi=90, delta=0.27, ndelta=0, offset=1.5)
