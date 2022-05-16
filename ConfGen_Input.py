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

    #def Input(self, systems = [[Sys1Top,Sys1Coord], [Sys2Top,Sys2Coord], [Sys3Top,Sys3Coord]]):
    def Input(self, systems = None):
        self.System1 = md.load("/home/teo/2022/ConfGen/test_site/testSystem/h2l_min.inpcrd",
                          top="/home/teo/2022/ConfGen/test_site/testSystem/h2l.prmtop")
        self.System2 = md.load("/home/teo/2022/ConfGen/test_site/testSystem/h2l_min.inpcrd",
                          top="/home/teo/2022/ConfGen/test_site/testSystem/h2l.prmtop")
        self.System3 = md.load("/home/teo/2022/ConfGen/test_site/testSystem/h2l_min.inpcrd",
                          top="/home/teo/2022/ConfGen/test_site/testSystem/h2l.prmtop")

    def Search(self):
        alignedSys1 = ConfGen.AlignSystem.alignToAxis(self.System1, NTermTop="Negative")
        alignedSys2 = ConfGen.AlignSystem.alignToAxis(self.System2, NTermTop="Positive")
        alignedSys3 = ConfGen.AlignSystem.alignToAxis(self.System3, NTermTop="Negative")

        System1_2 = ConfGen.Search.SearchSurface(fixedSystem = alignedSys1,
                                                 mobileSystem = alignedSys2,
                                                 vdWoffset = 0.5,
                                                 ndelta=0, phi=90, theta=90)
        Scores = []

        for FrameIx in range(System1_2.shape[0]):
            Scores.append(self.EvaluateScore(evalScore1, alignedSys1, System1_2[FrameIx]))
        bestSysIx = np.argmin(Scores)
        print ("{}: {}".format(np.min(Scores), bestSysIx))

        System1_2_3 = ConfGen.Search.SearchSurface(fixedSystem = System1_2[bestSysIx],
                                                   mobileSystem=alignedSys3,
                                                   vdWoffset=0.5,
                                                   ndelta=0, phi=90, theta=90)
        return (System1_2_3)

A_Obj = Main()
A_Obj.Input()
A_Obj.Search()