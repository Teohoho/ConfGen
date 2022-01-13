import mdtraj as md
import numpy as np
import datetime,sys

def evaluateScore(sys1, sys2, residues, score):
    """

    Parameters
    ----------
    sys1, sys2:     MDTrajTrajectory Object
                    MDTraj trajectory to use for the system whose fitness
                    needs computing
    residues:       list of lists of str
                    Residues to consider when computing scores. These can be
                    simple "1D" strings, if all combinations of residues are
                    wanted or nested ("2D") lists, if only pairs between
                    the elements of the 2 lists are wanted.
                    Ex: ["VAL","LEU","ILE"] means that all distances between
                    these residues will be computed and scored.
                    [["ASP", "GLU"]["LYS","ARG]] means that only distances between
                    ASP/GLU and LYS/ARG will be computed and scored.
    score:          list of ints
                    each list passed at the "residues" argument will be scored
                    by the score with the same index, from this list, so the length
                    of that argument needs to be the same as this one's


    Returns
    -------
    score:          tuple
                    Tuple containing:
                    bad charge contacts
                    good charge contacts
                    neighboring hydrophobics

    Notes
    -----
    Most of this script is taken from my old "ComputeFitness" script

    TODO: Add a check to see if the two systems either have the same number
    of frames, or if either of them has 1 frame
    """

    # First we run some checks on the input
    if not (isinstance(sys1, md.core.trajectory.Trajectory) and
            isinstance(sys2, md.core.trajectory.Trajectory)):
        raise TypeError("At least one of the inputs isn't an mdtraj trajectory object.")

    ## It's easier to work with both systems brought into one
    ## since MDTraj's distance computing library is "1000x faster
    ## than the native numpy implementation". To do this, both
    ## trajectories need to have the same number of frames

    if sys2.n_frames > sys1.n_frames:
        longerSys = sys2
        shorterSys = sys1
    else:
        longerSys = sys1
        shorterSys = sys2

    shorterSystemp = shorterSys

    for i in range(1,longerSys.n_frames):
        shorterSys = shorterSys.join(shorterSystemp)

    print (shorterSys.n_frames)
    fullSystem = longerSys.stack(shorterSys)
    print (fullSystem)

    sys.exit()

    # Generate residues matrix
    resList1 = resList2 = []
    resList = [resList1, resList2]
    sysList = [sys1,sys2]

    resIx = 0
    while resIx < 2:
        for residue in sysList[resIx].topology.residues:
            resList[resIx].append((str(residue)[:3]))
        resIx += 1

    #print (resList)
    # Compute contact score
    StartTime = datetime.datetime.now()
    # Transform the lists from residues to numpy arrays
    # so they're easier to work with
    for resIx in range(len(residues)):
        residues[resIx] = np.array(residues[resIx])

    for group in residues:
        # All residues cross
        if (group.ndim == 1):
            sele = ""
            for resname in group:
                sele = "{} {}".format(sele,resname)
            print (sele)
            sys1Indices = sys1.topology.select("resname {}".format(sele))
            print (sys1Indices)

    sys.exit()
