import mdtraj as md
import numpy as np
import datetime,sys, itertools

def evaluateScore(sys1, sys2, residues, score, cutoff=0.7):
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
    cutoff:         float
                    distance threshold to consider two particles as being "in contact".
                    Given in nm.


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

    if (len(residues) != len(score)):
       raise ValueError("Number of residue lists and number of scores need to be the same. "
                        "You provided {} and {}".format(len(residues), len(score)))

    ## It's easier to work with both systems brought into one
    ## since MDTraj's distance computing library is "1000x faster
    ## than the native numpy implementation". To do this, both
    ## trajectories need to have the same number of frames

    if sys2.n_frames > sys1.n_frames:
        tempSys = sys1
        for i in range(1,sys2.n_frames):
            sys1 = sys1.join(tempSys)
    else:
        tempSys = sys2
        for i in range(1,sys1.n_frames):
            sys2 = sys2.join(tempSys)

    #print (sys1.topology.residue(24))
    #print (sys2.topology.residue(0))
    fullSystem = sys1.stack(sys2)
    #print (fullSystem)
    #print (fullSystem.topology.residue(25))

    # # Generate residues matrix
    # resList1 = resList2 = []
    # resList = [resList1, resList2]
    # sysList = [sys1,sys2]
    #
    # resIx = 0
    # while resIx < 2:
    #     for residue in sysList[resIx].topology.residues:
    #         resList[resIx].append((str(residue)[:3]))
    #     resIx += 1

    #print (resList)
    # Compute contact score
    StartTime = datetime.datetime.now()


    # Transform the lists from residues to numpy arrays
    # so they're easier to work with
    for resIx in range(len(residues)):
        residues[resIx] = np.array(residues[resIx])

    # Compute indices array to see which distances we need to compute
    ixArray = np.zeros((sys1.topology.n_atoms, sys2.topology.n_atoms))
    print (ixArray.shape)

    for group in residues:
        # All residues cross
        if (group.ndim == 1):
            sele = ""
            for resname in group:
                sele = "{} {}".format(sele,resname)
            print (sele)
            # Generate index matrix
            sys1Indices = fullSystem.topology.select("resid < {} and resname {}".format(sys1.topology.n_residues, sele))
            sys2Indices = fullSystem.topology.select("resid >= {} and resname {}".format(sys2.topology.n_residues, sele))

            distPairs = []
            #for Ix1 in range(len(sys1Indices)):
            #    for Ix2 in range(Ix1,len(sys2Indices)):
            #        distPairs.append([Ix1,Ix2])
            distPairs = np.array([c for c in itertools.product(sys1Indices, sys2Indices)])

            dists = md.compute_distances(fullSystem, distPairs)[0]

        if (group.ndim >= 2):
            sele1 = sele2 = ""
            seles = [sele1,sele2]
            for groupIx in range(len(group)):
                for resname in group[groupIx]:
                    seles[groupIx] = "{} {}".format(seles[groupIx],resname)
            print(seles)
            sys1

