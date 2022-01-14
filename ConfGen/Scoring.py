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

    fullSystem = sys1.stack(sys2)

    # # Generate residues list
    resList = []
    for residue in fullSystem.topology.residues:
        resList.append((str(residue)[:3]))

    print (resList)
    resList1 = resList[0:sys1.n_residues]
    resList2 = resList[sys1.n_residues:]

    
    #print (resList)
    # Compute contact score
    StartTime = datetime.datetime.now()

    # First we compute the contact matrix between the two objects
    # for all frames
    contactMatrix = np.zeros((fullSystem.n_frames, sys1.n_residues, sys2.n_residues))
    #print (contactMatrix.shape)
    for i in range(0,sys1.n_residues):
        for j in range(0,sys2.n_residues):
            sys2ResId = j + sys1.n_residues
            #queryAtoms = fullSystem.topology.select("resid {}".format(i))
            #haystackAtoms = fullSystem.topology.select("resid {}".format(j))
            contactCheck = md.compute_contacts(fullSystem, contacts=[[i,sys2ResId]]
                                               #query_indices=queryAtoms,
                                               #haystack_indices=haystackAtoms
                                               )[0].flatten()
            #print (contactCheck.shape)
            contactMatrix[:,i,j] = contactCheck

    # Afterwards we can parse this array by any residues, cutoffs and scores we wish.

    # Transform the lists from residues to numpy arrays
    # so they're easier to work with
    for resIx in range(len(residues)):
        residues[resIx] = np.array(residues[resIx])

    # # Compute indices array to see which distances we need to compute
    # ixArray = np.zeros((sys1.topology.n_atoms, sys2.topology.n_atoms))
    # print (ixArray.shape)
    # print (sys1.topology.n_atoms)

    # TODO: REALLY GOTTA DOCUMENT THIS PART CAUSE IT'S UNINTELLIGIBLE RIGHT NOW JESUS CHRIST
    for group in residues:
        # All residues cross
        if (group.ndim == 1):
            sele = ""
            for i in range(fullSystem.n_residues):
                for resname in group:
                    indices = fullSystem.topology.select(sele)
                    #print (fullSystem.topology.select(sele))
                    if (len(indices) > 0):
                        ixList.append(indices)


            for indicesIx in range(len(ixList)):
                haystack = [x for i,x in enumerate(ixList) if i!=indicesIx]
                haystacklist = []
                for i in range(len(haystack)):
                    e = [x for x in haystack[i] ]
                    for j in e: #Only append if it's on the other object
                        if (j > sys1.topology.n_atoms):
                            haystacklist.append(j)

                query = ixList[indicesIx]
                #print (query)
                #print (haystacklist)
                neigbors = md.compute_neighbors(fullSystem, cutoff=cutoff,
                                                query_indices=ixList[indicesIx],
                                                haystack_indices=haystacklist)
                print (neigbors)
                sys.exit()

            # Now just compute
            sysIndices = fullSystem.topology.select("resname {}".format(sele))
            # We split these indices by system
            sys1Ix = []
            sys2Ix = []
            for i in range(len(sysIndices)):
                if (sysIndices[i] > sys1.topology.n_atoms):
                    sys1Ix.append(sysIndices[i] - sys1.topology.n_atoms)
                else:
                    sys2Ix.append(sysIndices[i])

            # # Mark the relevant pairs in the ixArray
            # for i in (sys1Ix):
            #     for j in (sys2Ix):
            #         ixArray[i][j] = 1

            print (ixArray)
            # Having the index matrix, we can now generate the appropriate
            # distance pairs
            print (sysIndices.shape)
            sys.exit()

            if sysIndices[i][j] == 1:
                distPairs.append(i, j+sysIndices.shape[0])

            print (distPairs)

            distPairs = []
            #for Ix1 in range(len(sys1Indices)):
            #    for Ix2 in range(Ix1,len(sys2Indices)):
            #        distPairs.append([Ix1,Ix2])
            #distPairs = np.array([c for c in itertools.product(sys1Indices, sys2Indices)])

            #dists = md.compute_distances(fullSystem, distPairs)[0]

        if (group.ndim >= 2):
            sele1 = sele2 = ""
            seles = [sele1,sele2]
            for groupIx in range(len(group)):
                for resname in group[groupIx]:
                    seles[groupIx] = "{} {}".format(seles[groupIx],resname)
            print(seles)
            sys1

