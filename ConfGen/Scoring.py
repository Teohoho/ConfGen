import mdtraj as md
import numpy as np
import datetime, sys, itertools
from operator import itemgetter

def evaluateScore(sys1, sys2, residues, score, cutoff=0.7):
    """

    Parameters
    ----------
    sys1, sys2:     MDTrajTrajectory Object
                    MDTraj trajectory to use for the system whose fitness
                    needs computing
    residues:       list of lists of str
                    Residues to consider when computing scores. Needs to have
                    2 "dimensions" (nested list)
                    [["ASP", "GLU"]["LYS","ARG]] means that distances between
                    ASP/GLU and LYS/ARG will be computed and scored.
                    [["ASP", "GLU"],["LYS","ARG],["VAL",ASP"]] means that
                    the following contacts will be considered:
                    ASP/LYS;ASP/ARG;ASP/VAL;ASP/ASP;
                    GLU/LYS;GLU/ARG;GLU/VAL;GLU/ASP



    score:          list of ints
                    each list passed at the "residues" argument will be scored
                    by the score with the same index, from this list, so the length
                    of that argument needs to be the same as this one's
    cutoff:         list of floats
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

    if (len(residues) != len(score) != len(cutoff)):
       raise ValueError("Number of residue lists, scores and cutoffs need to be the same. "
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

    # We also generate two lists of resnames, so we have a way to get back indices later
    resList1 = resList[0:sys1.n_residues]
    resList2 = resList[sys1.n_residues:]

    currentTime = datetime.datetime.now()
    print ("Contact Matrix succesfully computed. Took {} seconds.".format((currentTime - StartTime).total_seconds()))

    # Afterwards we can parse this array by any residues, cutoffs and assign any scores we wish.

    # To be able to rank the frames of a simulation we have to iterate through them
    framesScore = np.zeros(fullSystem.n_frames)

    # For each group in residues, we generate the pairs of residues whose contact we take into account
    for GroupIx in range(len(residues)):
        crossedRes = []
        for resIx in range(len(residues[GroupIx])):
            poppedRes = [x for i, x in enumerate(residues[GroupIx]) if i != resIx]
            allPoppedResList = []
            for i in range(len(poppedRes)):
                for j in poppedRes[i]:
                    allPoppedResList.append(j)
            for IX in range(len(residues[GroupIx][resIx])):
                for JX in range(len(allPoppedResList)):
                    crossedRes.append([residues[GroupIx][resIx][IX], allPoppedResList[JX]])
        # To account for the possibility of the same residue multiple times, we remove duplicate elements
        crossedRes.sort()
        crossedRes = list (crossedRes for crossedRes,_ in itertools.groupby(crossedRes))
        #print (crossedRes)

    # For each pair of residues computed above, we get the pair of indices that correspond
    # to that pair of aminoacids, in the contact Matrix
        for frameIx in range(framesScore.shape[0]):
            groupScore = 0 
            for crossIx in crossedRes:
                Index1 = []
                Index2 = []
                for i,element in enumerate(resList1):
                    if (element == crossIx[0]):
                        Index1.append(i)
                for i,element in enumerate(resList2):
                    if (element == crossIx[0]):
                        Index2.append(i)

                # Generate Matrix coordinate pairs:
                for dim1 in Index1:
                    for dim2 in Index2:
                        if (contactMatrix[frameIx][dim1][dim2] < cutoff[GroupIx]):
                            groupScore = groupScore + score[GroupIx]
            framesScore[frameIx] = framesScore[frameIx] + groupScore
    #print(framesScore)

    ## Having the frames and their respective scores, we can rank them
    ## Generate an array
    orderedFrames = []
    for frameIx in range(len(framesScore)):
        orderedFrames.append([fullSystem[frameIx],framesScore[frameIx]])

    sortedFrames = sorted(orderedFrames, key=itemgetter(1), reverse=True)

    # We then return the sortedFrames
    return (sortedFrames)