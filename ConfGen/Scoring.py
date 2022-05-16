import mdtraj as md
import numpy as np
import datetime, itertools
from operator import itemgetter

def evaluateScore(sys1, sys2, scoreFunction):
    """

    Parameters
    ----------
    sys1, sys2:     MDTrajTrajectory Object
                    MDTraj trajectory to use for the system whose fitness
                    needs computing

    scoreFunction:  list of ints
                    Scoring Function provided by the user, to use in ranking the
                    conformations

    Returns
    -------
    score:          list of tuple
                    List of tuples containing MDTraj Trajectory object of conformation
                    and its associated score


    Notes
    -----
    Most of this script is taken from my old "ComputeFitness" script

    TODO: Add a check to see if the two systems either have the same number
    TODO: of frames, or if either of them has 1 frame
    """

    # First we run some checks on the input
    if not (isinstance(sys1, md.core.trajectory.Trajectory) and
            isinstance(sys2, md.core.trajectory.Trajectory)):
        raise TypeError("At least one of the inputs isn't an mdtraj trajectory object.")

    ## Concatenate

    if sys2.n_frames > sys1.n_frames:
        tempSys = sys1
        for i in range(1,sys2.n_frames):
            sys1 = sys1.join(tempSys)
    else:
        tempSys = sys2
        for i in range(1,sys1.n_frames):
            sys2 = sys2.join(tempSys)

    fullSystem = sys1.stack(sys2)
    print (sys1.n_chains)





    # Compute contact score
    StartTime = datetime.datetime.now()

    # First we compute the contact matrix between the two objects
    # for all frames
    contactMatrix = np.zeros((fullSystem.n_frames, sys1.n_residues, sys2.n_residues))
    for i in range(0,sys1.n_residues):
        for j in range(0,sys2.n_residues):

            sys2ResId = j + sys1.n_residues
            #print(i, sys2ResId)
            contactCheck = md.compute_contacts(fullSystem, contacts=[[i,sys2ResId]]
                                               )[0].flatten()
            contactMatrix[:,i,j] = contactCheck

    ## DEBUG
    for frameix in range(fullSystem.n_frames):
        np.savetxt("temp_frame_{}.txt".format(frameix), contactMatrix[frameix], fmt="%5.3f")
    # We also generate two lists of resnames, so we have a way to get back indices later
    #for i in fullSystem.topology.residues:
    #    print (i)
    resList1 = resList[0:sys1.n_residues]
    resList2 = resList[sys1.n_residues:]

    currentTime = datetime.datetime.now()
    print ("Contact Matrix succesfully computed. Took {} seconds.".format((currentTime - StartTime).total_seconds()))

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
        crossedRes = list(crossedRes for crossedRes,_ in itertools.groupby(crossedRes))

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

    ## Having the frames and their respective scores, we can rank them
    orderedFrames = []
    for frameIx in range(len(framesScore)):
        orderedFrames.append([fullSystem[frameIx], framesScore[frameIx]])

    sortedFrames = sorted(orderedFrames, key=itemgetter(1), reverse=True)

    # We then return the sortedFrames
    return (sortedFrames)
