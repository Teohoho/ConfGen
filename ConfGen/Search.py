import mdtraj
import numpy as np
from scipy.spatial.transform import Rotation as R
import warnings, sys

def centroidPoints(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])

def circularOrbit(theta):
    CTM = np.array([[np.cos(theta), 0, (-1)*np.sin(theta)],
                    [0, 1, 0],
                    [np.sin(theta),0,np.cos(theta)]])

    return (CTM)

def SearchSurface(fixedSystem, mobileSystem, offset=0, theta=45, phi=45, delta=0.27, ndelta=0):
    """
    Parameters
    ----------

    fixedSystem:    mdtraj Obj
                    mdtraj Obj of the system we keep fixed, centered at (0,0,0)

    mobileSystem:   mdtraj Obj
                    mdtraj Obj of the system we intend to rotate


    offset:      float
                    how many nm to add to the offset, to account for
                    vdW interactions between the systems.

    theta:          float
                    increment by which to increase the rotation of
                    sys2 relative to the origin, in degrees

    phi:            float
                    increment by which to increase the rotation of
                    sys2 relative to itself, in degrees
    delta:          float
                    increment by which to translate sys2 relative
                    to the origin
    ndelta          int
                    how many times to translate sys2 by +/- delta
                    relative to the origin

    Returns
    -------
    conformations:  MDTraj.Trajectory
                    array of dimensions (M,N,3), where:
                    *)   M is the number of conformations found, which should be equal to
                    (360/theta + 360/phi + ndelta + 1)
                    *)   N is the number of particles in the system
                    *)   3 is the integer that comes after 2, but before 4.

    Notes
    -----
    By default, we initially increment along the z axis, then translate
    along the Y axis (by delta), then we rotate in the XZ plane
    (by both theta and phi). Afterwards we bring the mobile system
    closer to the origin so the contacts made between the two systems are tight.

    This class works on the assumption that the system was previously centered in
    the (0,0,0) coordinate. This can be done via the AlignSystem class, found
    in this module.
    """

    conformations = mobileSystem[0]
    tempSys = mobileSystem[0]

    # Get the coordinates of the point farthest from the origin in the XZ plane for
    # both systems
    centerChain = fixedSystem.topology.select("chainid {}".format(fixedSystem.n_chains-1))
    #print ("chainid {}".format(fixedSystem.n_chains-1))
    atomDistances_FS = []
    for atomIx in (centerChain):
        #print (atomIx)
        #print (np.array([fixedSystem.xyz[0][atomIx][0], fixedSystem.xyz[0][atomIx][2]]))
        atomDistances_FS.append(np.linalg.norm(np.array([fixedSystem.xyz[0][atomIx][0], fixedSystem.xyz[0][atomIx][2]])))
        #print (fixedSystem.topology.atom(atomIx))
        #print(np.linalg.norm(fixedSystem.xyz[0][atomIx]))
    maxDist_FS = np.max(atomDistances_FS)
    print ("maxDist_FS = {}".format(maxDist_FS))

    atomDistances_MS = np.zeros((mobileSystem.n_atoms))
    for atomIx in range(mobileSystem.n_atoms):
        atomDistances_MS[atomIx] = np.linalg.norm(np.array([mobileSystem.xyz[0][atomIx][0], mobileSystem.xyz[0][atomIx][2]]))
    maxDist_MS = np.max(atomDistances_MS)
    print ("maxDist_MS = {}".format(maxDist_MS))

    #print ("MaxDist_FS + MaxDist_MS + offset = {}".format(maxDist_FS + maxDist_MS + offset))
    #offset = maxDist_FS + maxDist_MS + vdWoffset
    #offset = vdWoffset


    print("\n{0} DEBUG INFO {0}".format(20 * "#"))
    #print (offset)
    #print("Farthest atom FS {}: {}nm".format(atomFar,maxDist_FS))
    #print("Farthest atom MS {}: {}nm".format(atomFar,maxDist_MS))
    #print("maxDist_FS + maxDist_MS: {}".format(maxDist_FS + maxDist_MS))
    #print("maxDist_FS + maxDist_MS + vdWoffset: {}".format(offset))

    # Set up number of rotations
    if (theta == 0):
        theta = 360
    if (phi == 0):
        phi = 360

    thetaIterations = np.int(360 / theta)
    phiIterations = np.int(360 / phi)

    print ("theta Iterations: {}".format(thetaIterations))
    print ("phi Iterations: {}".format(phiIterations))
    print ("delta Iterations: {}".format((2*ndelta)+1))
    print ("total Iterations: {}".format(thetaIterations*phiIterations*((2*ndelta)+1)))

    # Display warning if theta/phi are not divisors of 360
    if (360 % theta != 0 or 360 % phi != 0):
        warnings.warn("The theta/phi values are not divisors of 360!")

    #print ("chainid {}".format(fixedSystem.n_chains-1))
    #print ("chainid {}".format(fixedSystem.n_chains-1))
    #centerOfRotation = mdtraj.compute_center_of_mass(fixedSystem,
    #                                                 select="chainid {}".format(fixedSystem.n_chains-1))[0]
    #print(centerOfRotation)

    # It's easier to compute the Phi rotation now, then translate
    startPose = mobileSystem.xyz[0]
    phi_vector = R.from_rotvec(np.radians(phi) * np.array([0, 1, 0]))
    for phiIx in range(phiIterations):
        sys2 = startPose
        #print (sys2)
        for rotIx in range(phiIx):
            sys2 = np.array(phi_vector.apply(sys2))

        # Translate along the X axis
        #print ("sys2: {}".format(sys2))
        #print ("Offset: {}".format(offset))
        #sys2 = sys2 + centerOfRotation + [offset, 0, 0]
        sys2 = sys2 + [maxDist_MS + maxDist_FS + offset, 0, 0]
        #sys2 = sys2 + [offset, 0, 0]
        #print ("translated by: {}".format(maxDist_MS + maxDist_FS + offset))
        #print ("sys2: {}".format(sys2))
        #print ("sys2_CM: {}".format(centroidPoints(sys2)))

        # Translate by delta
        for i in range(-ndelta, ndelta + 1):
            sys3 = sys2
            sys3 = sys3 + (np.array([0, 1, 0]) * i * delta)



            for thetaIx in range(thetaIterations):
                # Define the circular TM
                sys4 = np.zeros(sys3.shape)
                #print ("theta * thetaIx: {}".format(theta * thetaIx))
                #print ("sin theta: {}".format(np.sin(theta * thetaIx)))
                #print ("cos theta: {}".format(np.cos(theta * thetaIx)))
                circularTM = circularOrbit(np.radians(theta * thetaIx))
                #print (circularTM)
                #movedPoint = np.matmul(circularTM, [0,0,0])
                #print (movedPoint)
                for atomIx in range(sys3.shape[0]):
                    sys4[atomIx] = np.matmul(circularTM, sys3[atomIx])
                #print (sys4)

                # In order to better approximate the surface of the fixed system, we
                # bring the mobile system closer to the origin by a distance equal to the minimum
                # distance between the two helices minus the vdWoffset
                #distMin = 99
                #for FS_AtomIX in range(fixedSystem.n_atoms):
                    #for MS_AtomIX in range(mobileSystem.n_atoms):
                #distFS_MS = np.linalg.norm(centroidPoints(fixedSystem.xyz[0]) - centroidPoints(sys4))
                #distFS_MS = np.linalg.norm([0,0,0] - centroidPoints(sys4))
                    #if (distFS_MS < distMin):
                        #distMin = distFS_MS
                #distMin = distFS_MS
                ## distMin = 0
                ## tempList = []
                ## for MS_AtomIX in range(len(sys4)):
                ##     tempList.append(np.linalg.norm([sys4[MS_AtomIX][0],
                ##                                     sys4[MS_AtomIX][2]]))
                ## tempList = np.array(tempList)
                ## pointC = sys4[np.argmin(tempList)]
                ## #print (np.argmin(tempList), pointC, np.min(tempList))
                ## distMin = distMin + np.min(tempList)
                ## #print("Dist min: {}".format(distMin))
                ## tempList = []
                ## #for FS_AtomIX in range(fixedSystem.n_atoms):
                ## for atomIx in centerChain:
                ##     #print ("atomIx centerChain: {}".format(atomIx))
                ##     tempList.append(np.linalg.norm([pointC[0] - fixedSystem.xyz[0][atomIx][0],
                ##                                     pointC[2] - fixedSystem.xyz[0][atomIx][2]]))
                ## tempList = np.array(tempList)
                ## pointD = fixedSystem.xyz[0][np.argmin(tempList)+centerChain[0]]
                ## #print(np.argmin(tempList), pointD, np.min(tempList))
                ## #print (np.argmin(tempList)+centerChain[0], pointD, np.min(tempList))
                ## #print ("pointD: {}".format(np.linalg.norm(pointD)))
                ## #print ("Dist min - point d: {}".format(distMin - np.linalg.norm(pointD)))
                ## distMin = distMin - np.linalg.norm([pointD[0],pointD[2]])
                ## distMin = distMin - offset
                ## #print ("DistMin - offset: {}".format(distMin))
                ## angle = np.radians(theta * thetaIx)
                ## sys4 = sys4 - [distMin*np.cos(angle), 0, distMin*np.sin(angle)]

                # Need to check if the mobile helix collides with any part of
                # the fixed helices. If yes, discard that conformation. If not,
                # save it.
                collision = 0
                for chainIx in range(fixedSystem.n_chains):
                    centerChain = fixedSystem.topology.select("chainid {}".format(chainIx))
                    #print (centerChain)
                    #print (centerChain[0], centerChain[-1])
                    #for FS_AtomIX in range(fixedSystem.n_atoms):
                    #    for MS_AtomIX in range(mobileSystem.n_atoms):
                    #        distFS_MS = np.linalg.norm(fixedSystem.xyz[0][FS_AtomIX] - sys4[MS_AtomIX])
                    distFS_MS = np.linalg.norm(centroidPoints(fixedSystem.xyz[0][centerChain[0] + 1:centerChain[-1]]) - centroidPoints(sys4))
                    #print ("{} - {} = {}".format(centroidPoints(fixedSystem.xyz[0]), centroidPoints(sys4), distFS_MS))
                    if (distFS_MS < 0.5):
                    #    pass
                        collision = 1
                        break

                if (collision == 0):
                    tempSys.xyz[0] = sys4
                    conformations = conformations.join(tempSys)

    conformations = conformations[1:]
    #print (conformations)
    print ("Conformations: {}".format(conformations))
    print("{0} DEBUG INFO END {0}\n".format(18 * "#"))
    return (conformations)