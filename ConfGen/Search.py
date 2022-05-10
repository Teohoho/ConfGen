import numpy as np
from scipy.spatial.transform import Rotation as R
import warnings, sys

def centroidPoints(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])

def ellipticalOrbit(theta,ratio):
    ETM = np.array([[np.cos(theta), 0, -np.sqrt(ratio)*np.sin(theta)],
                    [0,1,0],
                    [ratio**(-1)*np.sin(theta), 0, np.cos(theta)]])

    return (ETM)

def Search(fixedSystem, mobileSystem, vdWoffset, theta=45, phi=45, delta=0.27, ndelta=0):
    """
    Function that searches the conformational space by translating and rotating
    relative to the Y-axis

    Parameters
    ----------
    fixedSystem:    numpy.ndarray
                    array, shape (N,3), containing the coordinates of the
                    system we want to keep fixed, centered at (0,0,0)

    mobileSystem:   numpy.ndarray
                    array, shape (N,3), of the rotated positions
                    of the two systems. in nm!

    vdWoffset:       float
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
    ndelta:         int
                    how many times to translate sys2 by +/- delta
                    relative to the origin

    Returns
    -------
    conformations: numpy.ndarray
                   array of dimensions (M,N,3), where:
                   *)   M is the number of conformations found, which should be equal to
                   (360/theta + 360/phi + ndelta + 1)
                   *)   N is the number of particles in the system
                   *)   3 is the integer that comes after 2, but before 4.

    Notes
    -----
    By default, we initially increment along the z axis, then translate
    along the Y axis (by delta), then we rotate in the XZ plane
    (by both theta and phi)
    This class works on the assumption that the system was previously centered in
    the (0,0,0) coordinate. This can be done via the AlignSystem class, found
    in this module.
    r1 and r2 are the larger and respectively smaller distance between the (0,0,0)
    point and the points with the largest X/Z coordinate
    """

    conformations = np.zeros((0, mobileSystem.shape[0], 3))

    # Get the coordinates of the point with the largest X and Z
    # coordinate for FixedSystem and MobileSystem
    print ("\n{0} DEBUG INFO {0}".format(20*"#"))
    FS_X = fixedSystem[np.argmax(np.abs(fixedSystem[:,0]))]
    FS_Z = fixedSystem[np.argmax(np.abs(fixedSystem[:,2]))]
    MS_X = mobileSystem[np.argmax(np.abs(mobileSystem[:,0]))]
    MS_Z = mobileSystem[np.argmax(np.abs(mobileSystem[:,2]))]

    print ("FS_X: {} (index {})".format(FS_X,np.argmax(np.abs(fixedSystem[:,0]))))
    print ("FS_Z: {} (index {})".format(FS_Z,np.argmax(np.abs(fixedSystem[:,2]))))
    print ("MS_X: {} (index {})".format(MS_X,np.argmax(np.abs(mobileSystem[:,0]))))
    print ("MS_Z: {} (index {})".format(MS_Z,np.argmax(np.abs(mobileSystem[:,2]))))

    # Get the distance between the points and the origin and define
    # r1 and r2 (in XZ plane only!)
    FS_XDist = np.linalg.norm(np.array([FS_X[0],FS_X[2]]))
    FS_ZDist = np.linalg.norm(np.array([FS_Z[0],FS_Z[2]]))
    MS_XDist = np.linalg.norm(np.array([MS_X[0],MS_X[2]]))
    MS_ZDist = np.linalg.norm(np.array([MS_Z[0],MS_Z[2]]))

    print ("FS_XDist: {}".format(FS_XDist))
    print ("FS_ZDist: {}".format(FS_ZDist))
    print ("MS_XDist: {}".format(MS_XDist))
    print ("MS_ZDist: {}".format(MS_ZDist))

    FS_R1 = max([FS_XDist,FS_ZDist])
    FS_R2 = min([FS_XDist,FS_ZDist])
    MS_R1 = max([MS_XDist,MS_ZDist])
    MS_R2 = min([MS_XDist,MS_ZDist])

    print ("FS_R1: {}".format(FS_R1))
    print ("FS_R2: {}".format(FS_R2))
    print ("MS_R1: {}".format(MS_R1))
    print ("MS_R2: {}".format(MS_R2))

    LongAxis = FS_R1 + FS_R2
    ShortAxis = max([FS_R1 + MS_R2, MS_R1 + FS_R2])
    AxesRatio = LongAxis/ShortAxis

    print ("Long Axis: {}".format(LongAxis))
    print ("Short Axis: {}".format(ShortAxis))
    print ("Axes Ratio: {}".format(AxesRatio))

    offset = np.sum(LongAxis,vdWoffset)
    print ("Offset: {}".format(offset))
    print ("{0} DEBUG INFO END {0}\n".format(16*"#"))


    # Set up number of rotations
    thetaIterations = np.int(360/theta)
    phiIterations = np.int(360/phi)

    # Display warning if theta/phi are not divisors of 360
    if (360 % theta != 0 or 360 % phi != 0):
        warnings.warn("The theta/phi values are not divisors of 360!")

    # It's easier to compute the Phi rotation now, then translate
    startPose = mobileSystem
    phi_vector = R.from_rotvec(np.radians(phi) * np.array([0, 1, 0]))
    for phiIx in range(phiIterations):
        sys2 = startPose
        for rotIx in range(phiIx):
            sys2 = np.array(phi_vector.apply(sys2))

    # Translate along the selected axis
        sys2 = sys2 + np.array(np.array([1,0,0]) * offset)

    # Translate by delta
        for i in range(-ndelta, ndelta+1):
            sys3 = sys2
            sys3 = sys3 + (np.array([0,1,0])*i*delta)

            # Define a center point which we move in an elliptical orbit. We don't
            # move the entire protein to avoid the skewing that comes with an
            # elliptical motion
            centralPoint = np.array([1,0,0]) * offset

            for thetaIx in range(thetaIterations):
                # Define the elliptical TM
                sys4 = np.zeros(sys3.shape)
                ellipticalTM = ellipticalOrbit(np.radians(theta*thetaIx), AxesRatio)
                movedPoint = np.matmul(ellipticalTM,centralPoint)
                for atomIx in range(sys3.shape[0]):
                    sys4[atomIx] = sys3[atomIx] + movedPoint - centralPoint
                conformations = np.concatenate((conformations,np.array(sys4,ndmin=3)))

    return (conformations)
