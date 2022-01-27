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
    ETM = np.array([[np.cos(theta), 0, -np.sqrt(ratio**(-1))*np.sin(theta)],
                    [0,1,0],
                    [ratio*np.sin(theta), 0, np.cos(theta)]])

    return (ETM)

def Search(pos2, offset, axesRatio, theta=45, phi=45, delta=0.27, ndelta=0):
    """
    Function that searches the conformational space by translating and rotating
    relative to the Y-axis

    Parameters
    ----------

    pos2:       numpy.ndarray
                array, shape (N,3), of the rotated positions
                of the two systems. in nm!

    offset:     float
                measure of how much to offset the system (pos2)
                from the origin (in nm)

    axesRatio:  float
                ratio of the axes for the elliptical motion to be
                used

    theta:      float
                increment by which to increase the rotation of
                sys2 relative to the origin, in degrees
    phi:        float
                increment by which to increase the rotation of
                sys2 relative to itself, in degrees
    delta:      float
                increment by which to translate sys2 relative
                to the origin
    ndelta:     int
                how many times to translate sys2 by +/- delta
                relative to the origin

    Returns
    -------
    conformations: numpy.ndarray
                array of dimensions (M,N,3), where:
                M is the number of conformations found, which should be equal to
                (360/theta + 360/phi + ndelta + 1)
                N is the number of particles in the system
                3 is the integer that comes after 2, but before 4.

    Notes
    -----
    By default, we initially increment along the z axis, then translate
    along the Y axis (by delta), then we rotate in the XZ plane
    (by both theta and phi)
    This class works on the assumption that the system was previously centered in
    the (0,0,0) coordinate. This can be done via the AlignSystem class, found
    in this module.
    """

    conformations = np.zeros((0, pos2.shape[0], 3))

    # Set up number of rotations
    thetaIterations = np.int(360/theta)
    phiIterations = np.int(360/phi)

    # Display warning if theta/phi are not divisors of 360
    if (360 % theta != 0 or 360 % phi != 0):
        warnings.warn("The theta/phi values are not divisors of "
                      "360. A full rotation will not be achieved")


    # Define Z axis (translation)
    translation_axis = np.array([0, 0, 1])

    # It's easier to compute the Phi rotation now, then translate
    startPose = np.array(pos2)
    phi_vector = R.from_rotvec(np.radians(phi) * np.array([0, 1, 0]))
    for phiIx in range(phiIterations):
        sys2 = startPose
        for rotIx in range(phiIx):
            sys2 = np.array(phi_vector.apply(sys2))

    # Translate along the selected axis
        sys2 = np.array(sys2 + (translation_axis * offset))

    # Translate by delta
        for i in range(-ndelta, ndelta+1):
            sys3 = sys2
            sys3 = sys3 + (np.array([0,1,0])*i*delta)

            # Define a center point which we move in an elliptical orbit. We don't
            # move the entire protein to avoid the skewing that comes with an
            # elliptical motion
            centralPoint = np.array(translation_axis * offset)

            for thetaIx in range(thetaIterations):
                # Define the elliptical TM
                sys4 = np.zeros(sys3.shape)
                ellipticalTM = ellipticalOrbit(np.radians(theta*thetaIx), axesRatio)
                movedPoint = np.matmul(ellipticalTM,centralPoint)
                for atomIx in range(sys3.shape[0]):
                    sys4[atomIx] = sys3[atomIx] + movedPoint - centralPoint
                conformations = np.concatenate((conformations,np.array(sys4,ndmin=3)))

    return (conformations)
