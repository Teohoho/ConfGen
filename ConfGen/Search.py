import numpy as np
from scipy.spatial.transform import Rotation as R
import warnings
import mdtraj as md
import sys

def Search(pos1, pos2, outRoot, r, theta=45, phi=45, delta=0.27, ndelta=0, axis="z"):
    """
    Function that searches the conformational space.

    Parameters
    ----------

    pos1,pos2:  numpy.ndarray
               array, shape (N,3), of the rotated positions
               of the two systems. in nm!

    outRoot:    str
                root name for all outputted RST files

    r:          float
                initial distance between the two helices,
                which should be large enough to allow free
                rotation of the two helices, so at least the
                sum of the two helices' largest (CA-Sidechain_Atom)
                distance. This is calculated by the TrajectoryLoader
                function.
                TODO: EXPLAIN THIS BETTER, IT DOESN'T MAKE SENSE. FUCK!

    theta:      float
                increment by which to increase the rotation of
                sys2 relative to sys1, in degrees
    phi:        float
                increment by which to increase the rotation of
                sys2 relative to itself, in degrees
    delta:      float
                increment by which to translate sys2 relative
                to sys 1 on the axis they are aligned against
                (in nm)
    ndelta:     int
                how many times to translate sys2 by +/- delta
                relative to sys1
    axis:       str
                which axis to translate on

    Returns
    -------
    conformations: numpy.ndarray
                array of dimensions (M,N,3), where:
                M is the number of conformations found, which should be equal to ()
                #TODO: FILL THIS OUT
                N is the number of particles in the system
                3 is the integer that comes after 2, but before 4.

    Notes
    -----

    By default, we initially increment along the z axis, then translate
    along the Y axis (by delta), then we rotate in the XZ plane
    (by both theta and phi)

    """

    conformations = np.zeros((0,pos2.shape[0],3))

    #print (conformations)
    #print (conformations.shape)

    # Set up number of rotations
    thetaIterations = np.int(360/theta)
    phiIterations = np.int(360/phi)

    #print (thetaIterations,phiIterations)
    # Display warning if theta/phi are not divisors of 360
    if (360 % theta != 0 or 360 % phi != 0):
        warnings.warn("The theta/phi values are not divisors of "
                      "360. A full will not be achieved")

    # Define alignment axis
    axis = axis.lower()
    if axis not in ('x', 'y', 'z'):
        raise ValueError("Axis argument must be one of ['x', 'y', 'z']")
    axes_vectors = {"x": np.array([1, 0, 0]), "y":np.array([0, 1, 0]), "z":np.array([0, 0, 1])}
    translation_axis = axes_vectors[axis]

    # It's easier to just compute one rotation, in case Theta = Phi
    # So first we rotate along Phi
    startPose = np.array(pos2)
    phi_vector = R.from_rotvec(np.radians(phi) * np.array([0, 1, 0]))
    for phiIx in range(phiIterations):
        sys2 = startPose
        for rotIx in range(phiIx):
            sys2 = np.array(phi_vector.apply(sys2))

    # Translate along the selected axis
        sys2 = np.array(sys2 + (translation_axis * r))

    #print (startPose.shape)
    #print ((axes_vectors["y"]*1*delta).shape)
    # Translate by delta
        for i in range(-ndelta, ndelta+1):
            sys3 = sys2
            #print (axes_vectors["y"]*i*delta)
            #print (i)
            sys3 = sys3 + (axes_vectors["y"]*i*delta)
           # print (sys2[0])

        # Rotate by Theta
        # Define a rotation vector, which is the Y axis,
        # with a norm equal to theta (in radians)
            if (theta==phi):
                theta_vector = phi_vector
            else:
                theta_vector = R.from_rotvec(np.radians(theta) * np.array([0, 1, 0]))

            if thetaIterations > 1:
                for thetaIx in range(thetaIterations):
             #print ("Rotating around Y axis by {}".format(thetaIx*phi))
                    sys3 = np.array(theta_vector.apply(sys3))
                    conformations = np.concatenate((conformations,np.array(sys3,ndmin=3)))
            else:
                conformations = np.concatenate((conformations,np.array(sys3,ndmin=3)))


    #print (conformations.shape)
    return (conformations)
    #return conformations