import numpy as np
import mdtraj as md
from scipy.spatial.transform import Rotation as R
from math import *

def alignToAxis(system, NTermTop):
    """
    Parameters
    ----------
    system:  MDTrajTrajectory Object
            Trajectory objects of the system to be aligned

    NTermTop:  string
            Where the N-Terminus should be oriented, towards
            +Y (positive) or towards -Y (negative)

    Returns
    -------
    alignedPos: numpy ndarray
                Array of shape (N,3) corresponding to the system aligned to
                the Y axis
    longAxis:   float
                Length of long axis that makes up the ellipse that surrounds
                the system
    shortAxis:  float
                Length of short axis that makes up the ellipse that surrounds
                the system
    """
    ## Compute inertia tensor (inTens)
    inTens = md.compute_inertia_tensor(system)

    ## Compute eigenvalue/eigenvector of inTens
    eVal, eVec = np.linalg.eig(inTens)

    ## Translate the protein to (0, 0, 0)
    newPositions = system.xyz[0] - md.compute_center_of_mass(system)

    ## Apply the rotation object on the positions of the protein
    ## Multiplied by 10 since it's Angstroms, not nm
    rotObj = R.from_matrix(-eVec.T.reshape(1, 3, 3))
    newPositions = rotObj.apply(newPositions)

    ## To account for the topology of the system, we need to orient
    ## the helices so the C-terminal of helix N and the N-terminal
    ## of helix N+1 are on the same side of the XZ plane
    ## Get the location of the N-Term
    NTerm = newPositions[0][0]
    if (NTerm < 0):
        if (NTermTop.lower() == "negative"):
            rotVec = np.pi/2 * np.array([0,0,1])
        else:
            rotVec = -(np.pi)/2 * np.array([0,0,1])
    else:
        if (NTermTop.lower() == "negative"):
            rotVec = -(np.pi)/2 * np.array([0,0,1])
        else:
            rotVec = np.pi/2 * np.array([0,0,1])

    newPositions = R.from_rotvec(rotVec).apply(newPositions)
    return (newPositions)