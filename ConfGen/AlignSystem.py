import numpy as np
from scipy.spatial.transform import Rotation as R
from math import *

def centroidPoints(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])

def alignToAxis(system, axis="y", CenterAxisSele=None):
    """
    Parameters
    ----------
    system:  MDTrajTrajectory Object
            Trajectory objects of the system to be aligned
    axis: str
            which axis to align to (x,y,z)
    CenterAxisSele: list of strs
            By default, we compute the center axis of the helix as
            the vector that connects the centroid of the first 4
            residues' CA atoms and the centroid of the last 4 residues' CA atoms.
            However, this may not always be appropriate, so it's best to let the
            user define his own selections (MDTraj selection)
            between which to define the center axis.

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

    #Check imput
    if (CenterAxisSele is not None):
        if not (isinstance(CenterAxisSele, list)):
            raise ValueError("CenterAxisSele needs to be a list of two strings"
                             ", not {}".format(type(CenterAxisSele)))
        if (len(CenterAxisSele) != 2):
            raise ValueError("CenterAxisSele needs to be a list of two strings"
                             ", not {}".format(len(CenterAxisSele)))

    # Define alignment axis
    axis = axis.lower()
    if axis not in ('x', 'y', 'z'):
        raise ValueError("Axis argument must be one of ['x', 'y', 'z']")
    axes_vectors = {"x": np.array([1, 0, 0]), "y": np.array([0, 1, 0]), "z": np.array([0, 0, 1])}
    alignment_axis = axes_vectors[axis]

    # We assume that the system has a helical structure, and as such
    # assume that there is one line that goes through the center of
    # every helical turn. This is NOT true for a non-helical structure

    # Get indices for atoms that make up the first and last turns
    if (CenterAxisSele is None):
        FirstSele = "(resid 0 to 3) and name CA"
        LastSele = "(resid {} to {}) and name CA".format(system.n_residues - 4,
                                                         system.n_residues - 1)
    else:
        FirstSele = CenterAxisSele[0]
        LastSele = CenterAxisSele[1]

    FirstIndices = system.topology.select(FirstSele)
    LastIndices = system.topology.select(LastSele)
    BackboneIndices = system.topology.select("protein and backbone")

    FirstHelixPositions = []
    LastHelixPositions = []

    for i in range(len(FirstIndices)):
        FirstHelixPositions.append(system.xyz[0][FirstIndices[i]])
    for i in range(len(LastIndices)):
        LastHelixPositions.append(system.xyz[0][LastIndices[i]])

    FirstHelixPositions = np.array(FirstHelixPositions)
    LastHelixPositions = np.array(LastHelixPositions)

    center_axis_point1 = centroidPoints(FirstHelixPositions)
    center_axis_point2 = centroidPoints(LastHelixPositions)

    center_axis_unNorm = center_axis_point1 - center_axis_point2
    # Normalize
    center_axis = center_axis_unNorm/np.linalg.norm(center_axis_unNorm)

    # Compute Rotation Vector
    rot_ang = np.arccos(np.dot(center_axis, alignment_axis))
    rot_vec = np.cross(center_axis, alignment_axis)
    rot_vec = R.from_rotvec((rot_vec / np.linalg.norm(rot_vec)) * rot_ang)

    aligned_atoms = rot_vec.apply(system.xyz[0])

    # Compute the central point of the rotated system, using only the
    # backbone atoms, so there isn't any bias when working with large residues
    BackbonePositions = []
    for i in range(len(BackboneIndices)):
        BackbonePositions.append(aligned_atoms[BackboneIndices[i]])

    # We align to the "axis" axis, but first recompute the central axis
    BackbonePositions = np.array(BackbonePositions)

    centerpoint = centroidPoints(BackbonePositions)

    for AtomIx in range(aligned_atoms.shape[0]):
        X_offset = np.array([centerpoint[0], 0, 0])
        Y_offset = np.array([0,centerpoint[1],0])
        Z_offset = np.array([0, 0, centerpoint[2]])
        aligned_atoms[AtomIx] = ((aligned_atoms[AtomIx] - X_offset) - Y_offset) - Z_offset


    # In addition to centering the system on (0,0,0), we also rotate it so the longest side
    # chain is oriented along the Z axis.
    # First, we find the atom that is farthest away from the (0,0,0), on the XZ plane
    # We do this by computing the max(sqrt(x^2+z^2))
    #if (np.max(np.absolute(aligned_atoms[:,0])) > np.max(np.absolute(aligned_atoms[:,2]))):
        #print(np.argmax(np.absolute(aligned_atoms[:, 0])))
    #    print ("X is largest!")
    #    farthestAtomIx = np.argmax(np.absolute(aligned_atoms[:, 0]))
    #else:
        #print(np.argmax(np.absolute(aligned_atoms[:, 2])))
    #    print("Z is largest!")
    #    farthestAtomIx = np.argmax(np.absolute(aligned_atoms[:, 2]))
    farthestAtomIx = 0
    farthestAtomDistance = 0
    for atomIx in range(len(aligned_atoms)):
        atomDistance = sqrt(aligned_atoms[atomIx][0]**2 + aligned_atoms[atomIx][2]**2)
        if (atomDistance > farthestAtomDistance):
            farthestAtomIx = atomIx
            farthestAtomDistance = atomDistance

    print("Farthest atom found: {} (index {}) ({} nm) ".format(
          system.topology.atom(farthestAtomIx), farthestAtomIx, farthestAtomDistance))
    print(aligned_atoms[farthestAtomIx])
    #print ("Farthest atom from the centerpoint, in the XZ plane, is {} (index {})".
    #       format(aligned_atoms[farthestAtomIx], farthestAtomIx))

    # Compute the angle we need to rotate by to align said atom to the Z axis.
    farthestAtom=aligned_atoms[farthestAtomIx]
    #print (farthestAtom)
    farthestAtom = farthestAtom / np.linalg.norm(farthestAtom)
    rot_ang = np.arccos(abs(farthestAtom[2])/sqrt(farthestAtom[0]**2 + farthestAtom[2]**2))
    print (np.degrees(rot_ang))

    if (farthestAtom[2]>0):
        zsign = -1
    else:
        zsign= 1
    if (farthestAtom[0]>0):
        xsign = 1
    else:
        xsign = -1

    signage = zsign*xsign
    rot_vec = R.from_rotvec(np.array([0,1,0]) * (rot_ang * signage))
    aligned_atoms = rot_vec.apply(aligned_atoms)
    print (aligned_atoms[farthestAtomIx])

    longAxis=abs(aligned_atoms[farthestAtomIx][2])
    # Now we can get the short axis, by getting the atom that is the
    # farthest from the center on the x axis
    farthestAtomIx = np.argmax(np.absolute(aligned_atoms[:, 0]))
    print ("Farthest point on the X axis: {} (index {})".format(
        aligned_atoms[farthestAtomIx],farthestAtomIx))
    shortAxis=abs(aligned_atoms[farthestAtomIx][0])
    print(aligned_atoms[farthestAtomIx])

    # Long ellipsis axis:
    print ("Long axis: {} nm".format(longAxis))
    print ("Short axis: {} nm".format(shortAxis))

    return (aligned_atoms, longAxis, shortAxis)

