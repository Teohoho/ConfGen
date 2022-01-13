import numpy as np
from scipy.spatial.transform import Rotation as R
import mdtraj as md

import sys


def alignToAxis(system, axis="y"):
    """

    Parameters
    ----------
    system:  MDTrajTrajectory Object
          Trajectory objects of the system to be aligned
    axis: str
          which axis to align to (x,y,z)

    Returns
    -------
    alignedPos: numpy ndarray
               Array of shape (N,3) corresponding to the system aligned to
               the Y axis
    """

    # Define alignment axis
    axis = axis.lower()
    if axis not in ('x', 'y', 'z'):
        raise ValueError("Axis argument must be one of ['x', 'y', 'z']")
    axes_vectors = {"x": np.array([1,0,0]), "y":np.array([0,1,0]), "z":np.array([0,0,1])}
    alignment_axis = axes_vectors[axis]

    # We assume that the system has a helical structure, and as such
    # assume that there is one line that goes through the center of
    # every helical turn. This is NOT true for a non-helical structure
    # TODO: Think of a way to generalize this
    center_axis_point1 = md.compute_center_of_mass(system, "(resid 0 to 3) and name CA")
    center_axis_point2 = md.compute_center_of_mass(system, "(resid {} to {}) and name CA".
                                                   format(system.n_residues - 4,
                                                          system.n_residues - 1))

    center_axis_unNorm = center_axis_point1 - center_axis_point2
    # Normalize
    center_axis = center_axis_unNorm/np.linalg.norm(center_axis_unNorm)

    # Compute Rotation Vector
    rot_ang = np.arccos(np.dot(center_axis, alignment_axis))
    rot_vec = np.cross(center_axis, alignment_axis)
    rot_vec = R.from_rotvec((rot_vec / np.linalg.norm(rot_vec)) * rot_ang)

    #print (system.xyz[0])
    #print (rot_vec.apply(system.xyz[0]))
    aligned_atoms = rot_vec.apply(system.xyz[0])

    # For simplicity, we also translate the system so the
    # center of the axis is at (0,0,0)
    center_axis_unNorm = rot_vec.apply(center_axis_unNorm)
    #print (center_axis_unNorm)
    for AtomIx in range(aligned_atoms.shape[0]):
        offset = np.array([0, (center_axis_unNorm[0][1])/2, 0])
        aligned_atoms[AtomIx] = aligned_atoms[AtomIx] + offset
        #print(aligned_atoms[AtomIx])
        #aligned_atoms[AtomIx] = aligned_atoms[AtomIx]


    return (aligned_atoms)

