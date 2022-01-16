# Function that takes an array of positions and generates
# an INPCRD file or a DCD trajectory file
import os
import numpy as np
import mdtraj as md
import mdtraj.formats


def writeTraj(positions, out_FN, topology=None, b_factors=None):
    """

    Parameters
    ----------
    positions:      numpy.ndarray

            if Array of shape (N,3), write rst7
            if Array of shape (M,N,3), write DCD with M frames

    out_FN:         str
            RST7/DCD output file name

    topology:       MDTrajTrajectory.Topology Object
            Topology of system. Needed for writing PDBs

    cell_vectors:   numpy.ndarray
            numpy array, of shape (3,) containing the vectors of the
            periodic cell, in Angstrom

    cell_angles:    numpy.ndarray
            numpy array, of shape (3,) containing the periodic cell
            angles, in degrees
    b_factors:      list of floats.
            must be equal to number of frames in positions. b_factors to write
            to PDBs. useful for coloring in pymol/vmd

    Returns
    -------
    None

    Notes
    -----
    The INPCRD format states that each line has 6 fields, each one containing
    a float with 12 characters and 7 decimal places (12.7f)
    """

    # Open the out_FN file, after checking if it already exists
    # if (os.path.exists(out_FN)):
    #     raise ValueError("File {} already exists! Please choose a different name "
    #                      "or move to a new directory".format(out_FN))
    positions = positions * 10 # nm to Angstrom
    if ((positions.ndim == 2) or (positions.shape[0] == 1)):
        positions = positions.reshape((positions.shape[-2],positions.shape[-1]))
        if (out_FN[-4:]) in [".rst7", "inpcrd"]:
            with mdtraj.formats.AmberRestartFile(out_FN, "w") as f:
                f.write(positions)
            print ("AMBER RST7 file {} succesfully written!".format(out_FN))

        if (out_FN[-4:]) == ".pdb":
            # Since b_factors must be in the (-10, 100) range, we need to rescale them
            b_factors = np.interp(b_factors, (b_factors.min(), b_factors.max()), (0,99.9))
            currentFrame = md.formats.PDBTrajectoryFile(out_FN, "w")
            print (positions.shape)
            print (topology.n_atoms)
            currentFrame.write(positions, topology=topology,
                           bfactors=np.full(topology.n_atoms,
                                     fill_value = b_factors)
                                )
            print("PDB file {} succesfully written!".format(out_FN))

    elif positions.ndim == 3:
        with md.formats.DCDTrajectoryFile(out_FN, "w") as f:
            f.write(positions)
        print("CHARMM DCD file {} succesfully written!".format(out_FN))




