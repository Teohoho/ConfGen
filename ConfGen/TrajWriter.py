# Function that takes an array of positions and generates
# an INPCRD file or a DCD trajectory file
import os
import mdtraj as md
import mdtraj.formats


def writeTraj(positions, out_FN, cell_vectors=None, cell_angles=None):
    """

    Parameters
    ----------
    positions:      numpy.ndarray

            if Array of shape (N,3), write rst7
            if Array of shape (M,N,3), write DCD with M frames

    out_FN:         str
            RST7/DCD output file name
    cell_vectors:   numpy.ndarray
            numpy array, of shape (3,) containing the vectors of the
            periodic cell, in Angstrom
    cell_angles:    numpy.ndarray
            numpy array, of shape (3,) containing the periodic cell
            angles, in degrees

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
    if positions.ndim == 2:
        with mdtraj.formats.AmberRestartFile(out_FN, "w") as f:
            f.write(positions)
        # rst_out.write("AMBER-Style RST file, written by ConfGen.\n")
        # rst_out.write("{:5d}\n".format(positions.shape[0]))
        #
        # for AtomIx in range(positions.shape[0]):
        #     if (AtomIx % 2 == 0):
        #         rst_out.write("{:12.7f} {:12.7f} {:12.7f}".format(positions[AtomIx][0],
        #                                                 positions[AtomIx][1],
        #                                                 positions[AtomIx][2]))
        #     else:
        #         rst_out.write("{:12.7f} {:12.7f} {:12.7f}\n".format(positions[AtomIx][0],
        #                                                   positions[AtomIx][1],
        #                                                   positions[AtomIx][2]))
        
        
        print ("AMBER RST7 file {} succesfully written!".format(out_FN))
    elif positions.ndim == 3:
        with md.formats.DCDTrajectoryFile(out_FN, "w") as f:
            f.write(positions)
        print("CHARMM DCD file {} succesfully written!".format(out_FN))



