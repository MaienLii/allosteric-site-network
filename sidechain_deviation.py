import mdtraj as md
import numpy as np

def get_sidechain_atom_ids(top):
    """Discover the atom IDs for atoms that are in sidechains.

    Looks for atoms that are NOT named N, C, CA, O, HA, H, H1,
    H2, H3, or OXT.

    Parameters
    ----------
    top: md.Topology
        Topology object that supplies names for each atom.

    Returns
    -------
    sc_ids: list
        List of np.ndarray objects, each containing the atom ids belonging
        to each residue.
    """

    SELECTION = ('not (name N or name C or name CA or name O or '
                 'name HA or name H or name H1 or name H2 or name '
                 'H3 or name OXT)')

    sc_ids = []
    for i in range(top.n_residues):
        sstr = f'resid {i} and {SELECTION}'
        try:
            ids = top.select(sstr)
        except RecursionError:
            print("Failed with RecursionError on residue index", i, "with querystring:")
            print('"'+sstr+"'")
            import pickle

            with open('toppickle.top', 'wb') as f:
                pickle.dump(top, f)

            raise

        sc_ids.append(ids)

    return sc_ids

def calculate_sidechain_deviation(traj, sidechain_atom_ids, reference_frame=0):
    """Calculate the absolute deviation of sidechain atoms for each residue.

    Parameters
    ----------
    traj: md.Trajectory
        Loaded trajectory object with coordinates of each frame.
    sidechain_atom_ids: list
        List of arrays containing atom IDs of sidechain atoms for each residue.
    reference_frame: int
        The index of the frame to use as reference for calculating deviation.

    Returns
    -------
    absolute_deviation_array: np.ndarray
        Array containing the mean absolute deviation for each residue in each frame.
    """
    # Get the reference structure (can be the first frame or another specific frame)
    reference_coords = traj.xyz[reference_frame]

    # Initialize an array to store the absolute deviation of each residue's sidechain per frame
    absolute_deviation_array = np.zeros((traj.n_frames, traj.topology.n_residues))

    # Loop through each frame
    for j in range(traj.n_frames):
        # Get the coordinates of the current frame
        current_coords = traj.xyz[j]

        # Loop through the sidechain of each residue
        for i, atom_ids in enumerate(sidechain_atom_ids):
            if len(atom_ids) > 0:  # Ensure the current residue has sidechain atoms
                # Get the coordinates of the sidechain atoms for the current residue
                current_sidechain_coords = current_coords[atom_ids, :]
                reference_sidechain_coords = reference_coords[atom_ids, :]

                # Calculate the absolute deviation
                abs_deviation = np.abs(current_sidechain_coords - reference_sidechain_coords)
                # Calculate the norm of the deviation and store it in the array
                absolute_deviation_array[j, i] = np.linalg.norm(abs_deviation, axis=1).mean()  # Take the mean as the deviation for the residue

    return absolute_deviation_array

# Load the trajectory and topology files
traj = md.load('./traj.nc', top='protein.pdb')

# Get the list of sidechain atom IDs
sidechain_atom_ids = get_sidechain_atom_ids(traj.topology)

# Calculate the absolute deviation of each residue's sidechain per frame
absolute_deviation_array = calculate_sidechain_deviation(traj, sidechain_atom_ids)

# Save the absolute deviation array to a file
np.savetxt('sidechain_deviation.txt', absolute_deviation_array, fmt='%.6f', delimiter=' ')
