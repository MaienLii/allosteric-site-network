from enspara.cluster.hybrid import KHybrid
import numpy as np
import matplotlib.pyplot as plt

kybrid = KHybrid(metric='euclidean', n_clusters=None, cluster_radius=, kmedoids_updates=10, mpi_mode=None)

data1 = np.loadtxt('./sidechain_deviation.txt')

kybrid.fit(data1)

np.savetxt('assignment.txt', kybrid.result_.assignments, fmt='%d')





from enspara.msm.timescales import implied_timescales
from enspara.msm.builders import normalize

# Initialize data
data2 = np.loadtxt('./assignment.txt').reshape(1, -1).astype(np.int32)

n_times = 5
lag_times = [i for i in range(1, 100)]

# Calculate all implied timescales
its = implied_timescales(assigns=data2, lag_times=lag_times, method=normalize, n_times=n_times, sliding_window=True, trim=True)

# Plot
plt.figure(figsize=(8, 6))
for i in range(n_times):
    plt.plot(lag_times, its[:, i], label='Timescale {}'.format(i + 1))

plt.xlabel('Lag time / steps')
plt.ylabel('Timescales / steps')
plt.ylim(0, 50000)
plt.grid(True)
plt.savefig('implied_timescales.tiff', dpi=300)
plt.show()





from enspara.msm.msm import MSM
from scipy.spatial.distance import pdist, squareform

# Build the MSM model
msm = MSM(lag_time=40, method=normalize, trim=True)
msm.fit(data2)

# Obtain equilibrium probabilities of states as weights
equilibrium_probabilities = msm.eq_probs_

# Get the trimmed state mapping
trimmed_mapping = msm.mapping_.to_original

# Manually create a mapping from trimmed states to original frame numbers
state_to_frames = {}
for original_frame, state in enumerate(data2.flatten()):
    # Traverse the trimmed mapping to find the matching original state number
    for trimmed_state, original_state in trimmed_mapping.items():
        if original_state == state:
            if trimmed_state not in state_to_frames:
                state_to_frames[trimmed_state] = []
            state_to_frames[trimmed_state].append(original_frame)

# Initialize a list to store results
results1 = []

# Define an RMSD function
def rmsd(a, b):
    return np.sqrt(np.mean((a - b) ** 2))

# Iterate over each state
for state in range(len(equilibrium_probabilities)):
    if state in state_to_frames:
        frames = state_to_frames[state]
        weight = equilibrium_probabilities[state]

        if len(frames) > 0:
            # Extract SASA data for all frames in the state
            sasa_subset = data1[frames]

            # Calculate RMSD between each frame and all others
            rmsd_matrix = np.zeros((len(frames), len(frames)))
            for i in range(len(frames)):
                for j in range(i + 1, len(frames)):
                    rmsd_matrix[i, j] = rmsd(sasa_subset[i], sasa_subset[j])
                    rmsd_matrix[j, i] = rmsd_matrix[i, j]  # Symmetric matrix

            # Calculate the sum of RMSD for each frame to all other frames
            total_rmsd = rmsd_matrix.sum(axis=0)

            # Select the frame with the minimum total RMSD as the cluster center
            center_frame_index = np.argmin(total_rmsd)
            center_frame = frames[center_frame_index]  # Use the original frame number
            results1.append((center_frame, weight))
        else:
            # If frames is a single frame number
            center_frame = frames
            results1.append((center_frame, weight))

result1 = np.array(results1)
# Save results to file
np.savetxt('./weight1.txt', result1, delimiter=' ', fmt='%d %.6f', comments='')

res = result1[:, 1].reshape(-1, 1)

np.savetxt('./weight2.txt', res, delimiter=' ', fmt='%.6f')



results2 = []

# Define an RMSD function
def rmsd(a, b):
    return np.sqrt(np.mean((a - b) ** 2))

for state in range(len(equilibrium_probabilities)):
    if state in state_to_frames:
        frames = state_to_frames[state]

    if len(frames) > 0:
        # Extract SASA data for all frames in the state
        sasa_subset = data1[frames]

        # Calculate RMSD between each frame and all others
        rmsd_matrix = np.zeros((len(frames), len(frames)))
        for i in range(len(frames)):
            for j in range(i + 1, len(frames)):
                rmsd_matrix[i, j] = rmsd(sasa_subset[i], sasa_subset[j])
                rmsd_matrix[j, i] = rmsd_matrix[i, j]  # Symmetric matrix

        # Calculate the sum of RMSD for each frame to all other frames
        total_rmsd = rmsd_matrix.sum(axis=0)

        # Select the frame with the minimum total RMSD as the cluster center
        center_frame_index = np.argmin(total_rmsd)
        center_frame = frames[center_frame_index]  # Use the original frame number
        center_sasa = data1[center_frame]  # Extract the SASA data of the cluster center
        results2.append((center_frame, *center_sasa))
    else:
        # If frames is a single frame number
        center_frame = frames
        center_sasa = data1[center_frame]
        results2.append((center_frame, *center_sasa))

# Save results to file
result2 = np.array(results2)

np.savetxt('./repre1.txt', result2, delimiter=' ', fmt='%d' + ' %.6f' * data1.shape[1], comments='')

res2 = np.delete(result2, 0, axis=1)

np.savetxt('./repre2.txt', res2, delimiter=' ', fmt='%.6f')
