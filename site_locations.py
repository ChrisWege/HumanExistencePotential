""" 
Site locations are generated here and saved to the set path.
These site locations are synthetically generated for the purpose of testing the code. 
In a real application, these would be based on actual data.

For the example usage the provided code generates a xlsx-file with site locations that follow a 
logistic curve where the distribution of site location follows a gradient from west to east. 

Depending on preference grid definition and parameters can be changed. With that different configurations 
can be tested as well. 

"""

# UNDER CONSTRUCTION
import numpy as np
import pandas as pd
from scipy.special import expit

# Grid definition
nx = 50
ny = 50
lat = np.linspace(-25, -15, ny) 
lon = np.linspace(15, 25, nx)

# Parameters
n_presence = 100    # number of presence points
midpoint = 0.5      # gradient value at 50% presence
steepness = 10      # curve sharpness (higher = sharper)

np.random.seed(42)

# Generate random lat/lon points uniformly across domain
candidate_lats = np.random.uniform(lat.min(), lat.max(), n_presence * 10)
candidate_lons = np.random.uniform(lon.min(), lon.max(), n_presence * 10)

# Normalize lon to 0-1 range (this is the gradient)
gradient_vals = (candidate_lons - lon.min()) / (lon.max() - lon.min())

# Logistic probability that follows low gradient = high presence
prob_presence = expit(-steepness * (gradient_vals - midpoint))

# Accept/reject based on probability
accepted = np.random.random(len(candidate_lats)) < prob_presence
presence = np.column_stack([candidate_lats[accepted], candidate_lons[accepted]])[:n_presence]

# optional: print commands
print(f"Generated {len(presence)} presence points")
print(f"Lat range: {presence[:, 0].min():.4f} to {presence[:, 0].max():.4f}")
print(f"Lon range: {presence[:, 1].min():.4f} to {presence[:, 1].max():.4f}")

# Export to Excel
df = pd.DataFrame(presence, columns=["Latitude", "Longitude"])
df.to_excel("/data/hescor/akoepke/HEP-WHB/presence_locations.xlsx", index=False)
print("Excel file saved!")
