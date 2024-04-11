import polars as pl 


data = pl.read_csv("~/Téléchargements/proteins-2024-04-11.csv")


data = data.select([
    "pdbid", "resolution","topology_subunit","topology_show_in",
    "thickness",	"thicknesserror",	"subunit_segments",	"tilt","tilterror","gibbs","tau"
])


import matplotlib.pyplot as plt
from scipy.stats import pointbiserialr

# Convert boolean to numerical
data = data.with_columns(
    pl.col("thickness").replace({True : 1, False : 0})
)

# Plot
plt.figure(figsize=(10, 6))
plt.scatter(data["topology_show_in"], data["thickness"])
plt.xlabel('Topology Show In')
plt.ylabel('Thickness')
plt.title('Correlation between Topology Show In and Thickness')
plt.savefig("test.png")

# Calculate point biserial correlation
corr, _ = pointbiserialr(data["topology_show_in"], data["thickness"])
print('Point Biserial Correlation: %.3f' % corr)