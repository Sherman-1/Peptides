
import pandas as pd
import polars as pl 
import numpy as np

from scipy.spatial.distance import euclidean

import plotly.express as px
import plotly.io as pio

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors

import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def lighten_color(color, amount=0.5) -> str:

    try:
        c = mcolors.cnames[color]
    except:
        c = color
    c = mcolors.to_rgba(c)
    return mcolors.to_hex([c[0] + (1.0 - c[0]) * amount,
                           c[1] + (1.0 - c[1]) * amount,
                           c[2] + (1.0 - c[2]) * amount,
                           c[3]])
    
base_colormap = {
    
    "Bitopic_short" : "blue",
    "Bitopic_long" : "red",
    "Polytopic_short" : "green",
    "Polytopic_long" : "orange",
}
colormap = {}

for category, color in base_colormap.items():
    colormap[category] = color
    colormap[f'{category}_reversed'] = lighten_color(color, amount=0.5)

def draw_ellipse(ax, mean, cov, color, linestyle='--') -> None:

    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.sqrt(eigvals)
    
    ellipse = Ellipse(
        mean, width=2 * eigvals[0], height=2 * eigvals[1],
        angle=np.degrees(np.arctan2(*eigvecs[:, 0][::-1])), 
        edgecolor=color, facecolor='none', linestyle=linestyle
    )
    
    ax.add_patch(ellipse)


data = (
    pl.read_parquet("../dev/random_and_reversed_christos.parquet")
    .filter(
        (~pl.col("category").str.contains("_reversed")) & 
        (~pl.col("category").str.contains("Bitopic")) & 
        (pl.col("category") != "Polytopic_long")
    )
)

X = data.filter(pl.col("category") != "random").to_pandas().drop(columns = ["id", "category"])
y = data.filter(pl.col("category") != "random").to_pandas()["category"]
ids = data.filter(pl.col("category") != "random").to_pandas()["id"]




scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("Starting PCA ...")

pca = PCA(n_components = 3)

X_pca = pca.fit_transform(X_scaled)

df = pd.DataFrame(X_pca, columns = ['PC1', 'PC2', 'PC3'])
df['category'] = y
df['id'] = ids  

"""
random_pca = pca.transform(random_scaled)
random_df = pd.DataFrame(random_pca, columns = ['PC1', 'PC2', 'PC3'])
random_df['category'] = 'random'

reversed_pca = pca.transform(reversed_scaled)
reversed_df = pd.DataFrame(reversed_pca, columns = ['PC1', 'PC2', 'PC3'])
reversed_df['category'] = y_rev
"""

fig, ax = plt.subplots()
ax.grid(True)

ax.axvline(x=2, color='red', linestyle='--')

ax.axhline(y=-4, color = 'red', linestyle='--')



categories = df['category'].unique()

for category in categories:
    
    if "Polytopic" not in category:
        continue
    
    print(f"Processing category : {category}")
    
    
    subset = df[df['category'] == category]

    mean = subset[['PC1', 'PC2']].mean().values
    cov = np.cov(subset[['PC1', 'PC2']].T)

    color = colormap[category]
    
    sns.scatterplot(x='PC1', y='PC2', data=subset, alpha=0.1, color=color, label=f'Category {category}', ax=ax, s=10)
        
    draw_ellipse(ax, mean, cov, color=color)
    
handles, labels = ax.get_legend_handles_labels()
for handle in handles:
    handle.set_alpha(1)  
    
explained = pca.explained_variance_ratio_

ax.set_xlabel(f"PC1 ({explained[0]:.2f})")
ax.set_ylabel(f"PC2 ({explained[1]:.2f})")
ax.set_title("Small_vs_Folded")
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("short_vs_long.pdf")
plt.close()

# Extract datapoints that have PC1 > -5 and the ones that have it < -5
# and extract their ids from the data dataframe



population_1 = df[df['PC1'] < 2].id
population_2 = df[(df['PC1'] > 2) & (df['PC2'] < -4)].id
population_3 = df[(df['PC1'] > 2) & (df['PC2'] > -4)].id

population_1.to_csv("population_1.csv", index = False)
population_2.to_csv("population_2.csv", index = False)
population_3.to_csv("population_3.csv", index = False)

print(population_1.head())

#sns.scatterplot(data = random_df, x = 'PC1', y = 'PC2', color = "grey", ax = ax, s = 10, alpha = 0.5, label = "Random") 
"""
reversed_cat = reversed_df['category'].unique()

for category in reversed_cat:
    
    subset = reversed_df[reversed_df['category'] == category]

    mean = subset[['PC1', 'PC2']].mean().values
    cov = np.cov(subset[['PC1', 'PC2']].T)

    color = colormap[category]
    
    #sns.scatterplot(x='PC1', y='PC2', data=subset, alpha=0.1, color=color, label=f'Category {category}', ax=ax, s=10)
        
    draw_ellipse(ax, mean, cov, color=color, linestyle = ':')



print("PCA done and saved")

# Compute the euclidean distance between the original and reversed datapoints

def euclid(list1, list2):
    return np.sqrt(np.sum((np.array(list1) - np.array(list2)) ** 2))

full = pl.scan_parquet("../dev/random_and_reversed_christos.parquet")

data = (
    full
    .with_columns([
        pl.when(pl.col("category").str.contains("_reversed"))
        .then(pl.lit("reversed"))
        .otherwise(pl.lit("original"))
        .alias("type")
    ])
    .with_columns(
        coords = pl.concat_list(pl.exclude(["id","type","category"]))
    )
    .collect()
    .pivot(values = "coords", index = "id", columns = "type")
    .with_columns(
        euclidian_distance = pl.struct(['original', 'reversed']).map_elements(lambda row: euclid(row['original'], row['reversed']), return_dtype = pl.Float32)
    )
    .join(
        full
        .filter(
            (~pl.col("category").str.contains("_reversed"))
        )
        .select(["id", "category"])
        .collect()
    , on = "id", how = "inner")
    .select(["category", "euclidian_distance"])
    .to_pandas()
)

print(data.category.unique())

plt.figure(figsize = (10, 6))   
sns.kdeplot(data = data, x = "euclidian_distance", hue = "category", palette = colormap)
plt.savefig("euclidian_distance.png")

    
"""
