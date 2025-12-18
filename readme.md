
# Allpix2 Clustering Example ✅

A small collection of scripts and a notebook demonstrating how to cluster pixel hits produced by Allpix^2 simulations. The main scripts parse Allpix^2 ROOT output, run clustering on pixel hits (several algorithms supported), and print basic cluster characteristics while producing a hit-map histogram using ROOT.

---

## ⚙️ Requirements & Dependencies

Recommended: Python 3.8+ with the following packages installed:

- ROOT / PyROOT (used to read Allpix^2 ROOT files and plotting)
- numpy
- h5py
- scikit-learn
- particle

Typical install (conda recommended for ROOT):

```bash
# Create conda env (example)
conda create -n ap2cluster python=3.10 -y
conda activate ap2cluster

# Install ROOT from conda-forge
conda install -c conda-forge root -y

# Install Python packages
pip install numpy h5py scikit-learn particle
```

Note: `libAllpixObjects.so` (Allpix^2 object library) is required to interpret detector branches; it typically lives in the `lib/` folder of your Allpix^2 install or in `/opt/allpix-squared/lib/` when installed from packages.

---

## 🧭 Usage

Run the main processing script to process an Allpix^2 ROOT file and perform clustering:

```bash
python ClusteringAlgorithm_Allpix2.py -f path/to/file.root -d detectorName [-l /path/to/libAllpixObjects.so] [-a ALGO] [-dist DISTANCE] [-v]
```

Options:
- `-f` : Path to the Allpix^2 ROOT file (required)
- `-d` : Detector name (the detector branch inside the ROOT file, required)
- `-l` : Path to `libAllpixObjects.so` (optional if your library is in a standard location)
- `-a` : Clustering algorithm (choices: `agglomerative`, `iterative`, `dbscan`, `connected`) — default: `agglomerative`
- `-dist` : Distance threshold used by some algorithms (default: `1.5`)
- `-v` : Verbose mode

Example:

```bash
python ClusteringAlgorithm_Allpix2.py -f output_proton_fullsim.root -d detector1 -a dbscan -dist 2.0 -v
```

This will print per-cluster information for each event and draw a ROOT 2D histogram of hits (`TCanvas`).

---

## 🧪 Notebook

Open `ClusteringExample.ipynb` for an interactive demonstration and quick visualization of clustering results.

---

## 🔍 Internals

- `ClusteringAlgo.py` contains the clustering implementations and helper functions:
	- `clustering(points, algo=..., distance_threshold=...)` — unified wrapper
	- `describe_cluster(cluster_hits)` — compute and print basic observables (total charge, size, extents)

- `ClusteringAlgorithm_Allpix2.py` is a full example showing how to read Allpix^2 ROOT files, extract pixel hits, call the clustering wrapper, and visualize results.

---

## 🧠 Algorithms (what they do)

Below are concise descriptions of the clustering methods included in this repository and short notes on when to prefer each one.

- **Agglomerative (single linkage)**
	- Hierarchical clustering that merges points/clusters using the nearest-link (single linkage). In the code this uses scikit-learn's `AgglomerativeClustering` with `linkage='single'` and a `distance_threshold` to stop merges.
	- Good for: clusters where proximity (minimum pairwise distance) is the right notion of grouping and you don't want to pre-specify the number of clusters.
	- Parameters: `distance_threshold` (float). Complexity: typically O(n^2) or worse for naive implementations.

- **DBSCAN**
	- Density-based clustering (scikit-learn `DBSCAN`). Groups points into dense regions separated by lower-density areas; it can discover arbitrarily-shaped clusters and labels noise explicitly.
	- Good for: finding irregularly-shaped clusters and isolating noise/outliers.
	- Parameters: `eps` (in our wrapper we use `distance_threshold` → `eps`), `min_samples` (we use `min_samples=1` by default, meaning singletons are allowed). Complexity: O(n log n) with spatial indexing, or O(n^2) worst-case.

- **Iterative distance-based**
	- A simple greedy algorithm: pick a seed point, add all points within `distance_threshold`, then repeat growing the cluster until no more points are close enough.
	- Good for: small datasets, simple tests, or when you want an easy-to-understand threshold-based grouping.
	- Parameters: `distance_threshold`. Complexity: O(n^2) because distances between many point pairs may be compared.

- **Connected (grid-based connectivity)**
	- Grid-based connected-component clustering using 4-neighbor adjacency on integer pixel indices (pixels connected by up/down/left/right steps are in the same cluster).
	- Good for: pixelated detectors where exact grid connectivity matters (no tolerance for fractional offsets).
	- Complexity: O(n + V) where V is number of grid neighbors visited (essentially linear in visited pixels).


## 🚑 Troubleshooting

- "Cannot find `libAllpixObjects.so`": provide its path via `-l` or install Allpix^2 in a standard location.
- "PyROOT not found": install ROOT (conda-forge is easiest) and ensure your environment loads PyROOT.
- Missing Python packages: install via `pip install <package>`.

---



