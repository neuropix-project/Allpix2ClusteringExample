import ROOT
from ROOT import TFile, gDirectory, gSystem, TClass, TH2D, TCanvas

import os, sys
import argparse
import os.path as path
import math 

from particle import Particle

import h5py
import numpy as np

# load sklearn clustering algorithm
from sklearn.cluster import AgglomerativeClustering, DBSCAN
dist_threshold = 1.5
clusteringAlgo = AgglomerativeClustering(n_clusters=None, metric='euclidean', memory=None, connectivity=None, compute_full_tree='auto', linkage='single', distance_threshold=dist_threshold)



def describe_cluster(cluster_hits):
    """Given a list of hits in a cluster, compute and print some basic characteristics.

    Each hit is expected to be a tuple of the form (x, y, Q, t, id).
    """
    if not cluster_hits:
        return

    # Unpack the tuple list for convenience
    xs  = [h[0] for h in cluster_hits]
    ys  = [h[1] for h in cluster_hits]
    qs  = [h[2] for h in cluster_hits]

    total_charge   = sum(qs)
    size           = len(cluster_hits)
    size_x         = max(xs) - min(xs) + 1
    size_y         = max(ys) - min(ys) + 1

    # You can extend this dictionary with more observables if needed
    characteristics = {
        "total_charge": total_charge,
        "size": size,
        "size_x": size_x,
        "size_y": size_y,
    }

    print("Cluster summary:")
    charge = characteristics["total_charge"]
    unit = "e" if charge < 0 else "h"
    print("  total charge  :", charge, unit)
    print("  size (n hits):", characteristics["size"])
    print("  size_x       :", characteristics["size_x"])
    print("  size_y       :", characteristics["size_y"])

    # Also return the dictionary in case the caller wants to store it
    return characteristics

# ----------------------------------------------------------------------



# Iterative clustering function
def iterative_cluster(points, distance_threshold=1.5):
    """Simple iterative clustering algorithm based on distance threshold."""
    clusters = []
    cluster_labels = [-1] * len(points)

    def distance(p1, p2):
        return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

    current_cluster = 0
    for i, point in enumerate(points):
        if cluster_labels[i] != -1:
            continue
        cluster_labels[i] = current_cluster
        cluster_points = [point]
        changed = True
        while changed:
            changed = False
            for j, other_point in enumerate(points):
                if cluster_labels[j] == -1:
                    if any(distance(other_point, cp) <= distance_threshold for cp in cluster_points):
                        cluster_labels[j] = current_cluster
                        cluster_points.append(other_point)
                        changed = True
        current_cluster += 1

    return cluster_labels, current_cluster


def connected_cluster(points):

    """Cluster integer (x, y) points by 4-neighbor grid connectivity.

    Assumes points are pixel indices (integers). Two pixels belong to the
    same cluster if they are connected by steps of (±1, 0) or (0, ±1).

    Parameters
    ----------
    points : list of (x, y)

    Returns
    -------
    cluster_labels : list[int]
        Cluster assignment per point.
    n_clusters : int
        Number of clusters found.
    """
    # Map from coordinate to list of point indices (in case of duplicates)
    coord_to_indices = {}
    for idx, (x, y) in enumerate(points):
        coord_to_indices.setdefault((int(x), int(y)), []).append(idx)

    visited = set()
    cluster_labels = [-1] * len(points)
    cluster_id = 0

    # 4-neighbor offsets
    neighbors = [(1, 0), (-1, 0), (0, 1), (0, -1)]

    for coord, indices in coord_to_indices.items():
        if coord in visited:
            continue

        # BFS/DFS from this coordinate
        stack = [coord]
        visited.add(coord)

        # Assign all points at this coordinate to current cluster
        for idx in indices:
            cluster_labels[idx] = cluster_id

        while stack:
            cx, cy = stack.pop()
            for dx, dy in neighbors:
                nx, ny = cx + dx, cy + dy
                ncoord = (nx, ny)
                if ncoord in coord_to_indices and ncoord not in visited:
                    visited.add(ncoord)
                    stack.append(ncoord)
                    for idx in coord_to_indices[ncoord]:
                        cluster_labels[idx] = cluster_id

        cluster_id += 1

    return cluster_labels, cluster_id


def dbscan_cluster(points, distance_threshold=1.5):
    """Cluster points using DBSCAN algorithm.

    Parameters
    ----------
    points : list of (x, y)
    distance_threshold : float
        Maximum distance between two samples for one to be considered
        as in the neighborhood of the other.

    Returns
    -------
    cluster_labels : list[int]
        Cluster assignment per point.
    n_clusters : int
        Number of clusters found.
    """
    db = DBSCAN(eps=distance_threshold, min_samples=1)
    labels = db.fit_predict(points)

    # Remap labels so they are contiguous non-negative integers,
    # turning any noise label (-1) into its own cluster index.
    new_labels = [-1] * len(labels)
    label_map = {}
    next_label = 0
    for i, lab in enumerate(labels):
        if lab == -1:
            # Treat each noise point as its own singleton cluster
            new_labels[i] = next_label
            next_label += 1
        else:
            if lab not in label_map:
                label_map[lab] = next_label
                next_label += 1
            new_labels[i] = label_map[lab]

    return new_labels, next_label


def agglomerative_cluster(points, distance_threshold=1.5):
    """Cluster points using Agglomerative Clustering algorithm.

    Parameters
    ----------
    points : list of (x, y)
    distance_threshold : float
        Distance threshold to apply when forming clusters.

    Returns
    -------
    cluster_labels : list[int]
        Cluster assignment per point.
    n_clusters : int
        Number of clusters found.
    """
    clusteringAlgo.set_params(distance_threshold=distance_threshold)
    cluster_labels = clusteringAlgo.fit_predict(points) 
    return cluster_labels, clusteringAlgo.n_clusters_

# ------------ Unified clustering wrapper -------------

def clustering(points, algo="agglomerative", distance_threshold=dist_threshold):
    """Run a clustering algorithm on a list of (x, y) points.

    Parameters
    ----------
    points : list of (x, y)
        Coordinates of hits to cluster.
    algo : str
        Either "agglomerative" or "iterative" or "dbscan" or "connected".
    distance_threshold : float
        Threshold used by the iterative method.

    Returns
    -------
    clusters_indexes : list[int]
        Cluster assignment per point.
    n_clusters : int
        Number of clusters found.
    """

    if algo == "agglomerative":
        return agglomerative_cluster(points, distance_threshold=distance_threshold)

    elif algo == "iterative":
        return iterative_cluster(points, distance_threshold=distance_threshold)

    elif algo == "dbscan":
        return dbscan_cluster(points, distance_threshold=distance_threshold)

    elif algo == "connected":
        # Grid-based connected-component clustering (4-neighbor)
        return connected_cluster(points)

    else:
        raise ValueError(f"Unknown clustering algorithm: {algo}")


# ------------------------------------------------------
