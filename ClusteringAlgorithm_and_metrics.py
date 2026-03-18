# Create Plots to visualize MC Data from Allpix Squared


# ROOT imports
import ROOT
from ROOT import TFile, gDirectory, gSystem, TClass, TH2D, TCanvas

import os, sys
import argparse
import os.path as path
import math 

from particle import Particle

import h5py
import numpy as np

from ClusteringAlgo import clustering, describe_cluster

# argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-l", metavar='libAllpixObjects', required=False,
                    help="specify path to the libAllpixObjects library (generally in allpix-squared/lib/)) ")
parser.add_argument("-d", metavar='detector', required=True,
                    help="specify your detector name (generally detector1, dut, ...) ")
parser.add_argument("-f", metavar='rootfile', required=True, help="specify path to the rootfile to be processed ")
parser.add_argument("-a", metavar='algo', help="specify clustering algorithm to use ",
                    choices=['agglomerative', 'iterative', 'dbscan', 'connected'],
                    default='agglomerative')
parser.add_argument("-v", "--verbose", help="Toggle verbose settings", action="store_true")
parser.add_argument("-dist", "--distance_threshold", type=float, help="Distance threshold for iterative clustering", default=1.5)
# for hdf5 print
parser.add_argument("--metrics_out", type=str, default=None,
                    help="If set, write cluster metrics to this HDF5 file.")
parser.add_argument("--overwrite", action="store_true",
                    help="Allow overwriting --metrics_out if it already exists.")

args = parser.parse_args()


#------------Input variables------------- 

root_file_name = (str(args.f))
detector_name = (str(args.d))
clustering_algorithm = (str(args.a))
verbose = args.verbose
outDir = os.path.dirname(root_file_name)
dist_threshold = args.distance_threshold
#---------------------------------------- 

# initialize storage + counters
metrics_out = args.metrics_out
overwrite = args.overwrite

metrics_rows = []
clusters_total = 0
skew_skipped_count = 0

# don't override unless allowed
if metrics_out is not None and os.path.exists(metrics_out) and not overwrite:
    print(f"ERROR: {metrics_out} already exists. Use --overwrite or choose a new --metrics_out path.")
    sys.exit(1)


#------------Load Allpix2 libraries for object definitions------------- 

if args.l is not None:  # Try to find Allpix Library
    lib_file_name = (str(args.l))
    if (not os.path.isfile(lib_file_name)):
        print("WARNING: ", lib_file_name, " does not exist, exiting")
        exit(1)

elif os.path.isfile(path.abspath(path.join(__file__ ,"..","..","opt","allpix-squared","lib","libAllpixObjects.so"))): # For native installs
    lib_file_name = path.abspath(path.join(__file__ ,"..","..","opt","allpix-squared","lib","libAllpixObjects.so"))

elif os.path.isfile(path.join(path.sep, "opt","allpix-squared","lib","libAllpixObjects.so")): # For Docker installs
    lib_file_name = path.join(path.sep, "opt","allpix-squared","lib","libAllpixObjects.so")

else:
    print("WARNING: No Allpix Objects Library found, exiting")
    exit(1)


if (not os.path.isfile(lib_file_name)):
    print("WARNING: no allpix library found, exiting (Use -l to manually set location of libraries)")
    exit(1)

# load library and rootfile
gSystem.Load(lib_file_name)
#----------------------------------------------------------------------


#------------Input ROOT File------------- 

if (not os.path.isfile(root_file_name)):
    print("WARNING: " + root_file_name + " does not exist, exiting")
    exit(1)

rootfile = ROOT.TFile(root_file_name)
gDirectory.ls()

#------------------------------------------ 
# Check on detector existence
if not rootfile.GetDirectory("detectors/" + detector_name):
    print("\nDetector does not exist. Please choose one of the following detectors:")
    gDirectory.cd("detectors")
    list_of_keys = gDirectory.GetListOfKeys()
    for key in list_of_keys:
        print(key.GetName())
    exit(1)


# Open Data tree from Allpix2 output files
McParticle = rootfile.Get('MCParticle') # Tree with MC Particle information
McTrack = rootfile.Get('MCTrack') # Tree with MC Track information
PixelCharge = rootfile.Get('PixelCharge') # Tree with charge information on each pixel
PropCharge = rootfile.Get('PropagatedCharge') # Tree with propagated charge information
PixelHit = rootfile.Get('PixelHit') # Tree with pixel hit information

# Histogram for hit map
can = TCanvas()
h = TH2D("hit_histogram", "hit_histogram", 448, 0, 448, 512, 0, 512)


# loop on events (looping over any of the tree and fetching entries from the others)
for iev in range(0, PixelHit.GetEntries()):
    
    #fetch entries for each tree for event iev
    PixelHit.GetEntry(iev)
    PixelCharge.GetEntry(iev)
    McParticle.GetEntry(iev)
    McTrack.GetEntry(iev)

    # Fetch branches for each tree corresponding to the detector name
    PixelCharge_branch = PixelCharge.GetBranch(detector_name)
    PixelHit_branch = PixelHit.GetBranch(detector_name)
    McParticle_branch = McParticle.GetBranch(detector_name)
    McTrack_branch = McTrack.GetBranch("global")

    if (not PixelCharge_branch):
        Warning("WARNING: cannot find PixelCharge branch in the TTree with detector name: " + detector_name + ",  exiting")
        exit(1)
    if verbose: print(' processing event number {0}\r'.format(iev), "out of", PixelHit.GetEntries(), "events",)

    # assign AP2 vectors to branches
    br_pix_charge = getattr(PixelCharge, PixelCharge_branch.GetName())
    br_pix_hit = getattr(PixelHit, PixelHit_branch.GetName())
    br_mc_part = getattr(McParticle, McParticle_branch.GetName())
    br_mc_track = getattr(McTrack, McTrack_branch.GetName())


    # Perform  clustering on pixel hits
    if len(br_pix_hit)==0:
        continue
    elif len(br_pix_hit) == 1:
        clusters_indexes = [0]
        n_clusters = 1
        cluster_hits = [[] for _ in range(n_clusters)]

    else:
        # Build list of (x, y) points
        points = [[pix_hit.getPixel().getIndex().x(), pix_hit.getPixel().getIndex().y()] for pix_hit in br_pix_hit]

        # Call unified clustering function
        clusters_indexes, n_clusters = clustering(points, algo=clustering_algorithm, distance_threshold=dist_threshold)

        cluster_hits = [[] for _ in range(n_clusters)]
    
    #Assemble cluster arrays from indexes returned by algorithm
    for i,pix_hit in enumerate(br_pix_hit):
        cluster_hits[clusters_indexes[i]].append((pix_hit.getPixel().getIndex().x(), pix_hit.getPixel().getIndex().y(), pix_hit.getPixelCharge().getCharge(), pix_hit.getGlobalTime()))
        h.Fill(pix_hit.getPixel().getIndex().x(), pix_hit.getPixel().getIndex().y(), math.fabs(pix_hit.getPixelCharge().getCharge()))

    # ---------------- Metrics computation (per cluster) ----------------
    # these metrics can later be used to generate histograms/plots in order to analyze differences amongst particle type used
    # 4-neighbor connectivity check
    def _neighbor_fraction(xs_int, ys_int):
        coords = set(zip(xs_int, ys_int))
        if not coords:
            return float("nan")
        n_with_neighbor = 0
        for (x, y) in coords:
            if ((x + 1, y) in coords) or ((x - 1, y) in coords) or ((x, y + 1) in coords) or ((x, y - 1) in coords):
                n_with_neighbor += 1
        return n_with_neighbor / len(coords)

    #percentile without numpy errors on tiny lists
    def _percentile(sorted_vals, frac):
        # frac in [0,1], e.g. 0.90 for r90
        if not sorted_vals:
            return float("nan")
        if len(sorted_vals) == 1:
            return float(sorted_vals[0])
        idx = frac * (len(sorted_vals) - 1)
        lo = int(math.floor(idx))
        hi = int(math.ceil(idx))
        if lo == hi:
            return float(sorted_vals[lo])
        w = idx - lo
        return float((1.0 - w) * sorted_vals[lo] + w * sorted_vals[hi])

    #dominant axis (unweighted) and projections
    def _pca_axis_and_proj(xs, ys):
        # returns unit vector (ux, uy) and list of projections s
        n = len(xs)
        if n == 0:
            return (1.0, 0.0), []
        mx = sum(xs) / n
        my = sum(ys) / n
        dx = [x - mx for x in xs]
        dy = [y - my for y in ys]

        # 2x2 covariance (unweighted)
        sxx = sum(d * d for d in dx) / n
        syy = sum(d * d for d in dy) / n
        sxy = sum(dx[i] * dy[i] for i in range(n)) / n

        # eigenvector for largest eigenvalue of [[sxx, sxy],[sxy, syy]]
        # angle formula for 2D PCA:
        # theta = 0.5 * atan2(2*sxy, sxx - syy)
        theta = 0.5 * math.atan2(2.0 * sxy, (sxx - syy))
        ux = math.cos(theta)
        uy = math.sin(theta)

        # projections onto axis
        s = [dx[i] * ux + dy[i] * uy for i in range(n)]
        return (ux, uy), s

    # skewness of a list
    def _skewness(vals):
        n = len(vals)
        if n < 3:
            return float("nan")
        m = sum(vals) / n
        diffs = [v - m for v in vals]
        m2 = sum(d * d for d in diffs) / n
        if m2 == 0:
            return 0.0
        m3 = sum(d * d * d for d in diffs) / n
        return m3 / (m2 ** 1.5)

    # axis asymmetry using median split
    def _axis_asymmetry(proj_s):
        n = len(proj_s)
        if n == 0:
            return float("nan")
        s_sorted = sorted(proj_s)
        med = s_sorted[n // 2] if (n % 2 == 1) else 0.5 * (s_sorted[n // 2 - 1] + s_sorted[n // 2])
        n_pos = sum(1 for v in proj_s if v > med)
        n_neg = sum(1 for v in proj_s if v < med)
        # points exactly at median are ignored
        denom = n_pos + n_neg
        if denom == 0:
            return 0.0
        return abs(n_pos - n_neg) / denom

    # Compute metrics for each cluster in this event
    for cluster_id, clu in enumerate(cluster_hits):
        if not clu:
            continue

        clusters_total += 1

        xs = [int(hh[0]) for hh in clu]
        ys = [int(hh[1]) for hh in clu]
        qs = [float(hh[2]) for hh in clu]
        qabs = [abs(q) for q in qs]

        n_hits = len(xs)

        # Energy
        Q_cluster = sum(qabs)

        # Bounding box and fill factor
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        size_x = (xmax - xmin + 1)
        size_y = (ymax - ymin + 1)
        bbox_area = size_x * size_y
        fill_factor = (n_hits / bbox_area) if bbox_area > 0 else float("nan")

        # Charge-weighted center (weights = |Q|)
        wsum = sum(qabs)
        if wsum > 0:
            cx = sum(xs[i] * qabs[i] for i in range(n_hits)) / wsum
            cy = sum(ys[i] * qabs[i] for i in range(n_hits)) / wsum
        else:
            cx = sum(xs) / n_hits
            cy = sum(ys) / n_hits

        # Center residuals (distances to center)
        rs = [math.hypot(xs[i] - cx, ys[i] - cy) for i in range(n_hits)]
        mean_r = sum(rs) / n_hits
        rms_r = math.sqrt(sum(r * r for r in rs) / n_hits)
        r_sorted = sorted(rs)
        r90 = _percentile(r_sorted, 0.90)

        # Neighbor fraction (4-neighbor)
        neighbor_frac = _neighbor_fraction(xs, ys)

        # axis projections (unweighted) for axis metrics
        _, proj_s = _pca_axis_and_proj([float(x) for x in xs], [float(y) for y in ys])
        axis_asym = _axis_asymmetry(proj_s)

        # Skewness (skip if < 3 hits)
        if n_hits < 3:
            axis_skew = float("nan")
            skew_skipped = 1
            skew_skipped_count += 1
        else:
            axis_skew = _skewness(proj_s)
            skew_skipped = 0

        # Store row (one row per cluster)
        metrics_rows.append({
            "event": int(iev),
            "cluster_id": int(cluster_id),
            "Q_cluster": float(Q_cluster),
            "n_hits": int(n_hits),
            "size_x": int(size_x),
            "size_y": int(size_y),
            "bbox_area": int(bbox_area),
            "fill_factor": float(fill_factor),
            "cx": float(cx),
            "cy": float(cy),
            "mean_r": float(mean_r),
            "rms_r": float(rms_r),
            "r90": float(r90),
            "axis_asym": float(axis_asym),
            "axis_skew": float(axis_skew),
            "skew_skipped": int(skew_skipped),
            "neighbor_frac": float(neighbor_frac),
        })


    # Display cluster information
    if verbose and len(cluster_hits)   >= 1:
        for clu in cluster_hits:
            print("------------------")
            # Per-hit information (not necessary)
            for c in clu:
                print("X: ", c[0], " Y: ", c[1], " Q: ", c[2])

            # Cluster-level characteristics
            describe_cluster(clu)
            print("------------------")


# print number of events processed
print(" ----- processed events (pixelhit):" + str(PixelHit.GetEntries()) + " ----- ")

# ---------------- Write metrics to HDF5 ----------------
if metrics_out is not None:
    if len(metrics_rows) == 0:
        print("WARNING: metrics_out was set but no clusters were recorded (metrics_rows is empty).")
    else:
        # Store run-level data as file attributes
        with h5py.File(metrics_out, "w") as hf:
            hf.attrs["detector"] = detector_name
            hf.attrs["algo"] = clustering_algorithm
            hf.attrs["dist_threshold"] = float(dist_threshold)
            hf.attrs["rootfile"] = root_file_name
            hf.attrs["clusters_total"] = int(clusters_total)
            hf.attrs["skew_skipped_count"] = int(skew_skipped_count)

            # Define a fixed column order for the clusters table
            cols = [
                "event", "cluster_id",
                "Q_cluster", "n_hits",
                "size_x", "size_y", "bbox_area", "fill_factor",
                "cx", "cy",
                "mean_r", "rms_r", "r90",
                "axis_asym", "axis_skew", "skew_skipped",
                "neighbor_frac",
            ]

            # Build a numeric 2D array [n_clusters, n_cols]
            data = np.zeros((len(metrics_rows), len(cols)), dtype=np.float64)
            for r, row in enumerate(metrics_rows):
                data[r, 0]  = float(row["event"])
                data[r, 1]  = float(row["cluster_id"])
                data[r, 2]  = float(row["Q_cluster"])
                data[r, 3]  = float(row["n_hits"])
                data[r, 4]  = float(row["size_x"])
                data[r, 5]  = float(row["size_y"])
                data[r, 6]  = float(row["bbox_area"])
                data[r, 7]  = float(row["fill_factor"])
                data[r, 8]  = float(row["cx"])
                data[r, 9]  = float(row["cy"])
                data[r, 10] = float(row["mean_r"])
                data[r, 11] = float(row["rms_r"])
                data[r, 12] = float(row["r90"])
                data[r, 13] = float(row["axis_asym"])
                data[r, 14] = float(row["axis_skew"])
                data[r, 15] = float(row["skew_skipped"])
                data[r, 16] = float(row["neighbor_frac"])

            dset = hf.create_dataset("clusters", data=data, compression="gzip")
            dset.attrs["columns"] = np.array(cols, dtype=h5py.string_dtype(encoding="utf-8"))

        print(f"Wrote {len(metrics_rows)} cluster rows to {metrics_out}")

# Save histogram to numpy array
nx = h.GetNbinsX()
ny = h.GetNbinsY()

hit_map = np.zeros((ny, nx), dtype=np.float32)

for ix in range(1, nx + 1):
    for iy in range(1, ny + 1):
        hit_map[iy - 1, ix - 1] = h.GetBinContent(ix, iy)

# ---- Save to HDF5 if metrics_out is set ----
if metrics_out is not None:
    with h5py.File(metrics_out, "a") as hf:
        if "run_hit_map" in hf:
            del hf["run_hit_map"]   # allow overwrite cleanly
        hf.create_dataset("run_hit_map", data=hit_map, compression="gzip")

h.Draw("COLZ")