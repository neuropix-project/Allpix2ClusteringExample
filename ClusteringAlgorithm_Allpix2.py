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

args = parser.parse_args()


#------------Input variables------------- 

root_file_name = (str(args.f))
detector_name = (str(args.d))
clustering_algorithm = (str(args.a))
verbose = args.verbose
outDir = os.path.dirname(root_file_name)
dist_threshold = args.distance_threshold
#---------------------------------------- 




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
h = TH2D("hit_histogram", "hit_histogram", 512, 0, 512, 448, 0, 448)


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
    elif len(br_pix_hit)==1:
        clusters_indexes = [0]
        cluster_hits = [[]]
    else:
        # Build list of (x, y) points
        points = [[pix_hit.getPixel().getIndex().x(), pix_hit.getPixel().getIndex().y()] for pix_hit in br_pix_hit]

        # Call unified clustering function
        clusters_indexes, n_clusters = clustering(points, algo=clustering_algorithm, distance_threshold=dist_threshold)

        cluster_hits = [[] for _ in range(n_clusters)]
    
    #Assemble cluster arrays from indexes returned by algorithm
    for i,pix_hit in enumerate(br_pix_hit):
        cluster_hits[clusters_indexes[i]].append((pix_hit.getPixel().getIndex().x(), pix_hit.getPixel().getIndex().y(), pix_hit.getPixelCharge().getCharge(), pix_hit.getGlobalTime(), id))
        h.Fill(pix_hit.getPixel().getIndex().x(), pix_hit.getPixel().getIndex().y(), math.fabs(pix_hit.getPixelCharge().getCharge()))
        
    # Display cluster information
    if len(cluster_hits)   >= 1:
        for clu in cluster_hits:
            print("------------------")
            # Per-hit information (optional, keep or remove as you like)
            for c in clu:
                print("X: ", c[0], " Y: ", c[1], " Q: ", c[2])

            # Cluster-level characteristics
            describe_cluster(clu)
            print("------------------")


# print number of events processed
print(" ----- processed events (pixelhit):" + str(PixelHit.GetEntries()) + " ----- ")
h.Draw("COLZ")