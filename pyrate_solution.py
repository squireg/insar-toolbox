import os
import os.path
import subprocess
import sys
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET


CONF_FILE = r"""
# PyRate configuration file for GAMMA-format interferograms
#
#------------------------------------
# input/output parameters

# Directory for the (unwrapped) interferograms.
obsdir:       ${obs_dir}

# File containing the list of interferograms to use.
ifgfilelist:  ${ifg_list}

# The DEM file used in the InSAR processing
demfile:      ${dem_file}

# The DEM header file from GAMMA (*.par) or ROI_PAC (*.rsc).
demHeaderFile: ${dem_header_file}

# GAMMA only: The directory containing GAMMA slc.par header files for all epochs
slcFileDir:   ${slc_file_dir}

# GAMMA only: File listing the pool of available slc.par header files
slcfilelist: ${slc_file_list}

# Directory containing the coherence files. If not provided, obsdir will be used.
cohfiledir: ${coh_file_dir}

# File listing the pool of available coherence files.
cohfilelist: ${coh_file_list}

# Directory to write the outputs to
outdir:       {0}

# InSAR processing software: ROI_PAC = 0, GAMMA = 1
processor:    1

# No data averaging threshold for prepifg
noDataAveragingThreshold: 0.5

# The no data value in the interferograms
noDataValue:  0.0

# Nan conversion flag. Set to 1 if missing (0) phase values are converted to nan
nan_conversion: 1

#-----------------------------------
# Multi-threading parameters: used by linrate/timeseries/prepifg
# gamma prepifg runs in parallel on a single machine if parallel != 0
# parallel = 1, linrate/timeseries computation is done in parallel by rows
# parallel = 2, linrate/timeseries computation is done in parallel for each pixel
# parallel = 0, linrate/timeseries computation is done in serial pixel by pixel
parallel:  0
processes: 1

#------------------------------------
# Coherence masking options: used by process
# cohmask: 1 = ON, 0 = OFF
# cohthresh: coherence threshold value, between 0 and 1
cohmask:   ${cohmask}
cohthresh:  ${cohthresh}

#------------------------------------
# Interferogram multi-look and crop options
# ifgcropopt: 1 = minimum 2 = maximum 3 = customise 4 = all ifms already same size
# ifglksx/y: multi-look/subsampling factor in east and north direction respectively
# ifgxfirst,ifgyfirst: x,y of top-left corner
# ifgxlast,ifgylast: x,y of bottom-right corner
ifgcropopt:   3
ifglksx:      ${ifglksx}
ifglksy:      ${ifglksy}
ifgxfirst:    ${ifgxfirst}
ifgxlast:     ${ifgxlast}
ifgyfirst:    ${ifgyfirst}
ifgylast:     ${ifgylast}

#------------------------------------
# Reference pixel search options
# refx/y: Lon/Lat coordinate of reference pixel. If left blank then search for best pixel will be performed
# refnx/y: number of search grid points in x/y direction
# refchipsize: chip size of the data window at each search grid point
# refminfrac: minimum fraction of valid (non-NaN) pixels in the data window
refx:          ${refx}
refy:          ${refy}
refnx:         ${refnx}
refny:         ${refny}
refchipsize:   ${refchipsize}
refminfrac:    ${refminfrac}

#------------------------------------
# Reference phase calculation method
# refest: 1 = median of the whole interferogram
# refest: 2 = median within the window surrounding the chosen reference pixel
refest:        2

#------------------------------------
# Orbital error correction
# orbfit: ON = 1, OFF = 0
# orbfitmethod = 1: independent method; 2: network method
# orbfitdegrees: Degree of polynomial surface to fit (1 = planar; 2 = quadratic; 3 = part-cubic)
# orbfitlksx/y: additional multi-look factor for orbital correction
orbfit:        ${orbfit}
orbfitmethod:  2
orbfitdegrees: ${orbfitdegrees}
orbfitlksx:    10
orbfitlksy:    10

#------------------------------------
# APS correction using spatio-temporal filter
# apsest: ON = 1, OFF = 0
# Spatial low-pass filter parameters
# slpfmethod: filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 (greater than zero) in km for both butterworth and gaussian filters
# slpforder: order n for butterworth filter (default 1)
# slpnanfill: 1 for interpolation, 0 for zero fill
# slpnanfill_method: linear, nearest, cubic; only used when slpnanfill=1
# Temporal low-pass filter parameters
# tlpfmethod: 1 = Gaussian, 2 = Triangular, 3 = Mean filter
# tlpfcutoff: cutoff t0 for gaussian filter in year;
# tlpfpthr: valid pixel threshold;
apsest:         ${apsest}
slpfmethod:     2
slpfcutoff:     ${slpfcutoff}
slpforder:      1
slpnanfill:     1
slpnanfill_method:  cubic
tlpfmethod:   1
tlpfcutoff:   ${tlpfcutoff}
tlpfpthr:     1

#------------------------------------
# Time Series Calculation
# tscal: YES = 1, NO = 0
# tsmethod: Method for time series inversion (1 = Laplacian Smoothing; 2 = SVD)
# smorder: order of Laplacian smoothing operator (1 =  first-order difference; 2 = second-order difference)
# smfactor: smoothing factor for Laplacian smoothing (value provided is converted as 10**smfactor)
# ts_pthr: valid observations threshold for time series inversion
tscal:         ${tscal}
tsmethod:      2
smorder:       2
smfactor:      -0.25
ts_pthr:       ${ts_pthr}

#------------------------------------
# Stacking calculation
# pthr: minimum number of coherent ifg connections for each pixel
# nsig: n-sigma used as residuals threshold for iterativel least squares stacking
# maxsig: maximum residual used as a threshold for values in the rate map
nsig:          ${nsig}
pthr:          ${pthr}
maxsig:        10000
"""

# Tile dataset to grab pyrate tile data from.
insar_tiles = "${insar_tiles}"

# Namespaces for wfs/gml/etc
NS_GML1 = "http://www.opengis.net/gml"
NS_WFS1 = {
    "wfs": "http://www.opengis.net/wfs",
    "gml": NS_GML1,
    "insar": "http://csiro.au/insar"
}
PATH_WFS1 = "./gml:featureMembers//insar:S1_descending_frames_Data61"

NS_GML2 = "http://www.opengis.net/gml/3.2"
NS_WFS2 = {
    "wfs": "http://www.opengis.net/wfs/2.0",
    "gml": NS_GML2,
    "insar": "http://csiro.au/insar"
}
PATH_WFS2 = "./wfs:member//insar:S1_descending_frames_Data61"

tree = ET.parse(insar_tiles)
root = tree.getroot()
if root.tag == "{http://www.opengis.net/wfs}FeatureCollection":
    print("Found WFS v1 FeatureCollection at root.")
    NS = NS_WFS1
    PATH = PATH_WFS1
    NS_GML = NS_GML1
elif root.tag == "{http://www.opengis.net/wfs/2.0}FeatureCollection":
    print("Found WFS v2 FeatureCollection at root.")
    NS = NS_WFS2
    PATH = PATH_WFS2
    NS_GML = NS_GML2
else:
    print("Not a WFS FeatureCollection, giving up.")
    sys.exit(1)

for tile in root.findall(PATH, NS):
    tile_id = tile.get("{{{}}}id".format(NS_GML))
    relorb_element = tile.find("insar:RelOrbit", NS)
    relorb = relorb_element.text if relorb_element is not None else None
    frame_element = tile.find("insar:Frame", NS)
    frame = frame_element.text if frame_element is not None else None
    print("Found tile:", tile_id, relorb, frame)

sys.exit(0)

def run_pyrate_cmd(pyrate_cmd, config_file, *args):
    """Run pyrate_cmd using config_file and any extra args."""
    cmd = ["pyrate", pyrate_cmd, "-f", config_file, *args]
    print("Running PyRate:", cmd)
    return subprocess.run(cmd, check=True, text=True)


# Create a temporary directory for the config and output files.
with TemporaryDirectory() as temp_output_dir:
    print("Created PyRate output directory", temp_output_dir)

    # Create the config file with interpolated values
    conf_file = os.path.join(temp_output_dir, "pyrate_job.conf")
    with open(conf_file, 'w') as f:
        f.write(CONF_FILE.format(temp_output_dir))

    # Run PyRate
    try:
        run_pyrate_cmd("conv2tif", conf_file)
        run_pyrate_cmd("prepifg", conf_file)
        run_pyrate_cmd("process", conf_file)
        run_pyrate_cmd("merge", conf_file)
    except subprocess.CalledProcessError as ex:
        with open(os.path.join(temp_output_dir, "EXCEPTION.txt"), 'w') as f:
            f.write(f"PyRate job failed with exception.\nreturncode: {ex.returncode}\nfailed command: {ex.cmd}\n\nexception: {ex!r}")

    # Work around "cloud" bug in vgl by running in the output dir so we can use
    # relative filenames.
    os.chdir(temp_output_dir)
    # Upload results, files only
    for f in os.listdir(temp_output_dir):
        abs_f = os.path.join(temp_output_dir, f)
        if not os.path.isdir(abs_f):
            # Work around vgl "cloud" bug by using relative filenames
            subprocess.run(["cloud", "upload", f, f])
            # subprocess.run(["cloud", "upload", f, abs_f])
