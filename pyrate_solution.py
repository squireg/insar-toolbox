import os
import os.path
import subprocess
import sys
from glob import glob
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET
from datetime import date
import re

INSAR_ARD_DIR = "/g/data/dz56/INSAR_ARD/VV/INSAR_ANALYSIS/VICTORIA/S1/GAMMA/"
INSAR_ARD_VARIANT = "_VV_4rlks_eqa"
INSAR_COH_VARIANT = "_VV_4rlks_flat_eqa"
INSAR_INTERVAL_DIR_RE = re.compile("(\d\d\d\d)(\d\d)(\d\d)-(\d\d\d\d)(\d\d)(\d\d)")
INSAR_DEM_RE = re.compile("(\d\d\d\d)(\d\d)(\d\d).*")

# TODO: Use start/end dates to select interferograms, but how to resolve date
# ranges in input files where there are overlaps? E.g for the following three
# ranges, assuming start <= 20180606 and end >= 20180630, should we use [1,3],
# [2], [1,2,3]?
#
# 1. 20180606,20180618
# 2. 20180606,20180630
# 3. 20180618,20180630
#
# Answer from MattG: All of them. Use all ifgs with intervals that are completely within the specified interval.
#
# Dates will be parsed as ISO8601 strings, and assumed to match the timezone of
# the date strings used in the interferogram path/filename structures.
# 
START_DATE = "${startdate}"
END_DATE = "${enddate}"

# Tile dataset to grab pyrate tile data from.
insar_tiles = "${insar_tiles}"

CONF_FILE = r"""
# PyRate configuration file for GAMMA-format interferograms
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Optional ON/OFF switches - ON = 1; OFF = 0

# Coherence masking (PREPIFG)
cohmask:   ${cohmask}

# Orbital error correction (PROCESS)
orbfit:    ${orbfit}

# APS correction using spatio-temporal filter (PROCESS)
apsest:    ${apsest}

# Time series calculation (PROCESS)
tscal:     ${tscal}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Multi-threading parameters used by stacking/timeseries/prepifg
# gamma prepifg runs in parallel on a single machine if parallel = 1
# parallel: 1 = parallel, 0 = serial
parallel:  0
# number of processes
processes: 8

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Input/Output file locations
#
# File containing the list of interferograms to use.
ifgfilelist:   {ifgfilelist}

# The DEM file used in the InSAR processing
demfile:       {demfile}

# The DEM header file from GAMMA (*.par) or ROI_PAC (*.rsc).
demHeaderFile: {demheaderfile}

# File listing the pool of available header files (GAMMA: *slc.par, ROI_PAC: *.rsc)
hdrfilelist:   {hdrfilelist}

# File listing the pool of available coherence files.
cohfilelist:   {cohfilelist}

# Directory to write the outputs to
outdir:        {outdir}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PREPIFG parameters
#------------------------------------
# Input data format: ROI_PAC = 0, GAMMA = 1
processor:    1

# Coherence threshold value for masking, between 0 and 1
cohthresh:  ${cohthresh}

# Multi-look/subsampling factor in east (x) and north (y) dimension
ifglksx:      ${ifglksx}
ifglksy:      ${ifglksy}

# Cropping options
# ifgcropopt: 1 = minimum extent 2 = maximum extent 3 = crop 4 = no cropping
# ifgxfirst,ifgyfirst: longitude (x) and latitude (y) of north-west corner
# ifgxlast,ifgylast: longitude (x) and latitude (y) of south-east corner
ifgcropopt:   3
ifgxfirst:    ${ifgxfirst}
ifgyfirst:    ${ifgyfirst}
ifgxlast:     ${ifgxlast}
ifgylast:     ${ifgylast}

# No-data averaging threshold (0 = 0%; 1 = 100%)
noDataAveragingThreshold: 0.5

# The No-data value used in the interferogram files
noDataValue:  0.0

# Nan conversion flag. Set to 1 if missing No-data values are to be converted to NaN
nan_conversion: 1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PROCESS parameters
#------------------------------------
# Reference pixel search options

# refx/y: Lon/Lat coordinate of reference pixel. If left blank then search for best pixel will be performed
# refnx/y: number of search grid points in x/y image dimensions
# refchipsize: size of the data window at each search grid point
# refminfrac: minimum fraction of valid (non-NaN) pixels in the data window
refx:          ${refx}
refy:          ${refy}
refnx:         ${refnx}
refny:         ${refny}
refchipsize:   ${refchipsize}
refminfrac:    ${refminfrac}

#------------------------------------
# Reference phase correction method

# refest: 1 = median of the whole interferogram
# refest: 2 = median within the window surrounding the chosen reference pixel
refest:        2

#------------------------------------
# Orbital error correction

# orbfitmethod = 1: interferograms corrected independently; 2: network method
# orbfitdegrees: Degree of polynomial surface to fit (1 = planar; 2 = quadratic; 3 = part-cubic)
# orbfitlksx/y: additional multi-look factor for orbital correction
orbfitmethod:  2
orbfitdegrees: ${orbfitdegrees}
orbfitlksx:    10
orbfitlksy:    10

#------------------------------------
# APS spatial low-pass filter parameters

# slpfmethod: Spatial low-pass filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 (greater than zero) in km for both butterworth and gaussian filters
# slpforder: order n for butterworth filter (default 1)
# slpnanfill: 1 for interpolation, 0 for zero fill
# slpnanfill_method: linear, nearest, cubic; only used when slpnanfill=1
slpfmethod:     2
slpfcutoff:     ${slpfcutoff}
slpforder:      1
slpnanfill:     1
slpnanfill_method:  cubic

#------------------------------------
# APS temporal low-pass filter parameters

# tlpfmethod: 1 = Gaussian, 2 = Triangular, 3 = Mean filter
# tlpfcutoff: cutoff t0 for gaussian filter in year;
# tlpfpthr: valid pixel threshold;
tlpfmethod:   1
tlpfcutoff:   ${tlpfcutoff}
tlpfpthr:     1

#------------------------------------
# Time Series Calculation parameters

# tsmethod: Method for time series inversion (1 = Laplacian Smoothing; 2 = SVD)
# smorder: order of Laplacian smoothing operator (1 = first-order difference; 2 = second-order difference)
# smfactor: smoothing factor for Laplacian smoothing (value provided is converted as 10**smfactor)
# ts_pthr: valid observations threshold for time series inversion
tsmethod:      2
smorder:       2
smfactor:     -0.25
ts_pthr:       ${ts_pthr}

#------------------------------------
# Stacking calculation parameters

# pthr: minimum number of coherent ifg connections for each pixel
# nsig: n-sigma used as residuals threshold for iterative least squares stacking
# maxsig: maximum residual used as a threshold for values in the rate map
pthr:          ${pthr}
nsig:          ${nsig}
maxsig:        1000
"""

# Namespaces for wfs/gml/etc
NS_WFS1 = {
    "wfs": "http://www.opengis.net/wfs",
    "gml": "http://www.opengis.net/gml",
    "insar": "http://csiro.au/insar"
}
PATH_WFS1 = "./gml:featureMembers//insar:S1_descending_frames"

NS_WFS2 = {
    "wfs": "http://www.opengis.net/wfs/2.0",
    "gml": "http://www.opengis.net/gml/3.2",
    "insar": "http://csiro.au/insar"
}
PATH_WFS2 = "./wfs:member//insar:S1_descending_frames"


class PyrateException(Exception):
    pass


class InsarInterval(object):
    """Wrapper for an interferogram interval. """
    def __init__(self, name=None, path=None, start=None, end=None,
                 unw=None, tif=None):
        self.name = name
        self.path = path
        self.start = start
        self.end = end
        self.unw = unw
        self.tif = tif


class InsarTile(object):
    """Wrapper for an INSAR tile feature."""

    def __init__(self, tile, ns):
        self.tile = tile
        self.ns = ns
        self._intervals = []

        def _prop(element):
            elem = tile.find(element, namespaces=ns)
            if elem is not None:
                return elem.text
            return f"NO TILE ELEMENT {element}"

        self.gmlid = tile.get(f"{{{ns['gml']}}}id")
        self.relorb = _prop("insar:relorbit")
        self.frame = _prop("insar:frame")

    def scan_ard(self, ard_dir, variant, coh_variant):
        """Scan INSAR files for this tile and return number of interferograms."""
        # Check we have an interferograms directory
        ifgdir = os.path.join(ard_dir, "INT")
        if not os.path.isdir(ifgdir):
            raise PyrateException(f"Interferograms directory not found: {ifgdir}")

        # Find all interferogram directories in the ard_dir, and sort by date
        self._intervals = []
        with os.scandir(ifgdir) as it:
            for entry in it:
                if entry.is_dir():
                    self.load_interval(entry, variant, coh_variant)
        self._intervals = sorted(self._intervals, key=lambda x: x.start)

        # Find the DEM files for this tile
        demdir = os.path.join(ard_dir, "DEM")
        if not os.path.isdir(demdir):
            print("DEM directory not found in standard location.")
        else:
            print("Standard DEM directory found for tile:", demdir)
            # Find latest dem file
            latest = None
            dem_f = None
            dem_h = None
            for f in glob(os.path.join(demdir, f"*{variant}[._]dem.tif")):
                demdate = INSAR_DEM_RE.match(os.path.basename(f))
                if latest is None or demdate > latest:
                    latest = demdate
                    dem_f = f
            if not dem_f:
                print("No DEM file found in DEM directory.")
            else:
                # DEM header can be *.dem.par or *_dem.par, where the * is the
                # same as the stem from the DEM file, but the '.' or '_' before
                # 'dem.par' may be different to the pattern used in the DEM
                # file!
                #
                # Find the stem from dem_f, strip the trailing '.' or '_', and
                # append '[._]dem.par' to find the header file.
                stem = dem_f.rpartition("dem.tif")[0][:-1]
                for c in "._":
                    tryh = f"{stem}{c}dem.par"
                    if os.path.isfile(tryh):
                        dem_h = tryh
                        break
                if dem_h:
                    print("Found standard DEM file:", dem_f)
                    print("Found matching standard DEM header file:", dem_h)
                    self.dem_file = dem_f
                    self.dem_header = dem_h
                else:
                    print("No DEM header file found to match DEM file")

        # Find the SLC header files
        slcdir = os.path.join(ard_dir, "SLC")
        if not os.path.isdir(slcdir):
            print("SLC headers directory not found in standard location.")
        else:
            print("Found standard SLC directory for tile:", slcdir)
            self.slcs = glob(os.path.join(slcdir, "*", "*_VV.slc.par"))

    def load_interval(self, entry, variant, coh_variant):
        """Load the interval contained in DirEntry entry."""
        match = INSAR_INTERVAL_DIR_RE.match(entry.name)
        if match:
            start = date(*map(int, match.group(1,2,3)))
            end = date(*map(int, match.group(4,5,6)))
            interval = InsarInterval(entry.name, entry.path, start, end)
            unw = os.path.join(entry.path, f"{entry.name}{variant}.unw")
            if os.path.isfile(unw):
                interval.unw = unw
            tif = os.path.join(entry.path, f"{entry.name}{variant}.unw.tif")
            if os.path.isfile(tif):
                interval.tif = tif
            else:
                # Try the other tif filename style
                tif = os.path.join(entry.path, f"{entry.name}{variant}_unw.tif")
                if os.path.isfile(tif):
                    interval.tif = tif
            if interval.unw or interval.tif:
                coh = os.path.join(entry.path, f"{entry.name}{coh_variant}.cc.tif")
                if os.path.isfile(coh):
                    interval.coh = coh
                else:
                    # Try underscore style
                    coh = os.path.join(entry.path, f"{entry.name}{coh_variant}_cc.tif")
                    if os.path.isfile(coh):
                        interval.coh = coh
                    else:
                        print("Missing coherence file", coh, "for interferogram:", interval.tif if interval.tif else interval.unw)
                self._intervals.append(interval)
            else:
                print(f"No unwrapped interferogram or tiff found in {entry.path}")
        else:
            print("Directory", entry.name, "Does not match INSAR interval pattern.")

    def intervals(self, start=None, end=None):
        """Return list of intervals for this tile between start and end.

        If start is None return from the earliest available.
        If end is None return to the latest available.

        An interval is returned if it is falls within [start, end], not if it
        overlaps before or after the specified temporal range.

        """
        if self._intervals:
            if start is None:
                start = self._intervals[0].start
            if end is None:
                end = self._intervals[-1].end
        return [i for i in self._intervals if i.start >= start and i.end <= end]


class InsarTileFeatures(object):
    """Parse a feature collection of INSAR tiles."""

    def load(self, tiles_path):
        """Load and parse the collection from file at tiles_path."""
        self.tree = ET.parse(tiles_path)
        self.root = self.tree.getroot()
        if self.root.tag == "{http://www.opengis.net/wfs}FeatureCollection":
            print("Found WFS v1 FeatureCollection at root.")
            self.ns = NS_WFS1
            self.PATH = PATH_WFS1
        elif self.root.tag == "{http://www.opengis.net/wfs/2.0}FeatureCollection":
            print("Found WFS v2 FeatureCollection at root.")
            self.ns = NS_WFS2
            self.PATH = PATH_WFS2
        else:
            raise PyrateException("Not a WFS FeatureCollection, giving up.")

    def tiles(self):
        """Return list of all tiles in the collection."""
        if not self.root:
            raise Exception("PyrateTiles.tiles() called before load().")

        tiles = self.root.findall(self.PATH, self.ns)
        if tiles:
            return [InsarTile(tile, self.ns) for tile in tiles]
        raise PyrateException("No tile found in input dataset.")

    def first(self):
        """Return the first tile in the collection."""
        if not self.root:
            raise Exception("PyrateTiles.first() called before load().")

        tile = self.root.find(self.PATH, self.ns)
        if tile:
            return InsarTile(tile, self.ns)
        raise PyrateException("No tile found in input dataset.")


dataset = InsarTileFeatures()
try:
    # Load tile features from input dataset
    dataset.load(insar_tiles)

    # Grab the first tile
    # TODO: Support processing multiple tiles
    tile = dataset.first()
    tile_id = tile.gmlid
    print("Found tile", tile_id)
except Exception as ex:
    print("Loading INSAR tile features failed:", ex)
    sys.exit(1)

tile_dir = os.path.join(INSAR_ARD_DIR, f"T{tile.relorb}D_F{tile.frame}")
if not os.path.isdir(tile_dir):
    print("Interferograms directory from tile info not found")
    print(tile_dir, " is not a directory.")
    sys.exit(1)
print("INSAR ARD data from", tile_dir)

def run_pyrate_cmd(pyrate_cmd, config_file, *args):
    """Run pyrate_cmd using config_file and any extra args."""
    cmd = ["pyrate", pyrate_cmd, "-f", config_file, *args]
    print("Running PyRate:", cmd)
    return subprocess.run(cmd, check=True, text=True)


# Create a temporary directory for the config and output files.
with TemporaryDirectory() as temp_output_dir:
    print("Created PyRate output directory", temp_output_dir)
    tile.scan_ard(tile_dir, INSAR_ARD_VARIANT, INSAR_COH_VARIANT)
    intervals = tile.intervals()
    has_tifs = all(i.tif for i in intervals)
    if has_tifs:
        print("Found tifs for all interferograms, will skip conv2tif")

    # Generate the ifg list
    ifgfilelist = os.path.join(temp_output_dir, "interferograms.list")
    cohfile = os.path.join(temp_output_dir, "coh_list.txt")
    with open(ifgfilelist, 'w') as f:
        with open(cohfile, 'w') as g:
            for interval in intervals:
                f.write(interval.tif if has_tifs else interval.unw)
                f.write('\n')
                g.write(interval.coh)
                g.write('\n')

    # Find the DEM and its header file
    if not tile.dem_file or not tile.dem_header:
        print("No DEM file or header found in standard location.")
        print("I can't proceed without a DEM file and header specified manually.")
        sys.exit(1)

    # Write the header files list and coherence files list.
    if not tile.slcs:
        print("No SLC header files for tile.")
        print("I can't proceed without a list of header files.")
        sys.exit(1)
    hdrfile = os.path.join(temp_output_dir, "hdr_list.txt")
    with open(hdrfile, 'w') as f:
        for slc in tile.slcs:
            f.write(slc)
            f.write('\n')

    config = dict(ifgfilelist=ifgfilelist,
                  demfile=tile.dem_file,
                  demheaderfile=tile.dem_header,
                  hdrfilelist=hdrfile,
                  cohfilelist=cohfile,
                  outdir=temp_output_dir)

    # Create the config file with interpolated values
    conf_file = os.path.join(temp_output_dir, "pyrate_job.conf")
    with open(conf_file, 'w') as f:
        f.write(CONF_FILE.format(**config))

    # Run PyRate if we have a hopefully valid config
    if all(v is not None for v in config.values()):
        try:
            if not has_tifs:
                run_pyrate_cmd("conv2tif", conf_file)
            run_pyrate_cmd("prepifg", conf_file)
            run_pyrate_cmd("process", conf_file)
            run_pyrate_cmd("merge", conf_file)
        except subprocess.CalledProcessError as ex:
            with open(os.path.join(temp_output_dir, "EXCEPTION.txt"), 'w') as f:
                f.write(f"PyRate job failed with exception.\nreturncode: {ex.returncode}\nfailed command: {ex.cmd}\n\nexception: {ex!r}")

    # Upload results, files only
    # Work around "cloud" bug in vgl by running in the output dir so we can use
    # relative filenames.
    os.chdir(temp_output_dir)
    for f in os.listdir(temp_output_dir):
        abs_f = os.path.join(temp_output_dir, f)
        if not os.path.isdir(abs_f):
            # Work around vgl "cloud" bug by using relative filenames
            subprocess.run(["cloud", "upload", f, f])
            # subprocess.run(["cloud", "upload", f, abs_f])
