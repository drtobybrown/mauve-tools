# convert spectral cubes to fits using fits2idia
from make_moments import parse_cmd_line_args
from gistPipeline.initialise import _initialise
from pathlib import Path
import subprocess

# parse arguments from either the command line or user input
try:
    # galaxy name along with config and default directory files from the gist run
    dirPath, args = parse_cmd_line_args()
except:
    input_config_string = input("Please full path to the gist config file:")
    input_defaultDir_string = input("Please full path to the defaultDir file:")

print("Converting FITS to HDF5...")

# get the redshift defined in the gistPipeline config
config = _initialise.readMasterConfig(dirPath.configFile, 1)
config = _initialise.addPathsToConfig(config, dirPath)

# Read the FITS cube in Angstrom
productdir = Path(config["GENERAL"]["OUTPUT"])
print('config["GENERAL"]["OUTPUT"] =',config["GENERAL"]["OUTPUT"])
print('productdir =',productdir)
print('productdir.glob =',productdir.glob)
print('productdir.glob("*_LINEcube.fits") =',productdir.glob("*_LINEcube.fits"))
print('list(productdir.glob("*_LINEcube.fits")) =',list(productdir.glob("*_LINEcube.fits")))
linecubepath = list(productdir.glob("*_LINEcube.fits"))[0]
contcubepath = list(productdir.glob("*_CONTcube.fits"))[0]

# convert each cube from fits 2 hdf5
ret = [
    subprocess.run(
        ["fits2idia", "-o", str(fitsfn).replace(".fits", ".hdf5"), fitsfn],
        capture_output=True,
    )
    for fitsfn in [linecubepath, contcubepath]
]
