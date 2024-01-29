import optparse
import re
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import fits
from spectral_cube import SpectralCube
from SpectralCubeTools.spectral_cube_tools import identify_signal


def air_to_vac(wavelength):
    """
    Implements the air to vacuum wavelength conversion described in eqn 65 of
    Griesen 2006

    Stolen from specutils with thanks
    https://github.com/astropy/specutils/blob/0933a8e271ab1e333337087433cb683da8ab5a54/specutils/utils/wcs_utils.py#L374C5-L374C15
    """
    wlum = wavelength.to(u.um).value
    return (
        1 + 1e-6 * (287.6155 + 1.62887 / wlum**2 + 0.01360 / wlum**4)
    ) * wavelength


def integrated_intensity_uncertainty(msubcube, noisemap):
    """
    Calculate the uncertainty of the integrated intensity.

    This function takes in a masked subcube `msubcube` and a noise map `noisemap` and calculates the uncertainty of the integrated intensity using Eq. 2 in Brown et al. 2021.

    Parameters:
    - `msubcube` (ndarray): A masked subcube containing the data.
    - `noisemap` (ndarray): A noise map corresponding to the data.

    Returns:
    - `ndarray`: The uncertainty of the integrated intensity.

    Example:
    integrated_intensity_uncertainty(msubcube, noisemap)
    """

    # number of channels in each pixel mask
    N = np.count_nonzero(~np.isnan(msubcube.filled_data[:]), axis=0)

    """
    Eq. 2 in Brown et al. 2021
    u_i = sqrt(N) * sigma_t * delta_V
    """
    # channel width
    delta_v = msubcube.with_spectral_unit(u.Angstrom).spectral_axis.diff()[0]

    return np.sqrt(N) * noisemap * delta_v


def linewidth_uncertainty(msubcube, noisemap):
    """
    Calculate the uncertainty in the linewidth of a spectral cube.

    Parameters:
    - msubcube: A spectral cube representing the data.
    - noisemap: A noise map corresponding to the spectral cube.

    Returns:
    - The uncertainty in the linewidth of the spectral cube.
    """

    # delta_v_line = Nchan * channel width
    N = np.count_nonzero(~np.isnan(msubcube.filled_data[:]), axis=0)
    delta_v_line = N * msubcube.with_spectral_unit(u.km / u.s).spectral_axis.diff()[0]

    mom_s = msubcube.linewidth_sigma()

    # 1/nan and 1/zero (or close to it) is infinte which doesn't make sense in our case
    mom_s[mom_s.value <= 1e-3] = np.nan
    invalid_idx = ~np.isfinite(mom_s)
    invert_mom_s = 1 / mom_s

    mom_I, mom_I_u = moment_convienience(msubcube, noisemap, moment="I")

    """Eq. 4 in Brown et al. 2021
    u_sigv = (u_i / I) * (delta_V^2 / 8*sqrt(5)) * (1/sigma_v)
    """

    term1 = mom_I_u / mom_I
    term2 = delta_v_line**2 / (8 * np.sqrt(5))
    term3 = 1 / mom_s

    term3 = 1 / mom_s
    term3[invalid_idx] = np.nan

    return term1 * term2 * term3


def make_moments(
    cube,
    noisemap,
    line_list,
    line_widths,
    line_names,
    vz,
    brightest_line=6564.614 * u.Angstrom,
    galaxy_width=150 * u.km / u.s,
    moment_list=["I", "V", "S", "P"],
    quicklook=False,
    write_fits=False,
    outdir=".",
    method="singlethresh",  # or option "multithresh"
    run_id="",
):
    """
    Generate moment maps from a multiline cube.

    Args:
        cube: The input data cube.
        noisemap: The noise map of the cube.
        line_list: List of line frequencies.
        line_widths: List of line widths.
        line_names: List of line names.
        vz: The velocity of the target line.
        brightest_line: The rest frequency of the brightest line (default: 6564.614 Å).
        galaxy_width: The width of the galaxy in velocity units (default: 150 km/s).
        moment_list: List of moment types to generate (default: ["I", "V", "S", "P"]).
        quicklook: Whether to generate quicklook plots (default: False).
        write_fits: Whether to write the moment maps to FITS files (default: False).
        outdir: The output directory (default: current directory).
        method: The method to use for identifying signal in the cube (default: "singlethresh").
        run_id: The identifier for the run (default: empty string).
    """

    # Use the brightest line to identify the appropriate peak velocities, but ONLY
    # from a slab including +/- galaxy_width:
    brightest_cube = cube.with_spectral_unit(
        u.km / u.s, rest_value=brightest_line, velocity_convention="optical"
    ).spectral_slab(vz - galaxy_width, vz + galaxy_width)

    # velocity of the brightest pixel
    peak_velocity = brightest_cube.spectral_axis[brightest_cube.argmax(axis=0)]

    # Now loop over EACH line, extracting moments etc. from the appropriate region:
    # we'll also apply a transition-dependent width (line_widths) here because
    # fainter lines may not have peaks as far out as the bright line.
    hdu_list = [fits.PrimaryHDU(header=fits.Header())]

    for line_name, line_freq, line_width in zip(line_names, line_list, line_widths):
        # sub cube centered on target line
        subcube = cube.with_spectral_unit(
            u.km / u.s, rest_value=line_freq, velocity_convention="optical"
        ).spectral_slab(
            peak_velocity.min() - line_width, peak_velocity.max() + line_width
        )

        if method == "singlethresh":
            signalcube = singlethresh_signal_cube(
                subcube, brightest_cube, noisemap, sn_thresh=3, line_width=line_width
            )

        if method == "multithresh":
            noisecube = SpectralCube(
                data=np.repeat(noisemap.array[np.newaxis, :, :], len(subcube), axis=0)
                * noisemap.unit,
                wcs=subcube.wcs,
            )

            signalcube = identify_signal.find_signal_in_cube(
                cube=subcube,
                noisecube=noisecube,
                mask=None,
                nchan_hi=2,
                snr_hi=3,
                nchan_lo=2,
                snr_lo=2,
                prune_by_npix=None,
                prune_by_fracbeam=0.0,
                expand_by_npix=None,
                expand_by_fracbeam=0.0,
                expand_by_nchan=0,
                verbose=False,
            )

        # plt.figure()
        # plt.title("{0}: {1}".format(line_name, line_freq.to_string(format="latex")))
        # signalcube.mean(axis=(1, 2)).quicklook()

        # make & save the moment maps
        for moment in moment_list:
            mom, mom_u = moment_convienience(signalcube, noisemap, moment=moment)

            hdu_list = hdu_list + [
                fits.ImageHDU(
                    data=mom.hdu.data,
                    header=mom.hdu.header,
                    name=line_name + "_" + moment,
                ),
                fits.ImageHDU(
                    data=mom_u.value,
                    header=mom_u.header,
                    name=line_name + "_" + moment + "_unc",
                ),
            ]

            non_detection = np.isnan(mom.hdu.data).all()

            if quicklook and not non_detection:
                quicklook_fn = "{0}/moment_quicklooks/{1}_{4}_{2}_{3}.pdf".format(
                    outdir, run_id, line_name, moment, method
                )
                quicklook_moment(
                    mom, imshow_kwargs={"cmap": mom.default_cmap}, savefile=quicklook_fn
                )

        if non_detection:
            print("{line} not detected.".format(line=line_name))
    if write_fits:
        momfile = "{outdir}/{run_id}_{method}_moments.fits".format(
            outdir=outdir, run_id=run_id, method=method
        )
        fits.HDUList(hdu_list).writeto(momfile, overwrite=True)

    return


def moment_convienience(msubcube, noisemap, moment="I"):
    """
    Calculates a specific moment of a given subcube.

    Parameters:
    - msubcube: The subcube to calculate the moment on.
    - noisemap: The noise map used in the calculations.
    - moment: The specific moment to calculate. Default is "I" for integrated intensity.

    Returns:
    - mom: The calculated moment.
    - mom_u: The uncertainty associated with the calculated moment.
    """
    mom_u = None

    if moment == "I":  # integrated intensity
        mom = msubcube.with_spectral_unit(u.Angstrom).moment(order=0)
        mom.default_cmap = "viridis"
        mom_u = integrated_intensity_uncertainty(msubcube, noisemap)

    elif moment == "V":
        mom = msubcube.moment(order=1)
        mom.default_cmap = "coolwarm"
        mom_u = velocity_field_uncertainty(msubcube, noisemap)

    elif moment == "S":
        mom = msubcube.linewidth_sigma()
        mom.default_cmap = "coolwarm"
        mom_u = linewidth_uncertainty(msubcube, noisemap)

    elif moment == "P":
        mom = msubcube.with_spectral_unit(u.Angstrom).max(axis=0)
        mom.default_cmap = "magma_r"
        mom_u = noisemap.copy()
        mom_u[np.isnan(mom)] = np.nan

    # temp fix so mom_u is written to file
    if mom_u == None:
        mom_u = mom * np.nan

    return mom, mom_u


def parse_cmd_line_args():
    """
    Parse command-line arguments and capture them.

    Returns:
        dirPath (str): The path to the directory specified in the command-line arguments.
        args (list): A list of additional command-line arguments.
    """
    # Capture command-line arguments
    parser = optparse.OptionParser()
    parser.add_option(
        "--galaxy",
        dest="galaxy",
        type="string",
        help="Galaxy name.",
    )
    parser.add_option(
        "--config",
        dest="configFile",
        type="string",
        help="State the path of the config file.",
    )
    parser.add_option(
        "--default-dir",
        dest="defaultDir",
        type="string",
        help="File defining default directories for input, output, configuration files, and spectral templates.",
    )
    (dirPath, args) = parser.parse_args()
    return dirPath, args


def parse_line_freq_vac(lines_air):
    """
    Parses the frequencies of lines in air and converts them to vacuum wavelengths.

    Args:
        lines_air (List[str]): A list of strings representing the frequencies of lines in air.

    Returns:
        List[u.Quantity]: A list of `u.Quantity` objects representing the frequencies of lines
            converted to vacuum wavelengths.
    """

    vac_frequencies = [
        air_to_vac(float(re.findall(r"[-+]?(?:\d*\.*\d+)", line_air)[0]) * u.Angstrom)
        for line_air in lines_air
    ]
    return vac_frequencies


def parse_line_names(lines_air):
    """
    Parse the line names from the given list of lines.

    Parameters:
    - lines_air (list): A list of lines to parse the names from.

    Returns:
    - line_names (list): A list of line names parsed from the given lines.
    """

    line_names = []

    for line_air in lines_air:
        name, _, _ = line_air.partition(
            "."
        )  # format LINEXXXX.XX in where XXXX.XX is angstrom
        line_names.append(name)
    return line_names


def quicklook_moment(mom, imshow_kwargs={"cmap": "viridis"}, savefile=None):
    """
    plot quick look views of the moment maps for QA
    """
    plt.subplot(projection=mom.wcs)
    im = plt.imshow(mom.data, **imshow_kwargs)
    plt.colorbar(im, label=mom.unit.to_string(format="latex"))
    plt.xlabel("RA")
    plt.ylabel("Dec")

    if savefile:
        savefile = Path(savefile)
        plt.title(savefile.stem)
        # make directory if needed
        savefile.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(savefile, bbox_inches="tight")
        print("saved to {}".format(savefile))


def singlethresh_signal_cube(
    subcube, brightest_cube, noisemap, line_width, sn_thresh=3
):
    """
    Generate a masked subcube based on a single threshold value.

    Parameters:
        subcube (SpectralCube): The input subcube.
        brightest_cube (SpectralCube): The cube containing the brightest pixel.
        noisemap (float): The noise map.
        line_width (float): The width of the line.
        sn_thresh (float, optional): The signal-to-noise threshold. Defaults to 3.

    Returns:
        SpectralCube: The masked subcube.
    """

    # velocity of the brightest pixel
    peak_velocity = brightest_cube.spectral_axis[brightest_cube.argmax(axis=0)]

    spatial_mask = brightest_cube.max(axis=0) > (sn_thresh * noisemap)

    # this part makes a cube of velocities for masking work
    temp = subcube.spectral_axis
    velocities = np.tile(temp[:, None, None], subcube.shape[1:])
    # `velocities` has the same shape as `subcube`

    # now we use the velocities from the brightest line to create a mask region
    # in the same velocity range but with different rest frequencies (different
    # lines)
    mask = np.abs(peak_velocity - velocities) < line_width

    # Mask on a pixel-by-pixel basis with a 1-sigma cut
    signal_mask = subcube > noisemap

    # the mask is a cube, the spatial mask is a 2d array, but in this case
    # numpy knows how to combine them properly
    # (signal_mask is a different type, so it can't be combined with the others
    # yet - https://github.com/radio-astro-tools/spectral-cube/issues/231)
    msubcube = subcube.with_mask(mask & spatial_mask).with_mask(signal_mask)

    return msubcube


def vac_to_air(wavelength):
    """
    Griesen 2006 reports that the error in naively inverting Eqn 65 is less
    than 10^-9 and therefore acceptable.  This is therefore eqn 67

    Stolen from specutils with thanks
    https://github.com/astropy/specutils/blob/0933a8e271ab1e333337087433cb683da8ab5a54/specutils/utils/wcs_utils.py#L366C5-L372C15
    """
    wlum = wavelength.to(u.um).value
    nl = 1 + 1e-6 * (287.6155 + 1.62887 / wlum**2 + 0.01360 / wlum**4)
    return wavelength / nl


def velocity_field_uncertainty(msubcube, noisemap):
    """Eq. 3 in Brown et al. 2021
    u_vel = (delta_V / 2*sqrt(3)) * (u_i / I)
    """
    # delta_v_line = Nchan * channel width
    N = np.count_nonzero(~np.isnan(msubcube.filled_data[:]), axis=0)
    delta_v_line = N * msubcube.with_spectral_unit(u.km / u.s).spectral_axis.diff()[0]

    mom_I, mom_I_u = moment_convienience(msubcube, noisemap, moment="I")

    return (delta_v_line / (2 * np.sqrt(3))) * (mom_I_u / mom_I)


if __name__ == "__main__":
    from gistPipeline.initialise import _initialise

    # parse arguments from either the command line or user input
    try:
        # galaxy name along with config and default directory files from the gist run
        dirPath, args = parse_cmd_line_args()
    except:
        input_galaxy_string = input("Please enter the galaxy ID (e.g., NGC4321):")
        input_config_string = input("Please full path to the gist config file:")
        input_defaultDir_string = input("Please full path to the defaultDir file:")

    # format galaxy ID to upper case with no whitespace
    galaxy = dirPath.galaxy.upper().replace(" ", "")

    print("making moments for {galaxy}".format(galaxy=galaxy))

    # get the gistPipeline config
    config = _initialise.readMasterConfig(dirPath.configFile, 1)
    config = _initialise.addPathsToConfig(config, dirPath)

    run_id = config["GENERAL"]["RUN_ID"]

    # redshift used for gist analysis
    redshift = config["GENERAL"]["REDSHIFT"]

    # Read the FITS cube in Angstrom
    productdir = Path(config["GENERAL"]["OUTPUT"])
    linecubepath = list(productdir.glob("*_LINEcube.fits"))[0]
    cube = SpectralCube.read(str(linecubepath), target_chunksize=1e6, use_dask=True).with_spectral_unit(
        u.Angstrom
    )

    # Create a noise map from a line-free region.
    # found this range by eye-balling a spectrum
    # but should automate this in future:
    # s = cube.max(axis=(1,2))
    # s.quicklook()

    # rms
    noisemap = np.sqrt(
        (cube.spectral_slab(5600 * u.Angstrom, 5850 * u.Angstrom) ** 2).mean(axis=0)
    )

    # Lines to be analyzed (including brightest_line)
    lines_air = [
        "Hb4861.35",
        "OIII4958.91",
        "OIII5006.84",
        "NI5197.90",
        "NI5200.26",
        "NII5754.59",
        "HeI5875.61",
        "OI6300.30",
        "SIII6312.06",
        "OI6363.78",
        "NII6548.05",
        "Ha6562.79",
        "NII6583.45",
        "HeI6678.15",
        "SII6716.44",
        "SII6730.82",
    ]

    # line names which by convention are integers in air frequency
    # line frequencies which are in vac frequency (matching spectralCube)
    line_names = parse_line_names(lines_air)
    vac_freqs = parse_line_freq_vac(lines_air)

    line_widths = np.zeros(len(lines_air)) + (
        70 * u.km / u.s
    )  # just use 70 km/s for now, PPXF estimate coming

    brightest_line = 6564.614 * u.Angstrom  # H-alpha rest (VAC)

    # What is the maximum width spanned by the galaxy (in +/- km/s)
    width = 250 * u.km / u.s

    # Velocity center
    vz = redshift * 299792.458 * u.km / u.s

    # moment_list:
    # I -> Intensity
    # V -> Velocity field
    # S -> Linewidth sigma
    # P -> Peak value of spectrum
    moment_list = ["I", "V", "S", "P"]

    # methods:
    # singlethresh: apply simple S/N cut to cube
    # multithresh: signal mask is expanded to lower S/N threshold if
    # spatially adjacent signal spans multiple channels

    make_moments(
        cube=cube,
        noisemap=noisemap,
        line_list=vac_freqs,
        line_widths=line_widths,
        line_names=line_names,
        vz=vz,
        brightest_line=brightest_line,
        galaxy_width=width,
        moment_list=moment_list,
        quicklook=True,
        write_fits=True,
        outdir=str(productdir),
        method="multithresh",
        run_id=run_id,
    )
