import asdf
import numpy as np
from astropy.io import fits
from solid_waffle.correlation_run import run_ir_all
from solid_waffle.flat_simulator import simulate_flat
from solid_waffle.pyirc import get_num_slices, load_segment


def create_dummy_asdf(asdf_path, data_type="flat", frames=20, shape=(512, 512)):
    """
    Test function to make a 512x512 simulated asdf file.

    Parameters
    ----------
    asdf_path : str or pathlib.Path
        Path to the file which will be created
    data_type : str
        Type of file to be simulated - "dark" or "flat"
    frames : int
        number of frames in simulated asdf file. Default is 20
    shape : tuple (int, int)
        shape of data array in simulated asdf file. Default is (512, 512)

    Returns
    -------
    data : np.ndarray of float, shape (frames, shape[0], shape[1])
        data stored in simulated asdf file

    """
    rng = np.random.default_rng(42)  # fixed seed
    if data_type == "flat":
        data = 3000 + rng.normal(0, 50, size=(frames, *shape))
    else:
        data = 100 + rng.normal(0, 5, size=(frames, *shape))
    data = np.clip(data, 0, 65535)
    tree = {"roman": {"data": data}}
    with asdf.AsdfFile(tree) as af:
        af.write_to(asdf_path)
    return data


def test_run_asdf(tmp_path):
    """
    Test function to run solid-waffle on a simulated asdf 512x512 file.

    Parameters
    ----------
    tmp_path : str or pathlib.Path
        Directory in which to run the test.

    Returns
    -------
    None

    """

    frames, ny, nx = 20, 512, 512
    fill_value = 1000.0
    np.full((frames, ny, nx), fill_value)

    asdf_path = str(tmp_path / "test.asdf")
    data_out = create_dummy_asdf(asdf_path, data_type="flat", frames=frames, shape=(ny, nx))

    # Test get_num_slices
    ntslice = get_num_slices(2002, asdf_path)
    assert ntslice == frames

    # Test load_segment
    xyrange = [0, nx, 0, ny]
    tslices = [1, 2, 3]
    result = load_segment(asdf_path, 2002, xyrange, tslices, verbose=False)

    assert result.shape == (len(tslices), ny, nx)
    # Check the 65535 - data transformation was applied
    expected = np.clip(65535 - data_out, 0, 65535)
    assert np.allclose(result, expected[0:3])


def test_run(tmp_path):
    """
    Test function to make a 512x512 simulation and run solid-waffle.

    Parameters
    ----------
    tmp_path : str or pathlib.Path
        Directory in which to run the test.

    Returns
    -------
    None

    """

    temp_dir = str(tmp_path)
    print("using", temp_dir)

    # Make the simulation

    for k in range(8):
        illum = 300.0
        ty = "light"
        if k >= 4:
            illum = 0.0
            ty = "dark"

        nlstr = "quadratic 1.4" if k != 2 else "quartic -1.4 0 0"
        sim_cfg = (
            "FORMAT: 1001\n"
            "NREADS: 20\n"
            "SUBSTEPS: 3\n"
            "DT: 2.75\n"
            "GAIN: 1.5\n"
            f"ILLUMINATION: {illum}\n"
            "QE: 8e-1\n"
            f"RNGSEED: {201909+k}\n"
            "LEGACY\n"  # <-- legacy RNG is stable since numpy 1.16
            "RESET_E: 1.0e2\n"
            "NOISE: Gauss\n"
            "WAVEMODE: ir\n"
            "BFE: true\n"
            "L_IPC: true 0.01\n"
            f"NL: {nlstr}\n"
            f"OUTPUT: {temp_dir}/{ty}_{k+1:03d}.fits\n"
        )

        with open(temp_dir + "/sim_cfg.txt", "w") as f:
            f.write(sim_cfg)
        simulate_flat.run_config(temp_dir + "/sim_cfg.txt")

    # Now analyze it
    analyze_cfg = (
        "DETECTOR: Test_simulation\n"
        "LIGHT:\n"
        f"    {temp_dir}/light_001.fits\n"
        f"    {temp_dir}/light_002.fits\n"
        f"    {temp_dir}/light_003.fits\n"
        f"    {temp_dir}/light_004.fits\n"
        "DARK:\n"
        f"    {temp_dir}/dark_005.fits\n"
        f"    {temp_dir}/dark_006.fits\n"
        f"    {temp_dir}/dark_007.fits\n"
        f"    {temp_dir}/dark_008.fits\n"
        "FORMAT: 1001\n"
        "CHAR: Advanced 1 3 3 bfe\n"
        "NBIN: 4 4\n"
        "TIME:    1 10 12 20\n"
        "TIME2A:  1 2 4 20\n"
        "TIME2B:  1 2 4 20\n"
        "TIME3:   1 2 4 20\n"
        f"OUTPUT: {temp_dir}/analysis\n"
    )
    with open(temp_dir + "/analyze_cfg.txt", "w") as f:
        f.write(analyze_cfg)
    run_ir_all(temp_dir + "/analyze_cfg.txt")

    # Load the analysis
    data = np.loadtxt(temp_dir + "/analysis_summary.txt")
    print(">>", np.mean(data, axis=0))

    # outputs from the first run
    expected_outputs = np.array(
        [
            1.50000000e00,
            1.50000000e00,
            1.55543125e04,
            1.84032388e00,
            1.66757590e00,
            1.49223302e00,
            9.88118358e-03,
            1.01165409e-02,
            1.40400538e-06,
            6.55764439e02,
            2.68351544e-04,
            1.18591758e02,
            1.20411533e02,
            3.16166804e-07,
            -4.55713429e-07,
            2.01035345e-07,
            4.62063881e-08,
            6.72640364e-08,
            6.11719038e-07,
            -2.93234974e-07,
            9.45496306e-08,
            1.16587044e-07,
            4.83671986e-07,
            8.42673374e-08,
            4.60151153e-07,
            -1.68711037e-06,
            8.93645681e-08,
            1.54687202e-07,
            1.12356304e-07,
            -2.24261168e-07,
            1.69978299e-07,
            2.51569092e-07,
            3.37570881e-07,
            -1.48281734e-07,
            4.64848787e-08,
            2.94327429e-07,
            2.93979598e-07,
            -1.48579605e-07,
            6.75570972e-06,
            2.74537451e-06,
            2.67357466e-07,
        ]
    )

    # tolerances -- if anything changes by more than these amounts,
    # pre-commit should warn the user!
    tol = np.array(
        [
            1e-4,
            1e-4,
            2.5,
            0.01,
            0.01,
            0.01,
            1e-4,
            1e-4,
            5e-8,
            0.5,
            1e-4,
            0.1,
            0.1,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            1e-7,
            5e-7,
            5e-7,
            5e-7,
        ]
    )

    diff = np.amax(np.abs(np.mean(data, axis=0) - expected_outputs) / tol)
    print(diff)
    assert diff < 1.0

    # now with IPCSUB on
    with open(temp_dir + "/analyze_cfg.txt", "w") as f:
        f.write(analyze_cfg)
        f.write("IPCSUB: True\n")
    run_ir_all(temp_dir + "/analyze_cfg.txt")
    data = np.loadtxt(temp_dir + "/analysis_summary.txt")
    print(">>", np.mean(data, axis=0))

    expected_outputs[11:13] = np.array([1.22711440e02, 1.16214368e02])
    diff = np.amax(np.abs(np.mean(data, axis=0) - expected_outputs) / tol)
    print(diff)
    assert diff < 1.0

    # hot pixel test. first choose the hot pixels:
    nh = 400
    rate = 10.0 ** (1.5 + np.cos(np.linspace(0, 1000, nh)))
    x = np.zeros(nh, dtype=np.int32)
    y = np.zeros(nh, dtype=np.int32)
    # this part has alpha_adj = 0.01, alpha_diag = 0.001
    cube = np.outer(
        np.linspace(0, 19, 20), np.array([0.001, 0.01, 0.001, 0.01, 0.956, 0.01, 0.001, 0.01, 0.001])
    ).reshape((20, 3, 3))
    for j in range(nh):
        x[j] = 169 * j % 512
        y[j] = 257 * j % 512
        if x[j] <= 4 or x[j] >= 507 or y[j] <= 4 or y[j] >= 507:
            rate[j] = 0.0
    # here we add to the input files
    for k in range(8):
        ty = "light" if k < 4 else "dark"
        fn = f"{temp_dir}/{ty}_{k+1:03d}.fits"
        with fits.open(fn, mode="update") as fk:
            for j in range(nh):
                if rate[j] > 1e-2:
                    r = fk[0].data[:, y[j] - 1 : y[j] + 2, x[j] - 1 : x[j] + 2]
                    r[...] = np.clip(r + cube * rate[j], 0, 65535)
    # and run the analysis
    with open(temp_dir + "/analyze_cfg.txt", "w") as f:
        f.write(analyze_cfg)
        f.write("HOTPIX: 400 2400 0.1 0.08\n")
    run_ir_all(temp_dir + "/analyze_cfg.txt")

    # now see what we got
    hot = np.loadtxt(f"{temp_dir}/analysis_hot.txt")
    (nh2, ncol) = np.shape(hot)
    assert 0.2 < nh2 / nh < 0.4
    assert ncol == 18
    nrec = 0
    for j in range(nh):
        nrec += np.count_nonzero(np.logical_and(x[j] == hot[:, 0], y[j] == hot[:, 1]))
    assert 9 * nh2 // 10 <= nrec <= nh2

    # linear fits
    _x = np.ones((nh2, 2))
    _x[:, 1] = hot[:, -5] + 4 * (hot[:, -4] + hot[:, -1])
    _y = hot[:, -4]
    c, m = np.linalg.lstsq(_x, _y)[0]
    assert -1 < c < 1
    assert 0.009 < m < 0.011
    _y = hot[:, -1]
    c, m = np.linalg.lstsq(_x, _y)[0]
    assert -1 < c < 1
    assert 0.0005 < m < 0.0015

    # and with SLIDEMED and TIMEREF: 1
    with open(temp_dir + "/analyze_cfg.txt", "w") as f:
        f.write(analyze_cfg)
        f.write("IPCSUB: True\n")
        f.write("HOTPIX: 400 2400 0.1 0.08\n")
        f.write("HOTPIX SLIDEMED\n")
        f.write("NARROWFIG\n")
        f.write("TIMEREF: 1\n")
    run_ir_all(temp_dir + "/analyze_cfg.txt")
    hot = np.loadtxt(f"{temp_dir}/analysis_hot.txt")
    (nh2, ncol) = np.shape(hot)
    assert 0.2 < nh2 / nh < 0.4
    assert ncol == 18
    nrec = 0
    for j in range(nh):
        nrec += np.count_nonzero(np.logical_and(x[j] == hot[:, 0], y[j] == hot[:, 1]))
    assert 9 * nh2 // 10 <= nrec <= nh2

    # see what SLIDEMED did
    diff = np.amax(np.abs(np.mean(data, axis=0) - expected_outputs) / tol)
    print(diff)
    assert diff < 1.0

    # linear fits
    _x = np.ones((nh2, 2))
    _x[:, 1] = hot[:, -5] + 4 * (hot[:, -4] + hot[:, -1])
    _y = hot[:, -4]
    c, m = np.linalg.lstsq(_x, _y)[0]
    assert -1 < c < 1
    assert 0.009 < m < 0.011
    _y = hot[:, -1]
    c, m = np.linalg.lstsq(_x, _y)[0]
    assert -1 < c < 1
    assert 0.0005 < m < 0.0015

    # LOGTSPACE
    with open(temp_dir + "/analyze_cfg.txt", "w") as f:
        f.write(analyze_cfg)
        f.write("HOTPIX: 400 2400 0.1 0.08\n")
        f.write("HOTPIX LOGTSPACE\n")
    run_ir_all(temp_dir + "/analyze_cfg.txt")
    hot = np.loadtxt(f"{temp_dir}/analysis_hot.txt")
    (nh2, ncol) = np.shape(hot)
    assert 0.2 < nh2 / nh < 0.4
    assert ncol == 63
    nrec = 0
    for j in range(nh):
        nrec += np.count_nonzero(np.logical_and(x[j] == hot[:, 0], y[j] == hot[:, 1]))
    assert 9 * nh2 // 10 <= nrec <= nh2
    sl = hot[0, 3::5]
    sl *= 18.0 / sl[-1]
    for k in range(7):
        assert k + 0.8 < sl[k] < k + 1.2
    assert 8.8 < sl[7] < 9.2
    assert 10.8 < sl[8] < 11.2
    assert 12.8 < sl[9] < 13.2
    assert 14.8 < sl[10] < 15.2

    # FULLNL
    analyze_cfg = (
        "DETECTOR: Test_simulation\n"
        "LIGHT:\n"
        f"    {temp_dir}/light_001.fits\n"
        f"    {temp_dir}/light_002.fits\n"
        f"    {temp_dir}/light_003.fits\n"
        f"    {temp_dir}/light_004.fits\n"
        "DARK:\n"
        f"    {temp_dir}/dark_005.fits\n"
        f"    {temp_dir}/dark_006.fits\n"
        f"    {temp_dir}/dark_007.fits\n"
        f"    {temp_dir}/dark_008.fits\n"
        "FORMAT: 1001\n"
        "CHAR: Advanced 1 3 3 bfe\n"
        "NBIN: 4 4\n"
        "TIME:    1 10 12 20\n"
        "NLPOLY: 3 1 20\n"
        "FULLNL: True True True\n"
        f"OUTPUT: {temp_dir}/analysis\n"
    )
    with open(temp_dir + "/analyze_cfg.txt", "w") as f:
        f.write(analyze_cfg)
    run_ir_all(temp_dir + "/analyze_cfg.txt")

    data = np.loadtxt(temp_dir + "/analysis_summary.txt")
    d = np.mean(data, axis=0)
    assert 1.49 < d[5] < 1.51
    assert 0.0098 < d[6] < 0.0102
    assert 0.0098 < d[7] < 0.0102
    assert -0.001 < d[-3] < 0.001
    assert -2.05 < d[-2] * 1.0e6 < -1.95
    assert -0.1 < d[-1] * 1.0e10 < 0.1
