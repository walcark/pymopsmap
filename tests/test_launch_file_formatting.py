import textwrap
import re

from pymopsmap.classes import MicroParameters, Sphere, FixedPSD
from pymopsmap.mopsmap.launch_file_format import write_launching_file


def test_mopsmap_launch_file_format():
    def replace_temp_file(content: str) -> str:
        tmp_pattern = r"/tmp/pymopsmap-[^/]+"
        new_line = "/tmp/pymopsmap-TMP"
        new_content = re.sub(tmp_pattern, new_line, content)
        new_content = "\n".join(
            line.rstrip() for line in new_content.splitlines()
        ).strip()
        return new_content

    mp: MicroParameters = MicroParameters(
        wavelength=[0.500, 0.700, 0.900, 1.5],
        n_real=[1.00027] * 100,
        n_imag=[0.0] * 100,
        shape=Sphere(),
        psd=FixedPSD(radius=0.1, n=1.0),
    )

    # With rh
    path_dict = write_launching_file(mp, n_angles=2000, rh=0.0)
    file_path = path_dict["mopsmap"]

    with open(file_path, "r") as f:
        content = replace_temp_file(f.read())

    expected_content = textwrap.dedent(
        """
        scatlib '/home/kwalcarius/bin/mopsmap/optical_dataset'
        mode 1 shape sphere
        mode 1 size 0.1 1.0
        mode 1 refrac file '/tmp/pymopsmap-9870a6e0d3e04e3fb6c888408c0d53d1/ri.txt'
        output num_theta 2000
        wavelength from_refrac_file
        rH 0.0
        output integrated
        output netcdf '/tmp/pymopsmap-9870a6e0d3e04e3fb6c888408c0d53d1/output.nc'
        """
    )

    assert content == replace_temp_file(expected_content)

    # Without rh
    path_dict = write_launching_file(mp, n_angles=2000)
    file_path = path_dict["mopsmap"]

    with open(file_path, "r") as f:
        content = replace_temp_file(f.read())

    expected_content = textwrap.dedent(
        """
        scatlib '/home/kwalcarius/bin/mopsmap/optical_dataset'
        mode 1 shape sphere
        mode 1 size 0.1 1.0
        mode 1 refrac file '/tmp/pymopsmap-9870a6e0d3e04e3fb6c888408c0d53d1/ri.txt'
        output num_theta 2000
        wavelength from_refrac_file
        output integrated
        output netcdf '/tmp/pymopsmap-9870a6e0d3e04e3fb6c888408c0d53d1/output.nc'
        """
    )

    assert content == replace_temp_file(expected_content)
