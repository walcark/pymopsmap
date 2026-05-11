"""Unit tests for MicroParameters and DistrListPSD."""

import pytest
from pydantic import ValidationError

from pymopsmap.models.microparams import (
    DistrListPSD,
    DistrType,
    FixedPSD,
    MicroParameters,
    Sphere,
)


class TestScalarBroadcast:
    def test_scalar_n_real_broadcast(self):
        mp = MicroParameters(
            wavelength=[0.4, 0.55, 0.7],
            n_real=1.5,
            n_imag=0.001,
            shape=Sphere(),
            psd=FixedPSD(radius=0.1, n=1e6),
        )
        assert list(mp.n_real) == [1.5, 1.5, 1.5]
        assert list(mp.n_imag) == [0.001, 0.001, 0.001]

    def test_list_n_real_preserved(self):
        mp = MicroParameters(
            wavelength=[0.4, 0.55],
            n_real=[1.4, 1.5],
            n_imag=[0.001, 0.002],
            shape=Sphere(),
            psd=FixedPSD(radius=0.1, n=1e6),
        )
        assert list(mp.n_real) == [1.4, 1.5]

    def test_n_real_length_mismatch_raises(self):
        with pytest.raises(ValidationError):
            MicroParameters(
                wavelength=[0.4, 0.55],
                n_real=[1.4, 1.5, 1.6],
                n_imag=0.001,
                shape=Sphere(),
                psd=FixedPSD(radius=0.1, n=1e6),
            )


class TestDistrListPSD:
    def test_command_format(self):
        psd = DistrListPSD(
            radii=[0.1, 0.2, 0.5],
            concentrations=[1e4, 5e4, 1e3],
            distr_type=DistrType.DNDR,
        )
        cmd = psd.command
        assert cmd.startswith("size distr_list dndr")
        assert "0.1" in cmd
        assert (
            "1e+04" in cmd
            or "10000" in cmd
            or "1e4" in cmd
            or "10000.0" in cmd
        )

    def test_length_mismatch_raises(self):
        with pytest.raises(ValidationError):
            DistrListPSD(
                radii=[0.1, 0.2],
                concentrations=[1e4],
                distr_type=DistrType.DNDR,
            )

    def test_too_few_points_raises(self):
        with pytest.raises(ValidationError):
            DistrListPSD(
                radii=[0.1],
                concentrations=[1e4],
                distr_type=DistrType.DNDR,
            )

    def test_non_ascending_radii_raises(self):
        with pytest.raises(ValidationError):
            DistrListPSD(
                radii=[0.5, 0.1],
                concentrations=[1e4, 5e4],
                distr_type=DistrType.DNDR,
            )

    def test_distr_type_enum(self):
        for dt in DistrType:
            psd = DistrListPSD(
                radii=[0.1, 0.2],
                concentrations=[1e4, 5e4],
                distr_type=dt,
            )
            assert dt.value in psd.command

    def test_in_microparams(self):
        mp = MicroParameters(
            wavelength=[0.55],
            n_real=1.5,
            n_imag=0.01,
            shape=Sphere(),
            psd=DistrListPSD(
                radii=[0.05, 0.1, 0.2, 0.5],
                concentrations=[1e4, 5e4, 3e4, 1e4],
                distr_type=DistrType.DNDR,
            ),
        )
        assert isinstance(mp.psd, DistrListPSD)
