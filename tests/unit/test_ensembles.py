"""Turning a table of measured particles into the modes MOPSMAP runs."""

from __future__ import annotations

import numpy as np
import pytest

from pymopsmap.psd import FixedPSD, LognormalPSD
from pymopsmap.shapes import Sphere, SpheroidDistrFile
from pymopsmap.species import Mode


class TestFromParticles:
    def test_one_particle_on_a_grid_point_is_one_mode(self):
        modes = Mode.from_particles(
            radii=[0.5], aspect_ratios=[1.4], n_real=1.53, n_imag=0.0078
        )

        assert len(modes) == 1
        assert modes[0].shape.aspect_ratio == pytest.approx(1.4)
        assert modes[0].psd.radius == pytest.approx(0.5)
        assert modes[0].psd.n == pytest.approx(1.0)

    def test_a_ratio_between_two_points_is_split_between_them(self):
        """
        MOPSMAP would spread it the same way, over the eps grid it holds.

        Binning first is a compression, not an approximation: the
        contributions are the ones init_shape.f90 would build anyway.
        """
        modes = Mode.from_particles(
            radii=[0.5], aspect_ratios=[1.5], n_real=1.53, n_imag=0.0078
        )

        ratios = sorted(float(m.shape.aspect_ratio) for m in modes)
        assert ratios == pytest.approx([1.4, 1.6])
        assert sum(float(m.psd.n) for m in modes) == pytest.approx(1.0)
        assert [float(m.psd.n) for m in modes] == pytest.approx([0.5, 0.5])

    def test_particles_of_the_same_size_and_shape_are_one_mode(self):
        modes = Mode.from_particles(
            radii=[0.5, 0.5, 0.5],
            aspect_ratios=[1.4, 1.4, 2.0],
            n_real=1.53,
            n_imag=0.0078,
        )

        assert len(modes) == 2
        by_ratio = {float(m.shape.aspect_ratio): float(m.psd.n) for m in modes}
        assert by_ratio == pytest.approx({1.4: 2.0, 2.0: 1.0})

    def test_a_weight_per_particle_is_carried(self):
        modes = Mode.from_particles(
            radii=[0.5, 0.5],
            aspect_ratios=[1.4, 1.4],
            weights=[3.0, 1.0],
            n_real=1.53,
            n_imag=0.0078,
        )

        assert float(modes[0].psd.n) == pytest.approx(4.0)

    def test_sizes_and_ratios_are_clipped_where_the_data_set_stops(self):
        """
        Section 5.8 clips rather than drops: "particles with r > 47.5 um are
        modeled as r = 47.5 um and aspect ratios > 5 are set to 5".
        """
        modes = Mode.from_particles(
            radii=[80.0],
            aspect_ratios=[9.0],
            n_real=1.53,
            n_imag=0.0078,
            r_max=47.5,
        )

        assert float(modes[0].psd.radius) == pytest.approx(47.5)
        assert float(modes[0].shape.aspect_ratio) == pytest.approx(5.0)

    def test_a_ratio_of_one_stays_one(self):
        """MOPSMAP reads it from the sphere files; it is not an error."""
        modes = Mode.from_particles(
            radii=[0.5], aspect_ratios=[1.0], n_real=1.53, n_imag=0.0078
        )

        assert float(modes[0].shape.aspect_ratio) == pytest.approx(1.0)

    def test_a_ratio_below_one_is_refused(self):
        with pytest.raises(ValueError, match="at least one"):
            Mode.from_particles(
                radii=[0.5], aspect_ratios=[0.8], n_real=1.53, n_imag=0.0078
            )

    def test_the_lengths_must_match(self):
        with pytest.raises(ValueError, match="same length"):
            Mode.from_particles(
                radii=[0.5, 0.6],
                aspect_ratios=[1.4],
                n_real=1.53,
                n_imag=0.0078,
            )

    def test_every_mode_carries_the_microphysics_it_was_given(self):
        modes = Mode.from_particles(
            radii=[0.5, 1.0],
            aspect_ratios=[1.4, 2.0],
            n_real=1.53,
            n_imag=0.0078,
            density_dry=2.6,
            nonabs_fraction=0.5,
        )

        assert all(m.density_dry == 2.6 for m in modes)
        assert all(m.nonabs_fraction == 0.5 for m in modes)

    def test_the_total_weight_is_preserved_over_a_whole_table(self):
        rng = np.random.default_rng(0)
        count = 500
        modes = Mode.from_particles(
            radii=rng.uniform(0.1, 10.0, count).round(2),
            aspect_ratios=rng.uniform(1.0, 5.0, count),
            n_real=1.53,
            n_imag=0.0078,
        )

        assert sum(float(m.psd.n) for m in modes) == pytest.approx(count)


class TestFromIndexDistribution:
    def test_one_mode_per_bin_weighted(self):
        psd = LognormalPSD(rm=0.1, sigma=2.0, n=1e9, rmin=0.005, rmax=20.0)
        modes = Mode.from_index_distribution(
            n_imag=[0.001, 0.01, 0.1],
            weights=[0.5, 0.3, 0.2],
            n_real=1.53,
            shape=Sphere(),
            psd=psd,
        )

        assert [m.n_imag for m in modes] == pytest.approx([0.001, 0.01, 0.1])
        assert [float(m.psd.n) for m in modes] == pytest.approx(
            [5e8, 3e8, 2e8]
        )

    def test_the_weights_are_normalised(self):
        psd = LognormalPSD(rm=0.1, sigma=2.0, n=1e9, rmin=0.005, rmax=20.0)
        modes = Mode.from_index_distribution(
            n_imag=[0.001, 0.01],
            weights=[30.0, 10.0],
            n_real=1.53,
            shape=Sphere(),
            psd=psd,
        )

        assert sum(float(m.psd.n) for m in modes) == pytest.approx(1e9)

    def test_an_empty_bin_produces_no_mode(self):
        psd = FixedPSD(radius=0.5, n=1e9)
        modes = Mode.from_index_distribution(
            n_imag=[0.001, 0.01, 0.1],
            weights=[0.5, 0.0, 0.5],
            n_real=1.53,
            shape=Sphere(),
            psd=psd,
        )

        assert len(modes) == 2

    def test_it_works_on_a_distribution_file_shape(self):
        psd = LognormalPSD(rm=0.1, sigma=2.0, n=1e9, rmin=0.005, rmax=20.0)
        modes = Mode.from_index_distribution(
            n_imag=[0.001, 0.01],
            weights=[0.5, 0.5],
            n_real=1.53,
            shape=SpheroidDistrFile(distr_filename="ar_kandler"),
            psd=psd,
        )

        assert all(m.shape.type == "spheroid-distr-file" for m in modes)

    def test_the_lengths_must_match(self):
        psd = FixedPSD(radius=0.5, n=1.0)
        with pytest.raises(ValueError, match="same length"):
            Mode.from_index_distribution(
                n_imag=[0.001, 0.01],
                weights=[1.0],
                n_real=1.53,
                shape=Sphere(),
                psd=psd,
            )
