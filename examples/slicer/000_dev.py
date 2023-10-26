import xtrack as xt
import xpart as xp
import xobjects as xo
import xfields as xf

from pathlib import Path

import numpy as np

# line = xt.Line.from_json(
#     '../../../xtrack/test_data/sps_w_spacecharge/line_no_spacecharge_and_particle.json')
# line.particle_ref = xt.Particles(p0c=26e9, mass0=xt.PROTON_MASS_EV)
# line.build_tracker()
# tw = line.twiss()

# num_partilces_per_bunch = 100
# num_bunches = 3
# total_intensity_particles_bunch = 1e11

# beam = xp.generate_matched_gaussian_bunch(
#             num_particles=num_partilces_per_bunch * num_bunches,
#             total_intensity_particles=total_intensity_particles_bunch * num_bunches,
#             sigma_z=0.1, nemitt_x=2.5e-6, nemitt_y=2.5e-6, line=line)

# harmonic_number = 4620
# dz_bucket = tw.circumference / harmonic_number
# bunch_spacing_buckets = 5

# for ii in range(num_bunches):
#     beam.zeta[ii * num_partilces_per_bunch:(ii+1) * num_partilces_per_bunch] += (
#         ii * bunch_spacing_buckets * dz_bucket)


_configure_grid = xf.fieldmaps.interpolated._configure_grid


class UniformBinSlicer(xt.BeamElement):
    _xofields = {
        'z_min': xo.Float64,
        'num_slices': xo.Int64,
        'dzeta': xo.Float64,
        'i_bunch_0': xo.Int64,
        'num_bunches': xo.Int64,
        'bunch_spacing_zeta': xo.Float64,
        'particles_per_slice': xo.Float64[:],
    }

    _rename = {
        'z_min': '_z_min',
        'num_slices': '_num_slices',
        'dzeta': '_dzeta',
        'i_bunch_0': '_i_bunch_0',
        'num_bunches': '_num_bunches',
        'bunch_spacing_zeta': '_bunch_spacing_zeta',
        'particles_per_slice': '_particles_per_slice',
    }

    _extra_c_sources = [
        Path('uniform_bin_slicer.h')
    ]

    _per_particle_kernels = {
            'test_slice': xo.Kernel(
                c_name='UniformBinSlicer_slice',
                args=[
                    xo.Arg(xo.Int64, pointer=True, name='i_slice_for_particles'),
                    xo.Arg(xo.Int64, pointer=True, name='i_bunch_for_particles')
                ]),
        }

    def __init__(self, zeta_range=None, nbins=None, dzeta=None, zeta_grid=None,
                 num_bunches=None, i_bunch_0=None, bunch_spacing_zeta=None,
                 **kwargs):

        self._zeta_grid = _configure_grid('zeta', zeta_grid, dzeta, zeta_range, nbins)
        num_bunches = num_bunches or 0
        i_bunch_0 = i_bunch_0 or 0
        bunch_spacing_zeta = bunch_spacing_zeta or 0

        self.xoinitialize(z_min=self.zeta_grid[0], num_slices=self.num_slices,
                          dzeta=self.dzeta,
                          num_bunches=num_bunches, i_bunch_0=i_bunch_0,
                          bunch_spacing_zeta=bunch_spacing_zeta,
                          particles_per_slice=(num_bunches or 1) * self.num_slices, # initialization with tuple not working
                          **kwargs)
    @property
    def zeta_grid(self):
        """
        Array with the grid points (bin centers).
        """
        return self._zeta_grid


    @property
    def num_slices(self):
        """
        Number of bins
        """
        return len(self.zeta_grid)

    @property
    def dzeta(self):
        """
        Bin size in meters.
        """
        return self.zeta_grid[1] - self.zeta_grid[0]

    @property
    def num_bunches(self):
        """
        Number of bunches
        """
        return self._num_bunches

    @property
    def i_bunch_0(self):
        """
        Index of the first bunch
        """
        return self._i_bunch_0

    @property
    def particles_per_slice(self):
        """
        Number of particles per slice
        """
        return self._particles_per_slice.reshape(self.num_bunches, self.num_slices)

# Check in single-bunch mode

slicer = UniformBinSlicer(zeta_range=(-1, 1), nbins=3)
assert slicer.num_bunches == 0 # Single-bunch mode

p0 = xt.Particles(zeta  =[-2, -1.51, -1.49, -1, -0.51, -0.49, 0, 0.49, 0.51,  1, 1.49, 1.51,  2, 2.51],
                  weight=[10,  10,    10,    10, 10,    20,    20, 20,   30,  30,  30,   40, 40,   40])
p0.state[-1] = 0

p = p0.copy()

i_slice_expected    = [-1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1, -999]

ss = 0 * p.x
ctx= xo.ContextCpu()
i_slice_for_particles = p.particle_id * 0 - 999
i_bunch_for_particles = p.particle_id * 0 - 999
slicer.test_slice(particles=p, i_slice_for_particles=i_slice_for_particles,
                  i_bunch_for_particles=i_bunch_for_particles)

assert np.all(np.array(i_slice_expected) == i_slice_for_particles)

# Check in multi-bunch mode
bunch_spacing_zeta = 10.

p1 = p0.copy()
p2 = p0.copy()
p2.zeta += bunch_spacing_zeta
p2.weight *= 10
p3 = p0.copy()
p3.zeta += 2 * bunch_spacing_zeta
p3.weight *= 100
p4 = p0.copy()
p4.zeta += 3 * bunch_spacing_zeta
p4.weight *= 1000

p = xp.Particles.merge([p1, p2, p3, p4])

i_bunch_for_particles = p.particle_id * 0 - 999
i_slice_for_particles = p.particle_id * 0 - 999

slicer = UniformBinSlicer(zeta_range=(-1, 1), nbins=3, i_bunch_0=0,
                          num_bunches=4, bunch_spacing_zeta=bunch_spacing_zeta)
slicer.test_slice(particles=p, i_slice_for_particles=i_slice_for_particles,
                    i_bunch_for_particles=i_bunch_for_particles)

i_slice_expected  = np.array([
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -999, -999, -999, -999
])
i_bunch_expected  = np.array([
    -1, -1,    0,      0,  0,    0,     0,    0,    0, 0, 0,     0,   0,
     0,  0,    1,      1,  1,    1,     1,    1,    1, 1, 1,     1,  1,
     1,  1,    2,      2,  2,    2,     2,    2,    2, 2, 2,     2,  2,
     2,  2,    3,      3,  3,    3,     3,    3,    3, 3, 3,     3,  3,
    -999, -999, -999, -999
])

assert np.all(i_slice_for_particles == i_slice_expected)
assert np.all(i_bunch_for_particles == i_bunch_expected)
