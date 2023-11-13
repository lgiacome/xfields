import xtrack as xt
import xpart as xp
import xobjects as xo

from xfields import UniformBinSlicer

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


###############################################
# Check slice attribution (single-bunch mode) #
###############################################

slicer = UniformBinSlicer(zeta_range=(-1, 1), num_slices=3)
assert slicer.num_bunches == 0 # Single-bunch mode

p0 = xt.Particles(zeta  =[-2, -1.51, -1.49, -1, -0.51, -0.49, 0, 0.49, 0.51,  1, 1.49, 1.51,  2, 2.51],
                  weight=[10,  10,    10,    10, 10,    20,    20,  20,   30,  30,  30,   40, 40,   40],
                  x=     [0.,  1.,   2.,     3.,  4.,    5.,    6., 7.,    8.,  9.,  10.,  11., 12.,  13.],
                  y=     [13., 12.,  11.,    10., 9.,    8.,    7., 6.,    5.,  4.,  3.,   2.,  1.,   0.]
                  )
p0.state[-1] = 0

p = p0.copy()

i_slice_expected    = [-1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1, -999]

ss = 0 * p.x
ctx= xo.ContextCpu()
i_slice_particles = p.particle_id * 0 - 999
i_bunch_particles = p.particle_id * 0 - 9999
slicer.slice(particles=p, i_bunch_particles=i_bunch_particles,
                    i_slice_particles=i_slice_particles)

assert np.all(np.array(i_slice_expected) == i_slice_particles)
assert np.all(i_bunch_particles == -9999)


expected_particles_per_slice = np.array([30, 60, 90])
assert np.allclose(slicer.particles_per_slice, expected_particles_per_slice,
                     atol=1e-12, rtol=0)

##############################################
# Check slice attribution (multi-bunch mode) #
##############################################

bunch_spacing_zeta = 10.

p1 = p0.copy()
p2 = p0.copy()
p2.zeta -= bunch_spacing_zeta
p2.weight *= 10
p3 = p0.copy()
p3.zeta -= 2 * bunch_spacing_zeta
p3.weight *= 100
p4 = p0.copy()
p4.zeta -= 3 * bunch_spacing_zeta
p4.weight *= 1000

p = xt.Particles.merge([p1, p2, p3, p4])

i_bunch_particles = p.particle_id * 0 - 999
i_slice_particles = p.particle_id * 0 - 999

slicer = UniformBinSlicer(zeta_range=(-1, 1), num_slices=3, i_bunch_0=0,
                          num_bunches=4, bunch_spacing_zeta=bunch_spacing_zeta)
slicer.slice(particles=p,
                    i_slice_particles=i_slice_particles,
                    i_bunch_particles=i_bunch_particles)

i_slice_expected  = np.array([
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -1, -1,    0,      0,  0,    1,     1,    1,    2, 2, 2,    -1,  -1,
    -999, -999, -999, -999
])
# i_bunch_expected  = np.array([
#     -1, -1,    0,      0,  0,    0,     0,    0,    0, 0, 0,     0,   0,
#      0,  0,    1,      1,  1,    1,     1,    1,    1, 1, 1,     1,   1,
#      1,  1,    2,      2,  2,    2,     2,    2,    2, 2, 2,     2,   2,
#      2,  2,    3,      3,  3,    3,     3,    3,    3, 3, 3,     3,   3,
#     -999, -999, -999, -999
# ])

i_bunch_expected  = np.array([
     1,  1,    0,      0,  0,    0,     0,    0,    0, 0, 0,     0,   0,
     2,  2,    1,      1,  1,    1,     1,    1,    1, 1, 1,     1,   1,
     3,  3,    2,      2,  2,    2,     2,    2,    2, 2, 2,     2,   2,
     -1,  -1,    3,      3,  3,    3,     3,    3,    3, 3, 3,     3,   3,
    -999, -999, -999, -999
])

expected_particles_per_slice = np.array([
    [30, 60, 90],
    [300, 600, 900],
    [3000, 6000, 9000],
    [30000, 60000, 90000],
])

assert np.all(i_slice_particles == i_slice_expected)
assert np.all(i_bunch_particles == i_bunch_expected)
assert np.allclose(slicer.particles_per_slice, expected_particles_per_slice,
                   atol=1e-12, rtol=0)


#####################################
# Check moments (single-bunch mode) #
#####################################

slicer_single_bunch = UniformBinSlicer(zeta_range=(-1, 1), num_slices=3)

p = xt.Particles(zeta=[0.99, 1.0, 1.01],
                 weight=[1, 2, 1],
                 x = [99, 100, 101],
                 y = [201,200, 199])
slicer_single_bunch.slice(p)

assert slicer_single_bunch.bunch_spacing_zeta == 0

assert np.allclose(slicer_single_bunch.zeta_centers, np.array([-1, 0, 1]), rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.particles_per_slice, [0, 0, p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('x'), [0, 0, (p.x * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('y'), [0, 0, (p.y * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('zeta'), [0, 0, (p.zeta * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('xx'), [0, 0, (p.x**2 * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('yy'), [0, 0, (p.y**2 * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('zetazeta'), [0, 0, (p.zeta**2 * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('xy'), [0, 0, (p.x * p.y * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.sum('xzeta'), [0, 0, (p.x * p.zeta * p.weight).sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('x'), [0, 0, (p.x * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('y'), [0, 0, (p.y * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('xx'), [0, 0, (p.x**2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('yy'), [0, 0, (p.y**2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('zetazeta'), [0, 0, (p.zeta**2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('xy'), [0, 0, (p.x * p.y * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.mean('xzeta'), [0, 0, (p.x * p.zeta * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.cov('x', 'y'),
    slicer_single_bunch.mean('xy') - slicer_single_bunch.mean('x') * slicer_single_bunch.mean('y'),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.var('x'), slicer_single_bunch.cov('x', 'x'),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.var('x'), slicer_single_bunch.mean('xx') - slicer_single_bunch.mean('x')**2,
    rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.var('zeta'), slicer_single_bunch.mean('zetazeta') - slicer_single_bunch.mean('zeta')**2,
    rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.std('x'), np.sqrt(slicer_single_bunch.var('x')),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.std('y'), np.sqrt(slicer_single_bunch.var('y')),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_single_bunch.std('zeta'), np.sqrt(slicer_single_bunch.var('zeta')),
    rtol=0, atol=1e-12)

assert np.all(slicer_single_bunch.sum('xy') == slicer_single_bunch.sum('x_y'))
assert np.all(slicer_single_bunch.sum('x', 'y') == slicer_single_bunch.sum('x_y'))
assert np.all(slicer_single_bunch.mean('xy') == slicer_single_bunch.mean('x_y'))
assert np.all(slicer_single_bunch.mean('x', 'y') == slicer_single_bunch.mean('x_y'))
assert np.all(slicer_single_bunch.cov('xy') == slicer_single_bunch.cov('x_y'))
assert np.all(slicer_single_bunch.cov('x', 'y') == slicer_single_bunch.cov('x_y'))

# # Same parametrized

moms = ['x_px', 'x_y', 'x_py', 'x_delta',
        'px_y', 'px_py', 'px_delta',
        'y_py', 'y_delta',
        'py_delta']
for mm in moms:
    c1_name, c2_name = mm.split('_')

    p = xt.Particles(zeta=[0.99, 1.0, 1.01],
                    weight=[1, 2, 1],
                    **{c1_name: [99, 100, 101],
                       c2_name: [201,200, 199]})
    c1 = getattr(p, c1_name)
    c2 = getattr(p, c2_name)

    slicer_single_bunch.slice(p)

    assert np.allclose(slicer_single_bunch.zeta_centers, np.array([-1, 0, 1]), rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.particles_per_slice, [0, 0, p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum(c1_name), [0, 0, (c1 * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum(c2_name), [0, 0, (c2 * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum('zeta'), [0, 0, (p.zeta * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum(c1_name + c1_name), [0, 0, (c1**2 * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum(c2_name + c2_name), [0, 0, (c2**2 * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum(c1_name + c2_name), [0, 0, (c1 * c2 * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum('zetazeta'), [0, 0, (p.zeta**2 * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.sum(c1_name + 'zeta'), [0, 0, (c1 * p.zeta * p.weight).sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.mean(c1_name), [0, 0, (c1 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.mean(c2_name), [0, 0, (c2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.mean(c1_name + c1_name), [0, 0, (c1**2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.mean(c2_name + c2_name), [0, 0, (c2**2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.mean(c1_name + c2_name), [0, 0, (c1 * c2 * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.mean(c1_name + 'zeta'), [0, 0, (c1 * p.zeta * p.weight).sum() / p.weight.sum()], rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.cov(c1_name, c2_name),
        slicer_single_bunch.mean(c1_name + c2_name) - slicer_single_bunch.mean(c1_name) * slicer_single_bunch.mean(c2_name),
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.var(c1_name), slicer_single_bunch.cov(c1_name, c1_name),
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.var(c1_name), slicer_single_bunch.mean(c1_name + c1_name) - slicer_single_bunch.mean(c1_name)**2,
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_single_bunch.var('zeta'), slicer_single_bunch.mean('zetazeta') - slicer_single_bunch.mean('zeta')**2,
        rtol=0, atol=1e-12)

    assert np.all(slicer_single_bunch.sum(c1_name + c2_name) == slicer_single_bunch.sum(c1_name + '_' + c2_name))
    assert np.all(slicer_single_bunch.sum(c1_name, c2_name) == slicer_single_bunch.sum(c1_name + '_' + c2_name))
    assert np.all(slicer_single_bunch.mean(c1_name + c2_name) == slicer_single_bunch.mean(c1_name + '_' + c2_name))
    assert np.all(slicer_single_bunch.mean(c1_name, c2_name) == slicer_single_bunch.mean(c1_name + '_' + c2_name))
    assert np.all(slicer_single_bunch.cov(c1_name + c2_name) == slicer_single_bunch.cov(c1_name + '_' + c2_name))
    assert np.all(slicer_single_bunch.cov(c1_name, c2_name) == slicer_single_bunch.cov(c1_name + '_' + c2_name))

####################################
# Check moments (multi-bunch mode) #
####################################

slicer_multi_bunch = UniformBinSlicer(zeta_range=(-1, 1), num_slices=3,
                                        num_bunches=4, bunch_spacing_zeta=bunch_spacing_zeta)

slicer_multi_bunch_part = UniformBinSlicer(zeta_range=(-1, 1), num_slices=3, i_bunch_0=1,
                                        num_bunches=3, bunch_spacing_zeta=bunch_spacing_zeta)

p1 = xt.Particles(zeta=[0.99, 1.0, 1.01],
                 weight=[1, 2, 1],
                 x = [99, 100, 101],
                 y = [201,200, 199])
p2 = xt.Particles(zeta=np.array([-0.01, 0, 0.01]) - 2 * bunch_spacing_zeta,
                    weight=[1, 2, 1],
                    x = [99, 100, 101],
                    y = [201,200, 199])
p = xt.Particles.merge([p1, p2])

assert np.isclose(slicer_multi_bunch.bunch_spacing_zeta, 10, rtol=0, atol=1e-12)
assert np.isclose(slicer_multi_bunch_part.bunch_spacing_zeta, 10, rtol=0, atol=1e-12)

assert slicer_multi_bunch.num_bunches == 4
assert slicer_multi_bunch_part.num_bunches == 3
assert slicer_multi_bunch.i_bunch_0 == 0
assert slicer_multi_bunch_part.i_bunch_0 == 1

slicer_multi_bunch.slice(p)
slicer_multi_bunch_part.slice(p)

assert np.allclose(slicer_multi_bunch.zeta_centers, np.array([[-1, 0, 1], [-11, -10, -9], [-21, -20, -19], [-31, -30, -29]]), rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.particles_per_slice, [[0, 0, p1.weight.sum()], [0, 0, 0], [0, p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('x'), [[0, 0, (p1.x * p1.weight).sum()], [0, 0, 0], [0, (p2.x * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('y'), [[0, 0, (p1.y * p1.weight).sum()], [0, 0, 0], [0, (p2.y * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('zeta'), [[0, 0, (p1.zeta * p1.weight).sum()], [0, 0, 0], [0, (p2.zeta * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('xx'), [[0, 0, (p1.x**2 * p1.weight).sum()], [0, 0, 0], [0, (p2.x**2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('yy'), [[0, 0, (p1.y**2 * p1.weight).sum()], [0, 0, 0], [0, (p2.y**2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('zetazeta'), [[0, 0, (p1.zeta**2 * p1.weight).sum()], [0, 0, 0], [0, (p2.zeta**2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('xy'), [[0, 0, (p1.x * p1.y * p1.weight).sum()], [0, 0, 0], [0, (p2.x * p2.y * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.sum('xzeta'), [[0, 0, (p1.x * p1.zeta * p1.weight).sum()], [0, 0, 0], [0, (p2.x * p2.zeta * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.mean('x'), [[0, 0, (p1.x * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (p2.x * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.mean('y'), [[0, 0, (p1.y * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (p2.y * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.mean('xx'), [[0, 0, (p1.x**2 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (p2.x**2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.mean('yy'), [[0, 0, (p1.y**2 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (p2.y**2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.mean('xy'), [[0, 0, (p1.x * p1.y * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (p2.x * p2.y * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.mean('xzeta'), [[0, 0, (p1.x * p1.zeta * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (p2.x * p2.zeta * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.cov('x', 'y'),
    slicer_multi_bunch.mean('xy') - slicer_multi_bunch.mean('x') * slicer_multi_bunch.mean('y'),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.var('x'), slicer_multi_bunch.cov('x', 'x'),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.var('x'), slicer_multi_bunch.mean('xx') - slicer_multi_bunch.mean('x')**2,
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.var('zeta'), slicer_multi_bunch.mean('zetazeta') - slicer_multi_bunch.mean('zeta')**2,
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.std('x'), np.sqrt(slicer_multi_bunch.var('x')),
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.std('y'), np.sqrt(slicer_multi_bunch.var('y')),
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch.std('zeta'), np.sqrt(slicer_multi_bunch.var('zeta')),
                    rtol=0, atol=1e-12)

assert np.all(slicer_multi_bunch.sum('xy') == slicer_multi_bunch.sum('x_y'))
assert np.all(slicer_multi_bunch.sum('x', 'y') == slicer_multi_bunch.sum('x_y'))
assert np.all(slicer_multi_bunch.mean('xy') == slicer_multi_bunch.mean('x_y'))
assert np.all(slicer_multi_bunch.mean('x', 'y') == slicer_multi_bunch.mean('x_y'))
assert np.all(slicer_multi_bunch.cov('xy') == slicer_multi_bunch.cov('x_y'))
assert np.all(slicer_multi_bunch.cov('x', 'y') == slicer_multi_bunch.cov('x_y'))

# Check slicer_part
assert np.allclose(slicer_multi_bunch_part.zeta_centers, slicer_multi_bunch.zeta_centers[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.particles_per_slice, slicer_multi_bunch.particles_per_slice[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('x'), slicer_multi_bunch.sum('x')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('y'), slicer_multi_bunch.sum('y')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('zeta'), slicer_multi_bunch.sum('zeta')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('xx'), slicer_multi_bunch.sum('xx')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('yy'), slicer_multi_bunch.sum('yy')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('zetazeta'), slicer_multi_bunch.sum('zetazeta')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('xy'), slicer_multi_bunch.sum('xy')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.sum('xzeta'), slicer_multi_bunch.sum('xzeta')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.mean('x'), slicer_multi_bunch.mean('x')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.mean('y'), slicer_multi_bunch.mean('y')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.mean('xx'), slicer_multi_bunch.mean('xx')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.mean('yy'), slicer_multi_bunch.mean('yy')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.mean('xy'), slicer_multi_bunch.mean('xy')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.mean('xzeta'), slicer_multi_bunch.mean('xzeta')[1:],
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.cov('x', 'y'),
    slicer_multi_bunch_part.mean('xy') - slicer_multi_bunch_part.mean('x') * slicer_multi_bunch_part.mean('y'),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.var('x'), slicer_multi_bunch_part.cov('x', 'x'),
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.var('x'), slicer_multi_bunch_part.mean('xx') - slicer_multi_bunch_part.mean('x')**2,
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.var('zeta'), slicer_multi_bunch_part.mean('zetazeta') - slicer_multi_bunch_part.mean('zeta')**2,
    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.std('x'), np.sqrt(slicer_multi_bunch_part.var('x')),
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.std('y'), np.sqrt(slicer_multi_bunch_part.var('y')),
                    rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_part.std('zeta'), np.sqrt(slicer_multi_bunch_part.var('zeta')),
                    rtol=0, atol=1e-12)


# # Same parametrized

moms = ['x_px', 'x_y', 'x_py', 'x_delta',
        'px_y', 'px_py', 'px_delta',
        'y_py', 'y_delta',
        'py_delta']

for mm in moms:
    c1_name, c2_name = mm.split('_')

    p1 = xt.Particles(zeta=[0.99, 1.0, 1.01],
                    weight=[1, 2, 1],
                    **{c1_name: [99, 100, 101],
                       c2_name: [201,200, 199]})
    p2 = xt.Particles(zeta=np.array([-0.01, 0, 0.01]) - 2 * bunch_spacing_zeta,
                        weight=[1, 2, 1],
                        **{c1_name: [99, 100, 101],
                           c2_name: [201,200, 199]})

    p = xt.Particles.merge([p1, p2])

    slicer_multi_bunch.slice(p)
    slicer_multi_bunch_part.slice(p)

    c1_p1 = getattr(p1, c1_name)
    c2_p1 = getattr(p1, c2_name)
    c1_p2 = getattr(p2, c1_name)
    c2_p2 = getattr(p2, c2_name)

    assert np.allclose(slicer_multi_bunch.zeta_centers, np.array([[-1, 0, 1], [-11, -10, -9], [-21, -20, -19], [-31, -30, -29]]), rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.particles_per_slice, [[0, 0, p1.weight.sum()], [0, 0, 0], [0, p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum(c1_name), [[0, 0, (c1_p1 * p1.weight).sum()], [0, 0, 0], [0, (c1_p2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum(c2_name), [[0, 0, (c2_p1 * p1.weight).sum()], [0, 0, 0], [0, (c2_p2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum('zeta'), [[0, 0, (p1.zeta * p1.weight).sum()], [0, 0, 0], [0, (p2.zeta * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum(c1_name + c1_name), [[0, 0, (c1_p1**2 * p1.weight).sum()], [0, 0, 0], [0, (c1_p2**2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum(c2_name + c2_name), [[0, 0, (c2_p1**2 * p1.weight).sum()], [0, 0, 0], [0, (c2_p2**2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum(c1_name + c2_name), [[0, 0, (c1_p1 * c2_p1 * p1.weight).sum()], [0, 0, 0], [0, (c1_p2 * c2_p2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum('zetazeta'), [[0, 0, (p1.zeta**2 * p1.weight).sum()], [0, 0, 0], [0, (p2.zeta**2 * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.sum(c1_name + 'zeta'), [[0, 0, (c1_p1 * p1.zeta * p1.weight).sum()], [0, 0, 0], [0, (c1_p2 * p2.zeta * p2.weight).sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.mean(c1_name), [[0, 0, (c1_p1 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (c1_p2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.mean(c2_name), [[0, 0, (c2_p1 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (c2_p2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.mean(c1_name + c1_name), [[0, 0, (c1_p1**2 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (c1_p2**2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.mean(c2_name + c2_name), [[0, 0, (c2_p1**2 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (c2_p2**2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.mean(c1_name + c2_name), [[0, 0, (c1_p1 * c2_p1 * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (c1_p2 * c2_p2 * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.mean(c1_name + 'zeta'), [[0, 0, (c1_p1 * p1.zeta * p1.weight).sum() / p1.weight.sum()], [0, 0, 0], [0, (c1_p2 * p2.zeta * p2.weight).sum() / p2.weight.sum(), 0], [0, 0, 0]], rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.cov(c1_name, c2_name),
        slicer_multi_bunch.mean(c1_name + c2_name) - slicer_multi_bunch.mean(c1_name) * slicer_multi_bunch.mean(c2_name),
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.var(c1_name), slicer_multi_bunch.cov(c1_name, c1_name),
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.var(c1_name), slicer_multi_bunch.mean(c1_name + c1_name) - slicer_multi_bunch.mean(c1_name)**2,
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.var('zeta'), slicer_multi_bunch.mean('zetazeta') - slicer_multi_bunch.mean('zeta')**2,
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.std(c1_name), np.sqrt(slicer_multi_bunch.var(c1_name)),
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.std(c2_name), np.sqrt(slicer_multi_bunch.var(c2_name)),
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch.std('zeta'), np.sqrt(slicer_multi_bunch.var('zeta')),
                        rtol=0, atol=1e-12)

    assert np.all(slicer_multi_bunch.sum(c1_name + c2_name) == slicer_multi_bunch.sum(c1_name + '_' + c2_name))
    assert np.all(slicer_multi_bunch.sum(c1_name, c2_name) == slicer_multi_bunch.sum(c1_name + '_' + c2_name))
    assert np.all(slicer_multi_bunch.mean(c1_name + c2_name) == slicer_multi_bunch.mean(c1_name + '_' + c2_name))
    assert np.all(slicer_multi_bunch.mean(c1_name, c2_name) == slicer_multi_bunch.mean(c1_name + '_' + c2_name))
    assert np.all(slicer_multi_bunch.cov(c1_name + c2_name) == slicer_multi_bunch.cov(c1_name + '_' + c2_name))
    assert np.all(slicer_multi_bunch.cov(c1_name, c2_name) == slicer_multi_bunch.cov(c1_name + '_' + c2_name))

    # Check slicer_part
    assert np.allclose(slicer_multi_bunch_part.zeta_centers, slicer_multi_bunch.zeta_centers[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.particles_per_slice, slicer_multi_bunch.particles_per_slice[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum(c1_name), slicer_multi_bunch.sum(c1_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum(c2_name), slicer_multi_bunch.sum(c2_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum('zeta'), slicer_multi_bunch.sum('zeta')[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum(c1_name + c1_name), slicer_multi_bunch.sum(c1_name + c1_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum(c2_name + c2_name), slicer_multi_bunch.sum(c2_name + c2_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum(c1_name + c2_name), slicer_multi_bunch.sum(c1_name + c2_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum('zetazeta'), slicer_multi_bunch.sum('zetazeta')[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.sum(c1_name + 'zeta'), slicer_multi_bunch.sum(c1_name + 'zeta')[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.mean(c1_name), slicer_multi_bunch.mean(c1_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.mean(c2_name), slicer_multi_bunch.mean(c2_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.mean(c1_name + c1_name), slicer_multi_bunch.mean(c1_name + c1_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.mean(c2_name + c2_name), slicer_multi_bunch.mean(c2_name + c2_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.mean(c1_name + c2_name), slicer_multi_bunch.mean(c1_name + c2_name)[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.mean(c1_name + 'zeta'), slicer_multi_bunch.mean(c1_name + 'zeta')[1:],
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.cov(c1_name, c2_name),
        slicer_multi_bunch_part.mean(c1_name + c2_name) - slicer_multi_bunch_part.mean(c1_name) * slicer_multi_bunch_part.mean(c2_name),
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.var(c1_name), slicer_multi_bunch_part.cov(c1_name, c1_name),
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.var(c1_name), slicer_multi_bunch_part.mean(c1_name + c1_name) - slicer_multi_bunch_part.mean(c1_name)**2,
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.var('zeta'), slicer_multi_bunch_part.mean('zetazeta') - slicer_multi_bunch_part.mean('zeta')**2,
        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.std(c1_name), np.sqrt(slicer_multi_bunch_part.var(c1_name)),
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.std(c2_name), np.sqrt(slicer_multi_bunch_part.var(c2_name)),
                        rtol=0, atol=1e-12)
    assert np.allclose(slicer_multi_bunch_part.std('zeta'), np.sqrt(slicer_multi_bunch_part.var('zeta')),
                        rtol=0, atol=1e-12)

# Try selected moments
slicer_multi_bunch_mom = UniformBinSlicer(
    zeta_range=(-1, 1), num_slices=3,
    num_bunches=4, bunch_spacing_zeta=bunch_spacing_zeta,
    moments=['delta', 'xy', 'px_px'])

assert np.all(np.array(slicer_multi_bunch_mom.moments) == np.array(
    ['x', 'px', 'y', 'delta', 'x_y', 'px_px']))

p1 = xt.Particles(zeta=[0.99, 1.0, 1.01],
                 weight=[1, 2, 1],
                 x = [99, 100, 101],
                 y = [201,200, 199])
p2 = xt.Particles(zeta=np.array([-0.01, 0, 0.01]) + 2 * bunch_spacing_zeta,
                    weight=[1, 2, 1],
                    x = [99, 100, 101],
                    y = [201,200, 199])
p = xt.Particles.merge([p1, p2])

slicer_multi_bunch_mom.slice(p)
slicer_multi_bunch.slice(p)

assert np.allclose(slicer_multi_bunch_mom.particles_per_slice,
                     slicer_multi_bunch.particles_per_slice,
                     rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.zeta_centers,
                        slicer_multi_bunch.zeta_centers,
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.sum('x'),
                        slicer_multi_bunch.sum('x'),
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.sum('y'),
                        slicer_multi_bunch.sum('y'),
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.sum('px'),
                        slicer_multi_bunch.sum('px'),
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.sum('delta'),
                        slicer_multi_bunch.sum('delta'),
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.cov('x_y'),
                        slicer_multi_bunch.cov('x_y'),
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.cov('px_px'),
                        slicer_multi_bunch.cov('px_px'),
                        rtol=0, atol=1e-12)
assert np.allclose(slicer_multi_bunch_mom.var('px'),
                        slicer_multi_bunch.var('px'),
                        rtol=0, atol=1e-12)

p = xt.Particles(zeta=np.random.uniform(-1, 1, int(1e6)),
                 x = np.random.normal(0, 1, int(1e6)))
slicer_time = UniformBinSlicer(zeta_range=(-1, 1), num_slices=100,
                               moments=['x', 'xx'])