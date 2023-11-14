import matplotlib.pyplot as plt
import numpy as np
import time
from scipy.constants import c, e, m_p
from scipy.constants import e as qe
from numpy.fft import fft, ifft

from PyHEADTAIL.particles.slicing import UniformBinSlicer as PyHTUniformBinSlicer
from PyHEADTAIL.particles.particles import Particles as PyHTParticles
from PyHEADTAIL.impedances.wakes import CircularResonator as PyHTCircularResonator

from PyHEADTAIL.impedances.wakes import WakeField as PyHTWakeField
from PyHEADTAIL.machines.synchrotron import Synchrotron as PyHTSynchrotron

import xfields as xf
import xtrack as xt

# Machine settings

n_turns_wake = 3
flatten = False

n_macroparticles = 100000 # per bunch
intensity = 2.3e11

alpha = 53.86**-2

p0 = 7000e9 * e / c

accQ_x = 0.31
accQ_y = 0.32
Q_s = 2.1e-3
chroma = 0

h_RF = 600
bunch_spacing_buckets = 5
n_bunches = 12
n_slices = 250
plot_on = True
beta_x = 100
beta_y = 100

h_bunch = h_RF
circumference = 26658.883 / 35640 * h_RF

epsn_x = 2e-6
epsn_y = 2e-6
sigma_z = 0.09

machine = PyHTSynchrotron(
        optics_mode='smooth', circumference=circumference,
        n_segments=1, s=None, name=None,
        alpha_x=None, beta_x=beta_x, D_x=0,
        alpha_y=None, beta_y=beta_y, D_y=0,
        accQ_x=accQ_x, accQ_y=accQ_y, Qp_x=chroma, Qp_y=chroma,
        app_x=0, app_y=0, app_xy=0,
        alpha_mom_compaction=alpha, longitudinal_mode='linear',
        h_RF=np.atleast_1d(h_RF), p0=p0,
        charge=e, mass=m_p, wrap_z=False, Q_s=Q_s)
transverse_map = machine.transverse_map.segment_maps[0]





# Filling scheme

filling_scheme = [i*bunch_spacing_buckets for i in range(n_bunches)]

# Initialise beam
allbunches = machine.generate_6D_Gaussian_bunch(n_macroparticles, intensity,
                                                epsn_x, epsn_y, sigma_z=sigma_z,
                                                filling_scheme=filling_scheme,
                                                matched=False)

bucket_id_set = list(set(allbunches.bucket_id))
bucket_id_set.sort()

bucket_length = machine.circumference / h_RF
z_all = -allbunches.bucket_id * bucket_length + allbunches.z

amplitude = 1e-3
wavelength = 2
allbunches.x = amplitude * np.sin(2 * np.pi * z_all / wavelength)
allbunches.xp *= 0

# allbunches.x[z_all < 0] = 0

for b_id in bucket_id_set:
    mask = allbunches.bucket_id == b_id
    z_centroid = np.mean(allbunches.z[mask])
    z_std = np.std(allbunches.z[mask])
    mask_tails = mask & (np.abs(allbunches.z - z_centroid) > z_std)
    allbunches.x[mask_tails] = 0
    # if b_id != 0:
    #     allbunches.x[mask] = 0

beam = PyHTParticles(macroparticlenumber=allbunches.macroparticlenumber,
                 particlenumber_per_mp=allbunches.particlenumber_per_mp,
                 charge=allbunches.charge, mass=allbunches.mass,
                 circumference=allbunches.circumference, gamma=allbunches.gamma,
                 coords_n_momenta_dict=dict(x=allbunches.x.copy(),
                                            y=allbunches.y.copy(),
                                            xp=allbunches.xp.copy(),
                                            yp=allbunches.yp.copy(),
                                            z=z_all.copy(),
                                            dp=allbunches.dp.copy(),
                 ))

# Initialise wakes
slicer = PyHTUniformBinSlicer(n_slices, z_cuts=(-0.5*bucket_length, 0.5*bucket_length),
                          circumference=machine.circumference, h_bunch=h_bunch)

# pipe radius [m]
b = 13.2e-3
# length of the pipe [m]
L = 100000.
# conductivity of the pipe 1/[Ohm m]
sigma = 1. / 7.88e-10

# wakes = CircularResistiveWall(b, L, sigma, b/c, beta_beam=machine.beta)
wakes = PyHTCircularResonator(R_shunt=135e6, frequency=1.97e9*0.6, Q=31000/100, n_turns_wake=n_turns_wake)

# mpi_settings = 'circular_mpi_full_ring_fft'
# wake_field = WakeField(slicer, wakes, mpi=mpi_settings, Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)

R_s = wakes.R_shunt
Q = wakes.Q
f_r = wakes.frequency
omega_r = 2 * np.pi * f_r
alpha_t = omega_r / (2 * Q)
omega_bar = np.sqrt(omega_r**2 - alpha_t**2)
p0_SI = machine.p0


mpi_settings = False
mpi_settings = 'memory_optimized'
mpi_settings = 'linear_mpi_full_ring_fft'
wake_field = PyHTWakeField(slicer, wakes, mpi=mpi_settings)

# Wake full beam

n_buckets_slicer = max(filling_scheme) + 2
n_buckets_slicer = max(filling_scheme) + 1

slicer_full_beam = PyHTUniformBinSlicer(n_buckets_slicer * slicer.n_slices,
                                    z_cuts=((0.5 - n_buckets_slicer)*bucket_length, 0.5*bucket_length),
                                    circumference=machine.circumference, h_bunch=h_bunch)
slicer_full_beam.force_absolute = True

wakes_full_beam = PyHTCircularResonator(R_shunt=wakes.R_shunt, frequency=wakes.frequency, Q=wakes.Q, n_turns_wake=wakes.n_turns_wake)
wake_field_full_beam = PyHTWakeField(slicer_full_beam, wakes_full_beam, mpi=False)


plt.close('all')

fig0, (ax00, ax01) = plt.subplots(2, 1)
ax01.sharex(ax00)

x_at_wake_allbunches = []
xp_before_wake_allbunches = []
xp_after_wake_allbunches = []
slice_set_before_wake_allbunches = []
slice_set_after_wake_allbunches = []

x_at_wake_beam = []
xp_before_wake_beam = []
xp_after_wake_beam = []
slice_set_before_wake_beam = []
slice_set_after_wake_beam = []

n_turns = n_turns_wake
store_charge_per_mp = allbunches.charge_per_mp
store_particlenumber_per_mp = allbunches.particlenumber_per_mp

z_source_matrix_multiturn = np.zeros((n_bunches, n_slices, n_turns))
dipole_moment_matrix_multiturn = np.zeros((n_bunches, n_slices, n_turns))

from wakefield import Wakefield, TempResonatorFunction



particles = xt.Particles(
    mass0=xt.PROTON_MASS_EV,
    gamma0=beam.gamma,
    x=beam.x,
    px=beam.xp,
    y=beam.y,
    py=beam.yp,
    zeta=beam.z,
    delta=beam.dp,
    weight = beam.particlenumber_per_mp,
)


color_list = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
for i_turn in range(n_turns):

    needs_restore = False
    x_beam_before_wake = beam.x.copy()
    x_allbunches_before_wake = allbunches.x.copy()
    # if i_turn != 0:
    #     beam.x[:] = 0
    #     allbunches.x[:] = 0
    #     needs_restore = True

    allbunches.clean_slices()
    slice_set_before_wake_allbunches.append(allbunches.get_slices(slicer,
                        statistics=['mean_x', 'mean_xp', 'mean_y', 'mean_yp']))
    allbunches.clean_slices()

    # Continue with PyHEADTAIL
    beam.clean_slices()
    slice_set_before_wake_beam.append(beam.get_slices(slicer_full_beam,
                        statistics=['mean_x', 'mean_xp', 'mean_y', 'mean_yp']))
    beam.clean_slices()

    x_at_wake_allbunches.append(allbunches.x.copy())
    xp_before_wake_allbunches.append(allbunches.xp.copy())

    x_at_wake_beam.append(beam.x.copy())
    xp_before_wake_beam.append(beam.xp.copy())

    # Back to pyheadtail

    wake_field.track(allbunches)
    wake_field_full_beam.track(beam)

    allbunches.clean_slices()
    slice_set_after_wake_allbunches.append(allbunches.get_slices(slicer,
                        statistics=['mean_x', 'mean_xp', 'mean_y', 'mean_yp']))
    allbunches.clean_slices()

    beam.clean_slices()
    slice_set_after_wake_beam.append(beam.get_slices(slicer_full_beam,
                        statistics=['mean_x', 'mean_xp', 'mean_y', 'mean_yp']))
    beam.clean_slices()

    xp_after_wake_allbunches.append(allbunches.xp.copy())
    xp_after_wake_beam.append(beam.xp.copy())

    if needs_restore:
        beam.x[:] = x_beam_before_wake
        allbunches.x[:] = x_allbunches_before_wake

    transverse_map.track(allbunches)
    transverse_map.track(beam)

    z_centers = slice_set_before_wake_beam[-1].z_centers
    mean_x_slice = slice_set_before_wake_beam[-1].mean_x
    num_charges_slice = slice_set_before_wake_beam[-1].charge_per_slice/e
    dip_moment_slice = num_charges_slice * mean_x_slice

    for i_bunch in range(n_bunches):
        i_bucket = bucket_id_set[::-1][i_bunch]
        print(f"{i_bunch=} {i_bucket=}")
        i_z_bunch_center = np.argmin(np.abs(z_centers - (-i_bucket * bucket_length)))

        i_z_a_bunch = i_z_bunch_center - n_slices//2
        if i_z_a_bunch == -1: # Rounding error
            i_z_a_bunch = 0

        dipole_moment_matrix_multiturn[i_bunch, :, i_turn] = dip_moment_slice[
                                            i_z_a_bunch : i_z_a_bunch+n_slices]

        z_source_matrix_multiturn[i_bunch, :, i_turn] = z_centers[
                                        i_z_a_bunch : i_z_a_bunch+n_slices]

    if plot_on:
        ax00.plot(slice_set_after_wake_beam[-1].z_centers,
                slice_set_after_wake_beam[-1].mean_x, '-', color=color_list[i_turn],
                label=f"turn {i_turn}")

        ax01.plot(
            slice_set_after_wake_beam[-1].z_centers,
            slice_set_after_wake_beam[-1].mean_xp - slice_set_before_wake_beam[-1].mean_xp,
            '-', color=color_list[i_turn], label=f"turn {i_turn}")

if plot_on:
    ax00.legend()

z_source_matrix_multiturn = z_source_matrix_multiturn[:, :, ::-1] # last turn on top
dipole_moment_matrix_multiturn = dipole_moment_matrix_multiturn[:, :, ::-1] # last turn on top

print(f'Circumference occupancy {n_bunches * bunch_spacing_buckets/h_RF*100:.2f} %')


n_bunches_wake = 120 # Can be longer than filling scheme

wf = Wakefield(
    source_moments=['num_particles', 'x', 'y'],
    kick='px',
    scale_kick=None, # The kick is scaled by position of the particle for quadrupolar, would be None for dipolar
    function=TempResonatorFunction(R_shunt=wakes.R_shunt, frequency=wakes.frequency, Q=wakes.Q),
    zeta_range=(-0.5*bucket_length, 0.5*bucket_length), # These are [a, b] in the paper
    num_slices=n_slices, # Per bunch, this is N_1 in the paper
    bunch_spacing_zeta=bunch_spacing_buckets*bucket_length, # This is P in the paper
    num_bunches=n_bunches_wake, # This is N_S
    num_turns=n_turns_wake,
    circumference=circumference,
    log_moments=['px'],
    _flatten=flatten
)

slicer_after_wf = xf.UniformBinSlicer(
    zeta_range=(-0.5*bucket_length, 0.5*bucket_length),
    num_slices=n_slices,
    i_bunch_0=0, num_bunches=n_bunches_wake,
    bunch_spacing_zeta=bunch_spacing_buckets*bucket_length,
    moments=['x', 'px'])

betatron_map = xt.LineSegmentMap(
    length=circumference, betx=beta_x, bety=beta_y,
    qx=accQ_x, qy=accQ_y,
    longitudinal_mode='frozen')

line = xt.Line(elements=[wf, slicer_after_wf, betatron_map],
               element_names=['wf', 'slicer_after_wf', 'betatron_map'])
particle_ref = xt.Particles(p0c=particles.p0c, mass0=particles.mass0, q0=particles.q0)
line.build_tracker()

line.track(particles, num_turns=n_turns)

ax00.plot(wf.slicer.zeta_centers, wf.slicer.mean('x'), '.', color='b')
ax01.plot(wf.slicer.zeta_centers, slicer_after_wf.mean('px') - wf.slicer.mean('px'), '.', color='b')

if plot_on:
    plt.figure(200)
    plt.plot(wf.z_wake.T, wf.G_aux.T * (-e**2 / (p0_SI * c)), alpha=0.5)

t0 = time.time()
line.track(particles)
t1 = time.time()
print(f"Xsuite turn time {(t1-t0)*1e3:.2f} ms")

t0 = time.time()
wake_field.track(allbunches)
transverse_map.track(allbunches)
allbunches.get_slices(slicer,
                      statistics=['mean_x', 'mean_xp', 'mean_y', 'mean_yp'])
t1 = time.time()
print(f"PyHEADTAIL turn time {(t1-t0)*1e3} ms")

plt.show()

