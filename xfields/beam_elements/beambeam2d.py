# copyright ################################# #
# This file is part of the Xfields Package.   #
# Copyright (c) CERN, 2021.                   #
# ########################################### #

from sqlite3 import paramstyle
import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root

class BeamBeamBiGaussian2D(xt.BeamElement):

    _xofields = {

        'ref_shift_x': xo.Float64,
        'ref_shift_y': xo.Float64,

        'other_beam_shift_x': xo.Float64,
        'other_beam_shift_y': xo.Float64,

        'post_subtract_px': xo.Float64,
        'post_subtract_py': xo.Float64,

        # TODO this could become other_beam_q0, other_beam_beta0 (to be done also in 6D)
        'q0_other_beam': xo.Float64,
        'beta0_other_beam': xo.Float64,

        'other_beam_num_particles': xo.Float64,

        'other_beam_Sigma_11': xo.Float64,
        'other_beam_Sigma_13': xo.Float64,
        'other_beam_Sigma_33': xo.Float64,

        'min_sigma_diff': xo.Float64,

    }

    extra_sources= [
        _pkg_root.joinpath('headers/constants.h'),
        _pkg_root.joinpath('headers/sincos.h'),
        _pkg_root.joinpath('headers/power_n.h'),
        _pkg_root.joinpath('fieldmaps/bigaussian_src/complex_error_function.h'),
        '#define NOFIELDMAP', #TODO Remove this workaround
        _pkg_root.joinpath('fieldmaps/bigaussian_src/bigaussian.h'),
        _pkg_root.joinpath('beam_elements/beambeam_src/beambeam2d.h'),
    ]


    def __init__(self,
                    q0_other_beam=None,

                    other_beam_num_particles=None,

                    other_beam_Sigma_11=None,
                    other_beam_Sigma_13=None,
                    other_beam_Sigma_33=None,

                    ref_shift_x=0,
                    ref_shift_y=0,

                    other_beam_shift_x=0,
                    other_beam_shift_y=0,

                    post_subtract_px=0,
                    post_subtract_py=0,

                    min_sigma_diff=1e-10,

                    config_for_update=None,

                    **kwargs):

        if '_xobject' in kwargs.keys():
            self.xoinitialize(**kwargs)
            return

        # Collective mode (pipeline update)
        if config_for_update is not None:
            raise NotImplementedError
            # To be implemented based on 6d implementation

        params = self._handle_init_old_interface(kwargs)

        self.xoinitialize(**kwargs)

        if self.iscollective:
            if not isinstance(self._buffer.context, xo.ContextCpu):
                raise NotImplementedError(
                    'BeamBeamBiGaussian3D only works with CPU context for now')

        # Mandatory sigmas
        assert other_beam_Sigma_11 is not None, ("`other_beam_Sigma_11` must be provided")
        assert other_beam_Sigma_33 is not None, ("`other_beam_Sigma_33` must be provided")

        # Coupling between transverse planes
        if other_beam_Sigma_13 is None:
            other_beam_Sigma_13 = 0

        if np.abs(other_beam_Sigma_13) > 0:
            raise NotImplementedError(
                "Coupled case not tested yet.")

        assert other_beam_num_particles is not None, ("`other_beam_num_particles` must be provided")
        self.other_beam_num_particles = other_beam_num_particles

        self.other_beam_Sigma_11 = other_beam_Sigma_11
        self.other_beam_Sigma_13 = other_beam_Sigma_13
        self.other_beam_Sigma_33 = other_beam_Sigma_33

        assert q0_other_beam is not None
        self.q0_other_beam = q0_other_beam

        self.ref_shift_x = ref_shift_x
        self.ref_shift_y = ref_shift_y

        self.other_beam_shift_x = other_beam_shift_x
        self.other_beam_shift_y = other_beam_shift_y

        self.post_subtract_px = post_subtract_px
        self.post_subtract_py = post_subtract_py

        self.min_sigma_diff = min_sigma_diff

    def _handle_init_old_interface(self, kwargs):

        params = {}

        if 'n_particles' in kwargs.keys():
            params['other_beam_num_particles'] = kwargs['n_particles']
            del kwargs['n_particles']

        if 'q0' in kwargs.keys():
            params['q0_other_beam'] = kwargs['q0']
            del kwargs['q0']

        if 'beta0' in kwargs.keys():
            params['beta0_other_beam'] = kwargs['beta0']
            del kwargs['beta0']

        if 'mean_x' in kwargs.keys():
            params['other_beam_shift_x'] = kwargs['mean_x']
            del kwargs['mean_x']

        if 'mean_y' in kwargs.keys():
            params['other_beam_shift_y'] = kwargs['mean_y']
            del kwargs['mean_y']

        if 'sigma_x' in kwargs.keys():
            params['other_beam_Sigma_11'] = kwargs['sigma_x']**2
            del kwargs['sigma_x']

        if 'sigma_y' in kwargs.keys():
            params['other_beam_Sigma_33'] = kwargs['sigma_y']**2
            del kwargs['sigma_y']

        if 'd_px' in kwargs.keys():
            params['post_subtract_px'] = kwargs['d_px']
            del kwargs['d_px']

        if 'd_py' in kwargs.keys():
            params['post_subtract_py'] = kwargs['d_py']
            del kwargs['d_py']

        return params


