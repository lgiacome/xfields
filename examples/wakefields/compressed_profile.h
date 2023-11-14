#ifndef XFIELDS_COPRESSED_PROFILE_H
#define XFIELDS_COPRESSED_PROFILE_H

void CompressedProfile_track_local_particle(CompressedProfileData el,
                                           LocalParticle *part0){

                                           }

void CompressedProfile_interp_result(
    CompressedProfileData el, LocalParticle *part0,
    int64_t data_shape_0,
    int64_t data_shape_1,
    int64_t data_shape_2,
    /*gpuglmem*/ double* data,
    /*gpuglmem*/ int64_t* i_bunch_particles,
    /*gpuglmem*/ int64_t* i_slice_particles,
    /*gpuglmem*/ double* out
    ){

    int64_t const _N_S = CompressedProfileData_get__N_S(el);
    int64_t const _N_aux = CompressedProfileData_get__N_aux(el);

    //start_per_particle_block (part0->part)

        const int64_t ipart = part->ipart;
        const int64_t i_bunch = i_bunch_particles[ipart];
        const int64_t i_slice = i_slice_particles[ipart];

        const int64_t i_start_in_moments_data = (_N_S - i_bunch - 1) * _N_aux;

        const double rr = data[
            data_shape_2 * data_shape_1 * (data_shape_0 - 1) +
            // data.shape[0] * (i_turn) # i_turn is zero for the result
            i_start_in_moments_data + i_slice];

        out[ipart] = rr;

    //end_per_particle_block

    }




#endif // XFIELDS_COPRESSED_PROFILE_H