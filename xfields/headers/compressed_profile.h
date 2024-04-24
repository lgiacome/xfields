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
    /*gpuglmem*/ double* zeta_centers,
    /*gpuglmem*/ int64_t* i_bunch_particles,
    /*gpuglmem*/ int64_t* i_edge_particles,
    /*gpuglmem*/ double* out
    ){

    int64_t const _N_S = CompressedProfileData_get__N_S(el);
    int64_t const _N_aux = CompressedProfileData_get__N_aux(el);
    int64_t const num_turns = CompressedProfileData_get_num_turns(el);
    int64_t const num_slices = CompressedProfileData_get__N_1(el);
    const double dzeta = zeta_centers[1] - zeta_centers[0];

    //start_per_particle_block (part0->part)

        const int64_t ipart = part->ipart;
        const double zeta = LocalParticle_get_zeta(part);
        const int64_t i_bunch = i_bunch_particles[ipart];
        const int64_t i_edge = i_edge_particles[ipart];

        const int64_t i_start_in_moments_data = (_N_S - i_bunch - 1) * _N_aux;
        double val_left = 0;
        double val_right = 0;
        double zeta_left = 0;
        double zeta_right = 0;

        double rr = 0;

        for(int i_turn=0; i_turn<num_turns; i_turn++){
            int64_t turn_offset = data_shape_2 * (i_turn + data_shape_1*(data_shape_0-1));
            if(i_edge == 0){
                val_right =  data[i_start_in_moments_data + i_edge + turn_offset];
                val_left =  -val_right;
                zeta_left = zeta_centers[0] - dzeta;
                zeta_right = zeta_centers[i_edge];
            }else if(i_edge == num_slices){
                val_left =  data[i_start_in_moments_data + (i_edge-1) + turn_offset];
                val_right = -val_left;
                zeta_left = zeta_centers[i_edge-1];
                zeta_right = zeta_centers[num_slices-1] + dzeta;
            }else{
                val_left =  data[i_start_in_moments_data + (i_edge-1) + turn_offset];
                val_right =  data[i_start_in_moments_data + i_edge + turn_offset];
                zeta_left = zeta_centers[i_edge-1];
                zeta_right = zeta_centers[i_edge];
            }

            rr = rr + val_left*(zeta_right - zeta)/dzeta + val_right*(zeta - zeta_left)/dzeta;
        }
        out[ipart] = rr;

    //end_per_particle_block

    }




#endif // XFIELDS_COPRESSED_PROFILE_H
