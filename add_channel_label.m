function all_tracks_vec_out = add_channel_label(all_tracks_vec_in,channel)

all_tracks_vec_out = {};
for nwell = 1:length(all_tracks_vec_in)
    all_tracks_in = all_tracks_vec_in{nwell};
    phases = fieldnames(all_tracks_in);
    all_tracks_out = struct([]);
    for ph = 1:length(phases)
        all_tracks_in_phase = all_tracks_in.(phases{ph});
        all_tracks_out_phase = struct([]);
        for ntrack = 1:length(all_tracks_in_phase)
            track_in = all_tracks_in_phase(ntrack);
            track_fields = fieldnames(track_in);
                for nfield = 1:length(track_fields)
                    if sum(strcmp(track_fields{nfield},{'nf','nmi'})) == 1
                        all_tracks_out_phase(ntrack).(track_fields{nfield}).(channel) = track_in.(track_fields{nfield});
                    else
                        all_tracks_out_phase(ntrack).(track_fields{nfield}) = track_in.(track_fields{nfield});
                    end
                end
        end
        all_tracks_out(1).(phases{ph}) = all_tracks_out_phase;
    end
    all_tracks_vec_out{nwell} = all_tracks_out;
end

return