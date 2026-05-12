function [wm_gm_chunks, chunk_names, chunk_areas] = gm_wm_chunking(subject_ID, reference, cur_letter)

% OUTPUTS
% - wm_gm_chunks stores index (in cur_elec_contact_ind) of first contact in each chunk
% - chunk_names stores names of first contact in each chunk
% NOTE: 
% - last element of gm_wm_chunks is always numel(cur_elec_contact_ind) + 1
% - last element of chunk_names is always "END"

    [~, data_base_dir, ~, ~, ~, ~, ~] = tw_setup(subject_ID, reference);

    if strcmp(reference,'Ground')
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup.mat',data_base_dir,subject_ID,subject_ID,1), "elec_region");
    else
        load(sprintf('%s/%s/%s_MAD_SES%d_Setup_%s.mat',data_base_dir,subject_ID,subject_ID,1,reference), "channel_region");
        elec_region = channel_region;
    end

    [cur_elec_contact_ind, cur_elec_contact_names] = get_single_probe_contacts(reference, subject_ID, cur_letter);

    wm_gm_labels = strings(numel(cur_elec_contact_ind),1);
    for cnum = 1:numel(cur_elec_contact_ind)
        cur_region = elec_region{cur_elec_contact_ind(cnum)};
        if contains(cur_region, 'WM', 'IgnoreCase', true)
            wm_gm_labels(cnum) = "WM";
        elseif contains(cur_region, 'undefined', 'IgnoreCase', true) || contains(cur_region, 'out', 'IgnoreCase', true)
            wm_gm_labels(cnum) = "other";
        elseif contains(cur_region, 'bolt', 'IgnoreCase', true) || contains(cur_region, 'bone', 'IgnoreCase', true)
            wm_gm_labels(cnum) = "other";
        else
            wm_gm_labels(cnum) = "GM";
        end

    end

    wm_gm_chunks = [1];
    chunk_names = [cur_elec_contact_names(1)];
    chunk_areas = [wm_gm_labels(1)];

    cnum = 2;
    while cnum <= numel(cur_elec_contact_names)
        %contact = cur_elec_contact_ind(cnum);
        %prev_contact = cur_elec_contact_ind(cnum-1);
        cur_name = cur_elec_contact_names(cnum);
        cur_area = wm_gm_labels(cnum);
        prev_area = wm_gm_labels(cnum-1);

        if ~strcmpi(cur_area,prev_area)
            wm_gm_chunks = [wm_gm_chunks; cnum];
            chunk_names = [chunk_names; cur_name];
            chunk_areas = [chunk_areas; wm_gm_labels(cnum)];
        end

        cnum = cnum + 1;
    end

    wm_gm_chunks(end + 1) = numel(cur_elec_contact_ind) + 1;
    chunk_names(end + 1) = "END";
    chunk_areas(end + 1) = wm_gm_labels(end);

end