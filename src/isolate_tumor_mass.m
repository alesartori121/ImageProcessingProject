function final_tumor_mask = isolate_tumor_mass(candidate_mask)
% ISOLATE_TUMOR_MASS Keeps largest connected component and fills holes.
    extracted_mass = bwareafilt(candidate_mask, 1);
    final_tumor_mask = imfill(extracted_mass, 'holes');
end