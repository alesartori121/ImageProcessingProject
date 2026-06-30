function final_tumor_mask = isolate_tumor_mass(candidate_mask)
    extracted_mass = bwareafilt(candidate_mask, 1);
    final_tumor_mask = imfill(extracted_mass, 'holes');
end