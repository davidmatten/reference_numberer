



occupancy = {}
all_models_occupation = {}
for fldr in dirs:
    all_models_occupation[fldr] = {}
    print("now we will look at {}".format(fldr))
    if fldr == "Mature_gp120_seq1-99":
        continue
    if fldr == "HcBc2_HBX2_11":
        continue
    all_glycosylated = []
    all_not_glycosylated = []
    for i in range(1, 11):
        all_models_occupation[fldr][i] = {}
        print("repeat number: {}".format(i))
        full_path = os.path.join(root, fldr, str(i))

        positions_fn = os.path.join(full_path, "parsed_log.log")
        print(positions_fn)

        glycosylated_positions = []
        not_glyco_pos = []

        record_not = False
        record_glyco = True

        with open(positions_fn, "r") as fh:
            for line in fh:
                if "Number not glycosylated" in line:
                    record_not = True
                    record_glyco = False
                if record_glyco:
                    glycosylated_positions.append(line.strip())
                if record_not:
                    not_glyco_pos.append(line.strip())
        glycosylated_positions.pop()
        glycosylated_positions.remove(glycosylated_positions[0])
        glycosylated_positions = [int(i) for i in glycosylated_positions]
        all_models_occupation[fldr][i]['occupied'] = glycosylated_positions
        not_glyco_pos.remove(not_glyco_pos[0])
        not_glyco_pos = [int(i) for i in not_glyco_pos]
        all_models_occupation[fldr][i]['not_occupied'] = not_glyco_pos

        print(glycosylated_positions)
        print(not_glyco_pos)
        all_glycosylated += glycosylated_positions
        all_not_glycosylated += not_glyco_pos

    all_sites = all_glycosylated + all_not_glycosylated
    list_set_all_sites = list(set(all_sites))

    tmp_pos_prcnt_dct = {}
    for pos in list_set_all_sites:
        prcnt = (all_glycosylated.count(pos) / 10.0) * 100.0
        tmp_pos_prcnt_dct[pos] = prcnt
    occupancy[fldr] = {}
    occupancy[fldr]['modelled_percent'] = tmp_pos_prcnt_dct.copy()

    hxb2_seq = dct['HXB2']
    print(hxb2_seq)

    print(dct)
    print(fldr)
    sample_seq = dct[fldr]
    print(sample_seq)
    sample_seq_len = len(sample_seq)
    monomer_len = sample_seq_len / 3.0

    seq_pos, hxb2_pos = 0, 0
    in_hxb2_gap = False
    in_seq_insert_region = False
    in_hxb2_gap_counter, in_seq_gap_counter = 0, 0
    done_positions = []
    sample_hxb2_positions = []
    translate = {}
    for i in range(len(hxb2_seq)):
        if sample_seq[i] != "-":
            seq_pos += 1
            in_seq_insert_region = False
            in_seq_gap_counter = 0
        if sample_seq[i] == '-':
            in_seq_insert_region = True
            in_seq_gap_counter += 1
        if hxb2_seq[i] != '-':
            hxb2_pos += 1
            in_hxb2_gap = False
            in_hxb2_gap_counter = 0
        if hxb2_seq[i] == '-':
            in_hxb2_gap = True
            in_hxb2_gap_counter += 1

        if ((seq_pos in glycosylated_positions) or (seq_pos in not_glyco_pos)) and seq_pos not in done_positions:
            # get percent glycosylated here for that position
            done_positions.append(seq_pos)
            if in_hxb2_gap_counter != 0:
                ####################################################################################################sample_hxb2_positions.append(float(str(hxb2_pos) + "." + str(in_hxb2_gap_counter)))
                sample_hxb2_positions.append(str(str(hxb2_pos) + "." + str(in_hxb2_gap_counter)))
                translate[str(str(hxb2_pos) + "." + str(in_hxb2_gap_counter))] = seq_pos
            else:
                sample_hxb2_positions.append(str(str(hxb2_pos)))
                translate[str(str(hxb2_pos))] = seq_pos
    all_models_occupation[fldr]['translation'] = translate
    occupancy[fldr]['translate'] = translate.copy()
    occupancy[fldr]['positions'] = sample_hxb2_positions
