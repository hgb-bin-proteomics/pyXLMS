from pyteomics import mzid

def parse_scan_nr(spectrum_id):
    return int(spectrum_id.split("scan=")[1])

def read_mzid(file):
    csms = list()
    pyteomics_mzid = mzid.MzIdentML(file)
    for item in pyteomics_mzid:
        csm_id = None
        scan = None
        filename = None
        peptide_a = None
        pos_a = None
        peptide_b = None
        pos_b = None
        if "spectrumID" in item:
            scan = parse_scan_nr(item["spectrumID"])
        if "location" in item:
            filename = str(item["location"]).strip()
        if "SpectrumIdentificationItem" in item:
            for subitem in item["SpectrumIdentificationItem"]:
                if "rank" in subitem:
                    if int(subitem["rank"]) > 1:
                        continue
                if "cross-link spectrum identification item" in subitem:
                    if csm_id is None:
                        csm_id = int(float(subitem["cross-link spectrum identification item"]))
                        if "PeptideSequence" in subitem:
                            peptide_a = format_sequence(subitem["PeptideSequence"])
                        if "Modification" in subitem:
                            for mod in subitem["Modification"]:
                                if "name" in mod:
                                    if str(mod["name"]).strip().upper() in CROSSLINKERS:
                                        if "location" in mod:
                                            pos_a = int(mod["location"])
                    elif csm_id == int(float(subitem["cross-link spectrum identification item"])):
                        if "PeptideSequence" in subitem:
                            peptide_b = format_sequence(subitem["PeptideSequence"])
                        if "Modification" in subitem:
                            for mod in subitem["Modification"]:
                                if "name" in mod:
                                    if str(mod["name"]).strip().upper() in CROSSLINKERS:
                                        if "location" in mod:
                                            pos_b = int(mod["location"])
        if None not in [csm_id, scan, filename, peptide_a, pos_a, peptide_b, pos_b]:
            csms.append(create_csm())
    pyteomics_mzid.close()
    return csms
