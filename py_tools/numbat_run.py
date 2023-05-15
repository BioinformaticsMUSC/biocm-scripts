import os
import chevron

with open('set3.txt', 'r') as s:
    sample_names = [h.strip() for h in s.readlines() if h != 'script\n']

for samp in sample_names:

    data = {
        "job_id": f"{samp}_numbat_prep",
        "label": samp,
        "samples": samp,
        "bam_path": f"/mnt/sara_gbm/Set3/run/{samp}/outs/possorted_genome_bam.bam",
        "barcode_path": f"/mnt/sara_gbm/Set3/run/{samp}/outs/barcodes.tsv",
        "out_path": f"/mnt/Sara/Set3/numbat/{samp}"
    }

    with open('numbat_prep_template.txt', 'r') as f:
        a = chevron.render(f, data=data)

    with open(f'numbat/{samp}_numbat_prep.pbs', 'w') as outfile:
        outfile.write(a)
