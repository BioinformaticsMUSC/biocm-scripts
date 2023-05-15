#!/usr/bin/env python
def main(args):

    if not args.count and not args.multi:
        args.count = True

    if args.count and args.multi:
        raise RuntimeError("Please select ONE of either count (--count) or multi (--multi)")

    if args.count:

        #make sure output file ends with .pbs
        if not args.output.endswith(".pbs"):
            args.output = args.output + ".pbs"

        #fastqs required
        if not args.fastqs:
            raise RuntimeError("Please provide directory to fastq files (-f or --fastqs)")

        #ref path required
        if not args.transcriptome:
            raise RuntimeError("Please provide path to 10X reference (-t or --transcriptome)")

        # convert filepaths if needed
        res_dir = args.results_directory if args.results_directory else "/mnt"
        res_dir = res_dir.replace("/zfs/musc3", "/mnt")
        fastq_dir = args.fastqs.replace("/zfs/musc3", "/mnt")
        ref_dir = args.transcriptome.replace("/zfs/musc3", "/mnt")
        with open(args.output, "w") as outfile:
            #PBS header
            outfile.write(f"#PBS -N {args.sample_name}_count\n")
            outfile.write("#PBS -l select=1:ncpus=24:mem=200gb,walltime=48:00:00\n")
            outfile.write("#PBS -m abe\n\n")




            outfile.write(f"singularity exec -B /zfs/musc3:/mnt --pwd {res_dir} /zfs/musc3/singularity_images/biocm-cellranger_latest.sif \\\n")
            outfile.write("\tcellranger count \\\n")
            outfile.write(f"\t--id={args.sample_name}_count \\\n")
            outfile.write(f"\t--transcriptome={ref_dir} \\\n")
            outfile.write(f"\t--fastqs={fastq_dir} \\\n")
            outfile.write("\t--localmem=150")
        print(f"PBS script written and saved at {args.output}")

    if args.multi:
        # convert filepaths if needed
        res_dir = args.results_directory if args.results_directory else "/mnt"
        res_dir = res_dir.replace("/zfs/musc3", "/mnt")
        csv_dir = args.csv.replace("zfs/musc3", ".mnt")

        if not csv_dir.endswith(".csv"):
            csv_dir = csv_dir + ".csv"

        with open(args.output, "w") as outfile:
            #PBS header
            outfile.write(f"#PBS -N {args.sample_name}_multi\n")
            outfile.write("#PBS -l select=1:ncpus=24:mem=200gb,walltime=48:00:00\n")
            outfile.write("#PBS -m abe\n\n")

            outfile.write(f"singularity exec -B /zfs/musc3:/mnt --pwd {res_dir} /zfs/musc3/singularity_images/biocm-cellranger_latest.sif \\\n")
            outfile.write("\tcellranger multi \\\n")
            outfile.write(f"\t--id={args.sample_name}_count \\\n")
            outfile.write(f"\t--csv={args.csv} \\\n")

        print(f"PBS script written and saved at {args.output}")

    pass
if __name__ == "__main__":

    # setting the hyper parameters
    import argparse

    parser = argparse.ArgumentParser(description='Create Anndata file from prepared Seurat directory',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--count", action="store_true")
    parser.add_argument("--multi", action="store_true")
    parser.add_argument('-r', '--results_directory', default=None, required=False)
    parser.add_argument('-f', '--fastqs', default=None, required=False)
    parser.add_argument('-t', '--transcriptome', default=None, required=False)
    parser.add_argument('-s', '--sample_name', default=None, required=True)
    parser.add_argument('-c', '--csv', default=None, required=False)
    parser.add_argument('-o', '--output', default='cellranger_job.pbs')

    args = parser.parse_args()

    main(args)