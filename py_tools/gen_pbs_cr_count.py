import os
num_list = [78, 79,80,81]
pi = "tsao-hto"
transcriptome = "mouse"


save_dir = f"/Users/bryanwgranger/biocm/projects/{pi}/scripts/cluster"

fastq_path = f"/scratch1/bryangranger/{pi}/fastqs"
cluster_run_dir = f"/scratch1/bryangranger/{pi}/run"
expect_cells = "10000"
save_dir = f"/Users/bryanwgranger/biocm/projects/{pi}/scripts/cluster"

if not os.path.exists(save_dir):
    os.mkdir(save_dir)

for n in num_list:
    num = str(n)
    id = f"{pi}_{num}"
    filename = id + ".pbs"

    with open(os.path.join(save_dir, filename), 'w') as outfile:
        outfile.write(f"#PBS -N {id}\n")
        outfile.write("#PBS -l select=1:ncpus=16:mem=60gb:interconnect=1g,walltime=48:00:00\n")
        outfile.write("#PBS -m abe\n#PBS -M grangerb@musc.edu\n\n")
        outfile.write("source ~/.bashrc\n\n")
        outfile.write(f"cd {cluster_run_dir}\n\n")
        outfile.write(f"cellranger count --id={id} \\\n")
        if transcriptome == "human":
            outfile.write("--transcriptome=/scratch1/bryangranger/transcriptome/refdata-gex-GRCh38-2020-A \\\n")
        elif transcriptome == "mouse":
            outfile.write("--transcriptome=/scratch1/bryangranger/transcriptome/refdata-gex-mm10-2020-A \\\n")
        outfile.write(f"--fastqs={fastq_path} \\\n")
        outfile.write(f"--sample=7166-MR-{num} \\\n")
        outfile.write(f"--expect-cells={expect_cells}")