
file = '/Users/bryanwgranger/biocm/biocm-tools/r_tools/interactive_scripts/01_load_seurat_data_interactive.R'

with open(file, "r") as f:
    r = f.readlines()

glue_file = '/Users/bryanwgranger/Desktop/glue_file.txt'

with open(glue_file, 'w') as g:
    for line in r:
        line = line.replace('\n','\\n')
        line = line.replace("'", '"')
        new_line = "'" + line + "',\n"
        print(new_line)
        g.write(new_line)