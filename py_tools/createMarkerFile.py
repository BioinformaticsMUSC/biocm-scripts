#Generate Garnet marker file from FindAllMarkers output from Seurat
import pandas as pd

def createMarkerFile(args):
    df = pd.read_csv(args.input, index_col=[0])
    cell_types = set()
    with open(args.output, 'w') as outfile:
        for cluster in df['cluster'].unique():
            t_df = df[df['cluster'] == cluster]
            outfile.write(">" + cluster + "\nexpressed: ")
            t_df_over = t_df[t_df['avg_log2FC'] > 0].sort_values('avg_log2FC')[:args.n_markers]
            t_df_under = t_df[t_df['avg_log2FC'] < 0].sort_values('avg_log2FC')[:args.n_markers]
            if len(t_df_over) > 0:
                cell_types.add(cluster)
                for i, idx in enumerate(t_df_over.index):
                    outfile.write(t_df_over.loc[idx, 'gene'])
                    if len(t_df_over) > i+1:
                        outfile.write(', ')
                    else:
                        outfile.write('\n')
            if len(t_df_under) > 0:
                cell_types.add(cluster)
                outfile.write('not expressed: ')
                for i, idx in enumerate(t_df_under.index):
                    outfile.write(t_df_under.loc[idx, 'gene'])
                    if len(t_df_under) > i+1:
                        outfile.write(', ')
                    else:
                        outfile.write('\n')
            outfile.write('\n')
    print(f"Generated Garnett marker file: {args.output}")
    print("Included markers for:")
    for ct in cell_types:
        print("\t-" + ct)


if __name__ == "__main__":

    # setting the hyper parameters
    import argparse

    parser = argparse.ArgumentParser(description='Create Garnett Marker Files from FindAllMarkers',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', default=None, required=True)
    parser.add_argument('--output', default='garnet_marker_file.txt')
    parser.add_argument('--n_markers', default=10, type=int)

    args = parser.parse_args()

    createMarkerFile(args)