import pandas as pd
import numpy as np
import plotly.graph_objects as go


df = pd.read_csv('/Users/bryanwgranger/Documents/bioCM/scBCRseq/workflow/test/seurat_obj_metadata.csv', index_col=[0])

umap = pd.read_csv('/Users/bryanwgranger/Documents/bioCM/scBCRseq/workflow/test/seurat_umap.csv', index_col=[0])


df = df.merge(umap, left_on=df.index, right_on=umap.index)
print(df.columns)
print(list(range(len(df['Cell'].unique()))))
print(df['Cell'].unique())
df['Frequency'] = df['Frequency'].fillna(0)

cell_dict = {k:v for k,v in zip(df['Cell'].unique(), list(range(9)))}
print(cell_dict)



fig = go.Figure()


fig.add_trace(go.Scatter(x=df['UMAP_1'], y=df['UMAP_2'],
                         mode='markers',
                         hovertext=df['Cell'],
                         marker_color=[cell_dict[u] for u in df['Cell']],
                         marker_colorscale='viridis',
                         # marker_size= df['Frequency']

                         ))

# fig.add_trace(go.Scatter(x=df[df['sample'] == 'tumor']['UMAP_1'], y=df[df['sample'] == 'tumor']['UMAP_2'],
#                          mode='markers',
#                          hovertext=df['Cell'],
#                          marker_color='blue',
#                          name='tumor'
#
#                          ))


fig.show()

