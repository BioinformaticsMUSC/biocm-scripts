import pandas as pd
import numpy as np
import plotly.graph_objects as go

df = pd.read_csv('/Users/bryanwgranger/Documents/bioCM/scBCRseq/vdj_metadata.csv', index_col=[0])
print(len(df))
umap = pd.read_csv('/Users/bryanwgranger/Documents/bioCM/scBCRseq/vdj_umap.csv', index_col=[0])

df = df.merge(umap, left_on=df.index, right_on=umap.index)


clones_p = df[df['sample'] == 'pbmc']['CTstrict'].value_counts()
clones_t = df[df['sample'] == 'tumor']['CTstrict'].value_counts()

clones = pd.merge(left=clones_p, right=clones_t, left_on=clones_p.index, right_on=clones_t.index, suffixes=('_pbmc', '_tumor')) \
    .sort_values(by='CTstrict_tumor', ascending=False)

fig = go.Figure()

fig.add_trace(go.Bar(x=clones['key_0'],
                     y=clones['CTstrict_tumor'],
                     name='Tumor'))
fig.add_trace(go.Bar(x=clones['key_0'],
                     y=clones['CTstrict_pbmc'],
                     name="PBMC"))
fig.update_layout(template='ggplot2')

fig.show()

clones_h = clones.sort_values(by='CTstrict_tumor')
fig2 = go.Figure()

fig2.add_trace(go.Bar(y=clones_h['key_0'],
                     x=clones_h['CTstrict_tumor'],
                     name='Tumor',
                      orientation='h'))
fig2.add_trace(go.Bar(y=clones_h['key_0'],
                     x=clones_h['CTstrict_pbmc'],
                     name="PBMC",
                      orientation='h'))
fig2.update_layout(template='ggplot2')
fig2.show()




df = df.merge(clones, left_on='CTstrict', right_on='key_0', how='outer')
df['CTstrict_pbmc'] = df['CTstrict_pbmc'].fillna(0)
df['CTstrict_tumor'] = df['CTstrict_tumor'].fillna(0)
print(df)
print(df['CTstrict_tumor'].value_counts())

df['color'] = 'lightblue'

for i in df.index:
    if df.loc[i, 'CTstrict_pbmc'] > 0:
        df.loc[i, 'color'] = 'red'

print(df['color'].value_counts())

#
fig3 = go.Figure(go.Scatter(x=df['UMAP_1'], y=df['UMAP_2'], mode='markers',
                            marker_color=df['color'],
                            # marker_size=1,
                            hovertext=df['Cell']))
fig3.update_layout(template='plotly_white')
fig3.show()


df_overlaps = df[df['CTstrict_tumor'] > 0]

print(df_overlaps['Cell'].value_counts())

print(clones['key_0'].tolist())
