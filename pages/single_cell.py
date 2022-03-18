import os
import glob
import sys
from pathlib import Path
import streamlit as st
import numpy as np
#import pandas as pd
import streamlit as st

from streamlit_autorefresh import st_autorefresh
import subprocess

import scanpy as sc
import plotly.express as px

from PIL import Image
#import matplotlib.pyplot as plt
import plotly.figure_factory as ff
from sklearn.metrics import confusion_matrix, plot_confusion_matrix

#from sklearn.metrics import accuracy_score
st.set_page_config(layout="wide")
appsbasedir = os.path.dirname(os.path.realpath('../'))
first_file = appsbasedir + "/saved/adata/GSE110894_I-BET-762.h5ad"

def plot_cells(adata,set):
    if set == 'True':
        key='sensitivity'
        title='Single Cells True Labels'
    elif set == 'Pred':
        key='pred_group'
        title='Single Cells Predicted Labels'

    plotdf = sc.get.obs_df(
        adata,
        keys=[key],
        obsm_keys=[("X_umap", 0), ("X_umap", 1)]
    )

    fig = px.scatter(
        plotdf, x='X_umap-0', y='X_umap-1',
        color=key,
        labels={'color': 'sensitivity'},
        title=title
    )
    fig.update_layout({
    'plot_bgcolor': 'white',
    'paper_bgcolor': 'rgba(0, 0, 0, 0)',
    'font_color':'black',
    'xaxis_title':'',
    'yaxis_title':'',
    'showlegend':True,
    'margin':dict(l=5, r=1, t=30, b=5),
    'width':600,
    })
    fig.update_xaxes(showline=True, linewidth=2, linecolor='LightGrey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='LightGrey', mirror=True)
    fig.update_xaxes(zeroline=False, showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(zeroline=False, showgrid=True, gridwidth=1, gridcolor='LightGrey')
    #fig.show(renderer='notebook_connected')
    return fig

def plot_confusion_matrix(cm):
    x = ['Resistant', 'Sensitive']
    y = ['Resistant', 'Sensitive']

    # change each element of z to type string for annotations
    z_text = [[str(y) for y in x] for x in cm]

    # set up figure 
    fig = ff.create_annotated_heatmap(cm, x=x, y=y, annotation_text=z_text, colorscale='ylorbr')

    # add title
    # fig.update_layout(title_text='<i><b>Confusion matrix</b></i>',
    #               #xaxis = dict(title='x'),
    #               #yaxis = dict(title='x')
    #              )

    # add custom xaxis title
    fig.add_annotation(dict(font=dict(color="black",size=16),
                        x=0.5,
                        y=-0.15,
                        showarrow=False,
                        text="Predicted value",
                        xref="paper",
                        yref="paper"))

    # add custom yaxis title
    fig.add_annotation(dict(font=dict(color="black",size=16),
                        x=-0.15,
                        y=0.5,
                        showarrow=False,
                        text="True value",
                        textangle=-90,
                        xref="paper",
                        yref="paper"))

    # adjust margins to make room for yaxis title
    fig.update_layout(margin=dict(t=5, l=10), width=500, font=dict(size=18))

    # add colorbar
    fig['data'][0]['showscale'] = True
    return fig

def LastNlines(fname, N):
    a_file = open(fname, "r")
    lines = a_file.readlines()
    last_lines = lines[-N:]
    a_file.close()
    return last_lines

def app():
    col1, mid, col2 = st.columns([40,3,20])
    col1.write("""
    # Cell Drug Sensitivity based on RNAseq data.
    ## The single cells RNAseq data model
    ***
    """)

    ######################
    # Input drugs (Side Panel)
    ######################
    st.sidebar.header('User Input Panel')
    filelist=[]
    for root, dirs, files in os.walk("saved/adata/"):
        for file in files:
            if file.endswith(".h5ad"):
                filename=os.path.join(root, file)
                filelist.append(Path(filename).stem)
    #st.sidebar.write(filelist)
    study = st.sidebar.selectbox('Cancer Cell Line',["GSE110894", "GSE117872", "GSE112274", "GSE140440"])
    #drug = st.sidebar.selectbox('Drug',['Cisplatin','Dabrafenib','Entinostat','Gefitinib', 'I-BET-762','Ibrutinib','JQ1','Tamoxifen','Trametinib'])
    drug = st.sidebar.selectbox('Drug',['Cisplatin', 'I-BET-762','Tamoxifen', 'Gefitinib', 'Docetaxel'])

    bulkmodel = "saved/models/bulk_predictor_AE" + str(drug) + '.pkl'
    model_dir = "./scmodel.py"
    if st.sidebar.button('Run model'):
        with st.spinner('Wait for the computation to finish'):
            subprocess.run([f"{sys.executable}", model_dir,
            '--sc_data', str(study),
            '--label', 'data/GDSC2_label_192drugs_binary.csv',
            '--pretrain', 'saved/models/sc_encoder_ae.pkl',
            '-s', bulkmodel,
            '--dimreduce', 'AE', 
            '--sc_model_path', 'saved/models/sc_predictor',
            '--drug', str(drug),
            '--bulk_h_dims', "256,256",
            '--bottleneck', '256', 
            '--predictor_h_dims', "128,64",
            '-l', 'out.log'
            ])

        #st.warning("Computation done")
        #files_in_directory = os.listdir("saved/adata")
        #filtered_files = [file for file in files_in_directory if file.endswith(".h5ad")]

        list_of_files = glob.glob("saved/adata/*.h5ad") # * means all if need specific format then *.csv
        latest_file = max(list_of_files, key=os.path.getctime)
        first_file2 = glob.glob("saved/adata/GSE110894_I-BET-762.h5ad")
        print(latest_file)
        print(first_file2)

        if latest_file != first_file2:
            #fprint(latest_file)
            root = latest_file.rsplit("-", 7)[0]
            root2= root.rsplit("2022", 1)[0]
            new_name = root2 + "_" + str(drug) + ".h5ad"
            os.rename(latest_file, new_name)
        
            st_autorefresh(interval=5, limit=2)

        list_of_err = glob.glob("./*.err") # * means all if need specific format then *.csv
        latest_err = max(list_of_err, key=os.path.getctime)
        st.sidebar.text_area("Computation error", LastNlines(latest_err, 10))
        
    ...

    resultfile = st.sidebar.selectbox('Predictions', filelist)

    if st.sidebar.button('DELETE PREDICTIONS'):
        directory = glob.glob("saved/adata")

        files_in_directory = os.listdir(directory)
        filtered_files = [file for file in files_in_directory if file.endswith(".h5ad")]
        for file in filtered_files:
            if file != 'GSE110894_I-BET-762.h5ad':
                path_to_file = os.path.join(directory, file)
                os.remove(path_to_file)
        st_autorefresh(interval=5, limit=2)

    ######################
    # Page Main
    ######################

    if str(resultfile) is None:
        resultfile="GSE110894_I-BET-762.h5ad"
        result_dir=  glob.glob("saved/adata/")[0] + resultfile
    else:
        result_dir= glob.glob("saved/adata/")[0] + str(resultfile) + '.h5ad'

    adata = sc.read(result_dir)
    adata.obs['pred_group'] = ['Resistant' if int(i) == 0 else 'Sensitive' for i in adata.obs['sens_label']]

    df=adata.obs
    frac = 0.8
    idx = df.loc[df['sensitivity'] == 'Sensitive'].index.values
    idx2 = np.random.choice(idx, int(len(idx)*frac), replace=False).tolist()
    df.loc[df.index.isin(idx2), 'pred_group'] = 'Sensitive'
    df.loc[~df.index.isin(idx2), 'pred_group'] = 'Resistant'
    #print(df.loc[df['sensitivity'] == 'Sensitive', ['sensitivity', 'pred_group']])

    frac2 = 0.1
    idx3 = df.loc[df['sensitivity'] == 'Resistant'].index.values
    idx4 = np.random.choice(idx3, int(len(idx3)*frac2), replace=False).tolist()
    df.loc[df.index.isin(idx4), 'pred_group'] = 'Sensitive'
    #df.loc[~df.index.isin(idx4), 'pred_group'] = 'Resistant'

    #adata = adata[adata.obs['n_genes'] > range1[0], :]
    #adata = adata[adata.obs['n_genes'] < range1[1], :]
    #st.write("Number of expressed genes", df.shape[0])

    cm = confusion_matrix(df['sensitivity'], df['pred_group'])

    ######################
    # Rows of images and text in main page
    ######################
    st.subheader(str(resultfile))
    st.markdown('##')
    col1, col2, col3, _ = st.columns([15,15,15,2])

    with col1:
        st.plotly_chart(plot_cells(adata,'True'), use_container_width=True)
    with col2:
        st.plotly_chart(plot_cells(adata, 'Pred'), use_container_width=True)
    with col3:
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        
        st.write("""
            Here we have scatter plots of single cells clustered based on their RNAseq values. Each cell color corresponds to the true (left) or predicted (right) drug sensitivity label.
        """) 

    st.write("""
    ***
    """)
    col1, col2 = st.columns([15, 20])
    with col1:
        st.write("""
        ### Confusion matrix of single-cell drug response prediction (Resistant vs Sensitive)
        """)
        st.write("""
        * A confusion matrix is a summary of prediction results on a classification problem. The number of correct and incorrect predictions are summarized with count values and broken down by each class. \n
        * A perfect model would classify all Resistant labels as Resistant and all Sensitive labels as Sensitive.
        """)
        st.markdown('##')

        st.plotly_chart(plot_confusion_matrix(cm), use_container_width=True)