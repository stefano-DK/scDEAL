import os
import glob
import gdown
from pathlib import Path
import streamlit as st
import numpy as np
import pandas as pd
from streamlit_autorefresh import st_autorefresh
import subprocess
import seaborn as sns
import io
from pdf2image import convert_from_path
import scanpy as sc

from PIL import Image
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
from sklearn.metrics import confusion_matrix, plot_confusion_matrix
#import plotly.express as px
#from sklearn.metrics import accuracy_score
st.set_page_config(layout="wide")
appsbasedir = os.path.dirname(os.path.realpath(__file__))

stdout = io.StringIO()
st.set_option('deprecation.showPyplotGlobalUse', False)

# data download
if not os.path.exists('./data'):
    url1 = 'https://drive.google.com/drive/folders/1AH--DDw7_rDDySYlLEo_VgslIbZhFD_Q?usp=sharing'
    url2 = 'https://drive.google.com/drive/folders/1kTH2hQVNwlXQeV5hCMLVc4eqhrvvag7Q?usp=sharing'
    url3 = 'https://drive.google.com/drive/folders/1eDHumsD3Cbjd9mVkmnxjbEjNq3qXyg4l?usp=sharing'
    gdown.download_folder(url1, output='./data/other', quiet=False, use_cookies=False)
    gdown.download_folder(url2, output='./data/GSE117872', quiet=False, use_cookies=False)
    gdown.download_folder(url3, output='./data/GSE110894', quiet=False, use_cookies=False)

    subprocess.Popen(['mv data/other/*.* data && rm -r data/other'], shell=True)

if not os.path.exists('./saved/adata/GSE1108942022-02-24-10-55-33_I-BET-762.h5ad'):
    gdown.download("https://drive.google.com/file/d/1jZj25AEaeE2kCwwbhqkPgm_Qr4miHztN/view?usp=sharing", output='./saved/adata/GSE1108942022-02-24-10-55-33_I-BET-762.h5ad', quiet=False, use_cookies=False, fuzzy=True)

if not os.path.exists('./saved/models'):
    gdown.download_folder('https://drive.google.com/drive/folders/1HOldnGZ6ZL46bRej933AXf6JMcfV5Tjh?usp=sharing', output='./saved/models', quiet=False, use_cookies=False)


def plot_confusion_matrix(cm):
    x = ['Resistant', 'Sensitive']
    y = ['Sensitive', 'Resistant']

    # change each element of z to type string for annotations
    z_text = [[str(y) for y in x] for x in cm]

    # set up figure 
    fig = ff.create_annotated_heatmap(cm, x=x, y=y, annotation_text=z_text, colorscale='ylorbr')

    # add title
    fig.update_layout(title_text='<i><b>Confusion matrix</b></i>',
                  #xaxis = dict(title='x'),
                  #yaxis = dict(title='x')
                 )

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
                        x=-0.35,
                        y=0.5,
                        showarrow=False,
                        text="True value",
                        textangle=-90,
                        xref="paper",
                        yref="paper"))

    # adjust margins to make room for yaxis title
    fig.update_layout(margin=dict(t=50, l=200))

    # add colorbar
    fig['data'][0]['showscale'] = True
    return fig

######################
# Page Title
######################

st.write("""
# Cell Drug Sensitivity Web App
## This app is based on scDEAL, a Transfer Learning based app.

This app is based on the publication: *"Deep Transfer Learning of Drug Responses by Integrating Bulk and Single-cell RNAseq
data"* and the tool scDEAL. \n 
This app uses Deep Learning models trained on gene expression data from bulk cell lines treated with different 
drugs to learn the relationship between expression and drug exposure. This information is then transferred to train a Deep Learning model
that predicts the drug response of single cells based on their scRNA-seq data.

Predicts the drug sensitivity of single cells! Steps:
1. Select a drug and cancer cell line from the drop down menus;
2. Click on Run button to start the single-cell model training;
3. Explore model performance and predictions via the plots in the main area.

Info on available cancer cell lines:
* **GSE117872**: Cisplatin treated oral squamous cell carcinoma
* **GSE110894**: I-BET-762 treated acute myeloid leukemia
***
""")

######################
# Input drugs (Side Panel)
######################

st.sidebar.header('User Input Features')
filelist=[]
for root, dirs, files in os.walk("saved/adata/"):
      for file in files:
          if file.endswith(".h5ad"):
             filename=os.path.join(root, file)
             filelist.append(Path(filename).stem)
#st.sidebar.write(filelist)
study = st.sidebar.selectbox('Cancer Cell Line',["GSE110894", "GSE117872"])
drug = st.sidebar.selectbox('Drug',['Cisplatin','Dabrafenib','Entinostat','Gefitinib', 'I-BET-762','Ibrutinib','JQ1','Tamoxifen','Trametinib'])


bulkmodel = "saved/models/bulk_predictor_AE" + str(drug) + '.pkl'
model_dir = "./scmodel.py"
if st.sidebar.button('Run model'):
    subprocess.call(['python', model_dir,
    '--sc_data', str(study),
    '--pretrain', 'saved/models/sc_encoder_ae.pkl',
    '-s', bulkmodel,
    '--dimreduce', 'AE', 
    '--sc_model_path', 'saved/models/sc_predictor',
    '--drug', str(drug),
    '--bulk_h_dims', "256,256",
    '--bottleneck', '256', 
    '--predictor_h_dims', "128,64",
#    '-l', 'out.log'
    ])
...
#range1 = st.sidebar.slider('N genes', 0, 5000,(0,5000))
#st.write('Range:', range1)
#st.write('Drug', drug)
 
filelist2=[]
for root, dirs, files in os.walk("saved/figures"):
      for file in files:
             if file.endswith("_roc.pdf"):
                    filename=os.path.join(root, file)
                    filelist2.append(Path(filename).stem)

resultfile2 = st.sidebar.selectbox('Model performances', filelist2)

resultfile = st.sidebar.selectbox('Predictions', filelist)

if st.sidebar.button('REFRESH'):
    st_autorefresh(interval=1, limit=1)
    
    list_of_files = glob.glob(appsbasedir + "/saved/adata/*.h5ad") # * means all if need specific format then *.csv
    latest_file = max(list_of_files, key=os.path.getctime)
    new_name = os.path.splitext(latest_file)[0] + "_" + str(drug) + ".h5ad"
    os.rename(latest_file, new_name)


######################
# Page Title
######################

resultsfolder = appsbasedir + "/saved/adata/"

if str(resultfile) is None:
    resultfile="GSE1108942022-02-24-10-55-33_I-BET-762.h5ad"
    result_dir=appsbasedir + "/saved/adata/" + resultfile
else:
    result_dir=appsbasedir + "/saved/adata/" + str(resultfile) + '.h5ad'

adata = sc.read(result_dir)
adata.obs['pred_groups'] = ['Resistant' if int(i) == 0 else 'Sensitive' for i in adata.obs['sens_label']]

#adata = adata[adata.obs['n_genes'] > range1[0], :]
#adata = adata[adata.obs['n_genes'] < range1[1], :]
st.write("Number of expressed genes", adata.obs.shape[0])

cm = confusion_matrix(adata.obs['sensitivity'], adata.obs['pred_groups'])
#print(cm)
fig, ax = plt.subplots(figsize=(2,1.5))
plt.figure(figsize=(10,10))

sns.set(font_scale=0.3)
sns.heatmap(cm, ax=ax, annot=True, center=0, cmap="YlGnBu")


if str(resultfile2) is None:
    resultfile2  = "AEEntinostat2022-02-22-15-00-38_roc"
    image_dir = appsbasedir + "/saved/figures/" + resultfile2 + '.jpg'
else:
    image_dir = appsbasedir + "/saved/figures/" + str(resultfile2) + '.jpg'
    pdf = appsbasedir + "/saved/figures/" + str(resultfile2) + '.pdf'
    
    pages = convert_from_path(pdf, 500)
    
    for page in pages:
        page.save(image_dir, 'JPEG')

image = Image.open(image_dir)

######################
# First row of images and text
######################

col1, mid, col2 = st.columns([10,1,20])
with col1:
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.write("""
        ### Heatmap of the confusion matrix of model predictions compared to true labels
        """)
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')    
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.write("""
        ### ROC curve showing how the bulk model performs. At the top is the accuracy of the model.
        """)

with col2:
    st.write(plot_confusion_matrix(cm), width=5)
    st.image(image, width=500)