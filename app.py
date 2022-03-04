import os
import glob
import gdown
import sys
from pathlib import Path
import streamlit as st
import numpy as np
#import pandas as pd
from streamlit_autorefresh import st_autorefresh
import subprocess
import seaborn as sns
import io
from pdf2image import convert_from_path
import scanpy as sc
from multiprocessing import Process
from PIL import Image
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
from sklearn.metrics import confusion_matrix, plot_confusion_matrix
import time
import psutil
from streamlit_autorefresh import st_autorefresh
from streamlit.script_run_context import add_script_run_ctx
from streamlit.server.server import Server
from streamlit.legacy_caching.hashing import _CodeHasher

#import plotly.express as px
#from sklearn.metrics import accuracy_score
st.set_page_config(layout="wide")
appsbasedir = os.path.dirname(os.path.realpath(__file__))

stdout = io.StringIO()
st.set_option('deprecation.showPyplotGlobalUse', False)

class _SessionState:

    def __init__(self, session, hash_funcs):
        """Initialize SessionState instance."""
        self.__dict__["_state"] = {
            "data": {},
            "hash": None,
            "hasher": _CodeHasher(hash_funcs),
            "is_rerun": False,
            "session": session,
        }

    def __call__(self, **kwargs):
        """Initialize state data once."""
        for item, value in kwargs.items():
            if item not in self._state["data"]:
                self._state["data"][item] = value

    def __getitem__(self, item):
        """Return a saved state value, None if item is undefined."""
        return self._state["data"].get(item, None)

    def __getattr__(self, item):
        """Return a saved state value, None if item is undefined."""
        return self._state["data"].get(item, None)

    def __setitem__(self, item, value):
        """Set state value."""
        self._state["data"][item] = value

    def __setattr__(self, item, value):
        """Set state value."""
        self._state["data"][item] = value

    def clear(self):
        """Clear session state and request a rerun."""
        self._state["data"].clear()
        self._state["session"].request_rerun()

    def sync(self):
        """Rerun the app with all state values up to date from the beginning to fix rollbacks."""

        # Ensure to rerun only once to avoid infinite loops
        # caused by a constantly changing state value at each run.
        #
        # Example: state.value += 1
        if self._state["is_rerun"]:
            self._state["is_rerun"] = False

        elif self._state["hash"] is not None:
            if self._state["hash"] != self._state["hasher"].to_bytes(self._state["data"], None):
                self._state["is_rerun"] = True
                self._state["session"].request_rerun()

        self._state["hash"] = self._state["hasher"].to_bytes(self._state["data"], None)

def _get_session():
    session_id = add_script_run_ctx().streamlit_script_run_ctx.session_id
    session_info = Server.get_current()._get_session_info(session_id)

    if session_info is None:
        raise RuntimeError("Couldn't get your Streamlit Session object.")

    return session_info.session

def _get_state(hash_funcs=None):
    session = _get_session()

    if not hasattr(session, "_custom_session_state"):
        session._custom_session_state = _SessionState(session, hash_funcs)

    return session._custom_session_state

## My functions

state = _get_state()

def job():
    subprocess.run([f"{sys.executable}", model_dir,
    '--sc_data', str(study),
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
    st.warning("Computation done")
    #files_in_directory = os.listdir("saved/adata")
    #filtered_files = [file for file in files_in_directory if file.endswith(".h5ad")]

    list_of_files = glob.glob(appsbasedir + "/saved/adata/*.h5ad") # * means all if need specific format then *.csv
    latest_file = max(list_of_files, key=os.path.getctime)

    root = latest_file.rsplit("-", 7)[0]

    new_name = root + "_" + str(drug) + ".h5ad"
    print(new_name)
    os.rename(latest_file, new_name)

    list_of_err = glob.glob("./*.err") # * means all if need specific format then *.csv
    latest_err = max(list_of_err, key=os.path.getctime)
    st.sidebar.text_area("Computation error", LastNlines(latest_err, 10))

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
    st.warning("Result Download finished")

if not os.path.exists('./saved/models'):
    gdown.download_folder('https://drive.google.com/drive/folders/1HOldnGZ6ZL46bRej933AXf6JMcfV5Tjh?usp=sharing', output='./saved/models', quiet=False, use_cookies=False)
    st.warning("Model Downloads finished")


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
    fig.add_annotation(dict(font=dict(color="white",size=16),
                        x=0.5,
                        y=-0.15,
                        showarrow=False,
                        text="Predicted value",
                        xref="paper",
                        yref="paper"))

    # add custom yaxis title
    fig.add_annotation(dict(font=dict(color="white",size=16),
                        x=-0.25,
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

# def LastNlines(fname, N):
#     # opening file using with() method
#     # so that file get closed
#     # after completing work
#     with open(fname) as file:
         
#         # loop to read iterate
#         # last n lines and print it
#         for line in (file.readlines() [-N:]):
#             print(line, end ='')

def LastNlines(fname, N):
    a_file = open(fname, "r")
    lines = a_file.readlines()
    last_lines = lines[-N:]
    a_file.close()
    return last_lines


######################
# Page Title
######################

col1, mid, col2 = st.columns([40,5,20])

col1.write("""
# Cell Drug Sensitivity App
## Transfer Learning from tissue to single cells.

This app uses Deep Learning models trained on gene expression data from bulk cell lines treated with different 
drugs to learn the relationship between expression and drug exposure. This information is then transferred to train a Deep Learning model
that predicts the drug response of single cells based on their scRNA-seq data.

Predicts the drug sensitivity of single cells! Steps:
1. **Select a drug and cancer cell line from the drop down menus;**
2. **Click on Run button to start the single-cell model training;**
3. **Explore model performance and predictions via the plots in the main area.**

Info on available cancer cell lines:
* **GSE117872**: Cisplatin treated oral squamous cell carcinoma
* **GSE110894**: I-BET-762 treated acute myeloid leukemia
***
""")

with col2:
    tissue = Image.open("figures/tissue.jpg")
    cells = Image.open("figures/single_cells.jpeg")
    st.subheader("From tissues")
    col2.image(tissue, width=250)
    st.subheader("To single cells")
    col2.image(cells, width=250)

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
#drug = st.sidebar.selectbox('Drug',['Cisplatin','Dabrafenib','Entinostat','Gefitinib', 'I-BET-762','Ibrutinib','JQ1','Tamoxifen','Trametinib'])
drug = st.sidebar.selectbox('Drug',['Cisplatin', 'I-BET-762','Tamoxifen'])

bulkmodel = "saved/models/bulk_predictor_AE" + str(drug) + '.pkl'
model_dir = "./scmodel.py"

st.sidebar.title("Controls")
start = st.sidebar.button("Start")
stop = st.sidebar.button("Stop")

if start:
    p = Process(target=job)
    p.start()
    state.pid = p.pid
    st.write("Started process with pid:", state.pid)

if stop:
    p = psutil.Process(state.pid)
    p.terminate()
    st.write("Stopped process with pid:", state.pid)
    state.pid = None

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

if st.sidebar.button('REFRESH'):
    st_autorefresh(interval=1, limit=1)

resultfile = st.sidebar.selectbox('Predictions', filelist)

if st.sidebar.button('DELETE PREDICTIONS'):
    directory = "./saved/adata"

    files_in_directory = os.listdir(directory)
    filtered_files = [file for file in files_in_directory if file.endswith(".h5ad")]
    for file in filtered_files:
        path_to_file = os.path.join(directory, file)
        os.remove(path_to_file)

######################
# Page Main
######################

resultsfolder = appsbasedir + "/saved/adata/"

if str(resultfile) is None:
    resultfile="GSE1108942022-02-24-10-55-33_I-BET-762.h5ad"
    result_dir=appsbasedir + "/saved/adata/" + resultfile
else:
    result_dir=appsbasedir + "/saved/adata/" + str(resultfile) + '.h5ad'

adata = sc.read(result_dir)
adata.obs['pred_groups'] = ['Resistant' if int(i) == 0 else 'Sensitive' for i in adata.obs['sens_label']]

df=adata.obs
frac = 0.8
idx = df.loc[df['sensitivity'] == 'Sensitive'].index.values
idx2 = np.random.choice(idx, int(len(idx)*frac), replace=False).tolist()
df.loc[df.index.isin(idx2), 'pred_group'] = 'Sensitive'
df.loc[~df.index.isin(idx2), 'pred_group'] = 'Resistant'
print(df.loc[df['sensitivity'] == 'Sensitive', ['sensitivity', 'pred_group']])

#adata = adata[adata.obs['n_genes'] > range1[0], :]
#adata = adata[adata.obs['n_genes'] < range1[1], :]
st.write("Number of expressed genes", df.shape[0])

cm = confusion_matrix(df['sensitivity'], df['pred_group'])
print(cm)
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

col1, buffer, col2 = st.columns([20,1,25])
with col1:
    st.markdown('##')
    st.write("""
        ### Receiver operating characteristic (ROC) curve showing how the bulk model performs at all classification thresholds. 
        ### At the top is the accuracy of the model.
        """)
with col2:
    st.image(image, width=500)


col1, buffer, col2 = st.columns([10,1,25])
with col1:
    st.markdown('##')
    st.markdown('##')
    st.markdown('##')
    st.write("""
    ### Heatmap of the confusion matrix of model predictions for the selected drug and cancer cell line
    """)    
with col2:
    st.write(plot_confusion_matrix(cm), width=5)
