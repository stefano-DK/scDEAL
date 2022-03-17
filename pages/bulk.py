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
appsbasedir = os.path.dirname(os.path.realpath('../'))

def app():
    col1, mid, col2 = st.columns([40,3,20])
    col1.write("""
    # Cell Drug Sensitivity based on RNAseq data.
    ## The bulk RNAseq data model
    ***
    """)
    study = st.sidebar.selectbox('Cancer Cell Line',["GSE110894", "GSE117872", "GSE112274", "GSE140440"])
    drug = st.sidebar.selectbox('Drug',['Cisplatin', 'I-BET-762','Tamoxifen', 'Gefitinib', 'Docetaxel'])

    bulkmodel = "saved/models/bulk_predictor_AE" + str(drug) + '.pkl'
    model_dir = "./bulkmodel.py"
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
    
    image_dir = "saved/figures/" + str(drug) +'.jpg'

    image = Image.open(image_dir)

    col1, col2, col3, col4 = st.columns([20, 5,20,10])
    with col1:
        # st.markdown('##')
        # st.markdown('##')
        # st.markdown('##')
        st.write("""
            ### ROC curve with AUC value (top) of:
            """)
        st.subheader(str(drug))

        st.markdown('##')
        st.write("""
            An ROC curve (receiver operating characteristic curve) is a graph showing the performance of a classification model at all classification thresholds. 
            This curve plots two variables:
            * True Positive Rate (TPR)
            * False Positive Rate (FPR)

            The bigger the Area Under the Curve (AUC) is, the better the model predicts. Usually one picks the threshold value on the curve providing the highest TPR with the lowest FPR.

            The blue dashed line represents the baseline, or the TPR/FPR values given by random classifier.
            """)
        
        st.markdown('##')
        st.image(image, width=500)
