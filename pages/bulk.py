import os
import glob
import sys
from pathlib import Path
import streamlit as st
import numpy as np
#import pandas as pd
import streamlit as st
from pdf2image import convert_from_path

from streamlit_autorefresh import st_autorefresh
import subprocess

import scanpy as sc
import plotly.express as px

from PIL import Image
#import matplotlib.pyplot as plt
import plotly.figure_factory as ff
from sklearn.metrics import confusion_matrix, plot_confusion_matrix
appsbasedir = os.path.dirname(os.path.realpath('./scDEAL'))

def app():
    col1, mid, col2 = st.columns([40,3,20])
    col1.write("""
    # Cell Drug Sensitivity based on RNAseq data.
    ## The bulk RNAseq data model
    ***
    """)
    study = st.sidebar.selectbox('Cancer Cell Line',["GSE110894", "GSE117872", "GSE112274", "GSE140440"])
    drug = st.sidebar.selectbox('Drug',['Cisplatin', 'I-BET-762','Tamoxifen', 'Gefitinib', 'Docetaxel'])
    dim_red = st.sidebar.selectbox('Dim. Red.',['AE', 'VAE'])

    bulkmodel = "saved/models/bulk_predictor_"
    model_dir = "./bulkmodel.py"
    if str(drug) == 'I-BET-762':
        drug2 = 'I.BET.762'
    else:
        drug2 = str(drug)

    if st.sidebar.button('Run model'):
        with st.spinner('Wait for the computation to finish'):
            subprocess.run([f"{sys.executable}", model_dir,
            '--drug', drug2,
            '-e', 'saved/models/bulk_encoder_ae_256.pkl',
            '--label', 'data/GDSC2_label_192drugs_binary.csv',
            '-p', bulkmodel,
            '--epochs', '100',
            '--dimreduce', str(dim_red),
            '--encoder_h_dims', "256,256",
            '--predictor_h_dims', "128,64",
            '--bottleneck', '256',
            '-l', 'out.log'
            ])
    
    directory = appsbasedir + "/saved/figures/"

    files_in_directory = os.listdir(directory)
    filtered_files = [file for file in files_in_directory if file.endswith(".pdf")]

    for file in filtered_files:
        path_to_file = os.path.join(directory, file)
        if "2022" in file:
            if "_roc" in file:
                root= file.rsplit("2022", 1)[0]
                image_file = appsbasedir + "/saved/figures/" + root + '_roc.jpg'
                new_name = root + "_roc.pdf"
            if "_prc" in file:
                root= file.rsplit("2022", 1)[0]
                image_file = appsbasedir + "/saved/figures/" + root + '_prc.jpg'
                new_name = root + "_prc.pdf"
            
            path_to_new_file = os.path.join(directory, new_name)
            os.rename(path_to_file, path_to_new_file)
        
            pages = convert_from_path(path_to_new_file, 500)
            
            for page in pages:
                page.save(image_file, 'JPEG')

    st_autorefresh(interval=5, limit=2)

    filelist=[]
    for root, dirs, files in os.walk("saved/figures"):
        for file in files:
            if file.endswith(".jpg"):
                if "_roc" in file:
                    #filename=os.path.join(root, file)
                    filelist.append(file.rsplit("_", 1)[0])
    
    #print(filelist)
    bulk_pred = st.sidebar.selectbox('Bulk Predictors', filelist)
    

    image_dir = "saved/figures/" + str(bulk_pred) +'_roc.jpg'
    image = Image.open(image_dir)

    image_dir2 = "saved/figures/" + str(bulk_pred) +'_prc.jpg'
    image2 = Image.open(image_dir2)

    filelist2=[]
    for root, dirs, files in os.walk("saved/figures"):
        for file in files:
            if file.endswith("prc.jpg"):
                filename=os.path.join(root, file)
                filelist2.append(Path(filename).stem)


    col1, col2 = st.columns([20, 20])
    with col1:
        # st.markdown('##')
        # st.markdown('##')
        # st.markdown('##')
        st.write("""
            ### ROC curve with AUC value (top) of:
            """)
        st.subheader(str(drug))
        st.write("""
            A ROC curve (receiver operating characteristic curve) is a graph showing the performance of a classification model at all classification thresholds. 
            This curve plots two variables:
            * True Positive Rate (TPR)
            * False Positive Rate (FPR)

            The bigger the Area Under the Curve (AUC) is, the better the model predicts. Usually one picks the threshold value on the curve providing the highest TPR with the lowest FPR.

            The blue dashed line represents the baseline, or the TPR/FPR values given by random classifier.
            ***
            """)
    
    col1, col2, col3= st.columns([20, 20, 40])

    with col1:
        st.image(image, width=500)
    with col2:
        st.image(image2, width=500)
