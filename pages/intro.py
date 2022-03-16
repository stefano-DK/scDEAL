import streamlit as st
from PIL import Image

def app():
    col1, mid, col2 = st.columns([40,3,20])
    col1.write("""
    # Cell Drug Sensitivity based on RNAseq data.
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
    * **GSE140440**: Docetaxel treated prostate cancer
    * **GSE112274**: Gefitinib treated lung cancer
    ***
    """)

    with col2:
        tissue = Image.open("figures/tissue.jpg")
        cells = Image.open("figures/single_cells.jpeg")
        st.subheader("From tissues...")
        col2.image(tissue, width=250)
        st.subheader("... to single cells")
        col2.image(cells, width=250)