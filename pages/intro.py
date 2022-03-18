import streamlit as st
from PIL import Image

def app():
    col1, _, col3, _ = st.columns([40,3,20, 20])
    with col1:
        st.write("""
        # Cell Drug Sensitivity based on RNAseq data.
        ### scRNA-seq plays an important role in various aspects of tumor research. It reveals the heterogeneity of tumor cells and monitors the progress of tumor development, thereby preventing further cellular deterioration
        ***
        """)

        st.write("""
        This app uses Deep Learning models trained on gene expression data from bulk cell lines treated with different 
        drugs to learn the relationship between expression and drug exposure. This information is then transferred to train a Deep Learning model
        that predicts the drug response of single cells based on their scRNA-seq data.

        Predicts the drug sensitivity of single cells! Steps:
        1. **Select a drug from the drop down menu and train a bulk cells drug response predictor;**
        2. **Explore the model performance on the dedicated page;**
        3. **Select a cancer cell line and drug to train a model to predict single-cell drug sensitivity based on scRNAseq data;**
        4. **Explore model performance and predictions via the plots in the main area.**

        Info on available cancer cell lines:
        * **GSE117872**: Cisplatin treated oral squamous cell carcinoma
        * **GSE110894**: I-BET-762 treated acute myeloid leukemia
        * **GSE140440**: Docetaxel treated prostate cancer
        * **GSE112274**: Gefitinib treated lung cancer
        ***
        """)

    with col3:
        tissue = Image.open("figures/tissue.jpg")
        cells = Image.open("figures/single_cells.jpeg")
        st.subheader("From bulk cells ...")
        st.image(tissue, width=250)
        st.subheader("... to single cells")
        st.image(cells, width=250)


    col1, col2, col3, _ = st.columns([18,10,30,10])
    with col1:
        st.write("""
    #### 1. The goal is to predict the drug response (sensitive vs resistant) of single cancer cells.

    However you don't have previous drug exposure data, but you do have gene expression data.
    """)
        cancer = Image.open("figures/cancer_cells.jpeg")
        st.image(cancer, width=400)
    
    with col2:
        arrow = Image.open("figures/right-arrow.png")
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.image(arrow, width=150)

    with col3:
        st.write("""
    #### 2. You have available from experiments with cancer cell lines both drug sensitivity and gene expression data.

    Cancer cell lines offer abundance of data of all kind including omics data, drug sensitivity data and others. These can be used to train deep learning and machine learning models.
    """)
        colon = Image.open("figures/cell-culture-plates.tiff")
        st.markdown('##')
        st.image(colon, width=400)
    
    _, _, _, col4 = st.columns([5,5,8,20])
    with col4:
        down_arrow = Image.open("figures/down-arrow2.png")
        st.image(down_arrow, width=120)
    
    col1, col2, col3, col4 = st.columns([18,10,30,10])
    with col1:
        transfer_learning = Image.open("figures/transfer_learning.png")
        st.write("""
    #### 4. Transfer parameter values the tissue model to a single cells model.
    Once a tissue cells drug sensitivity predictor model has been trained, one can transfer its parameter values (weights) to a model that will predict the sensitivity of single cells. 
    
    This new model is trained on single cells gene expression data.
    """)
        st.image(transfer_learning, width=400)

    with col2:
        left_arrow = Image.open("figures/left-arrow.png")
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.markdown('##')
        st.image(left_arrow, width=150)


    with col3:
        NN = Image.open("figures/NN.jpg")
        st.write("""
    #### 3. Train a tissue data drug sensitivity predictor.
    With the available gene expression and drug senstivity data, you can train a Deep Learning model (Neural Network) 
    
    to predict the drug response based on the gene expression pattern of tissue cells.
    """)
        st.markdown('##')
        st.markdown('##')
        st.image(NN, width=500)

    _, col2, _, _ = st.columns([3,5,20,20])
    with col2:
        down_arrow = Image.open("figures/down-arrow2.png")
        st.markdown('##')
        st.markdown('##')
        st.image(down_arrow, width=150)


    col1, _, _, _ = st.columns([18,10,30,10])
    with col1:
        transfer_learning = Image.open("figures/transfer_learning.png")
        st.write("""
    #### 5. Use values from any single cell RNAseq dataset to predict drug response.
    Once we have a single-cell drug sensitivity predictor, we have use gene expression data from any scRNAseq experiment to predict drug response.

    Bare in mind, the prediction is drug-specific!
    """)
        gene_expression = Image.open("figures/gene-expression.png")
        st.markdown('##')
        st.image(gene_expression, width=500)