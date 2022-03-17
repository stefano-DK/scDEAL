import os
import gdown
from pathlib import Path
import streamlit as st
import numpy as np
import subprocess
#import pandas as pd

# Custom imports 
from multipage import MultiPage
from pages import single_cell, bulk, intro

#from sklearn.metrics import accuracy_score

appsbasedir = os.path.dirname(os.path.realpath(__file__))
first_file = appsbasedir + "/saved/adata/GSE110894_I-BET-762.h5ad"

# data download
if not os.path.exists('./data'):
    url1 = 'https://drive.google.com/drive/folders/1AH--DDw7_rDDySYlLEo_VgslIbZhFD_Q?usp=sharing'
    gdown.download_folder(url1, output='./data/other', quiet=False, use_cookies=False)
    subprocess.Popen(['mv data/other/*.* data && rm -r data/other'], shell=True)
if not os.path.exists('./data/GSE117872'):  
    url2 = 'https://drive.google.com/drive/folders/1kTH2hQVNwlXQeV5hCMLVc4eqhrvvag7Q?usp=sharing'
    gdown.download_folder(url2, output='./data/GSE117872', quiet=False, use_cookies=False)
if not os.path.exists('./data/GSE110894'):
    url3 = 'https://drive.google.com/drive/folders/1eDHumsD3Cbjd9mVkmnxjbEjNq3qXyg4l?usp=sharing'
    gdown.download_folder(url3, output='./data/GSE110894', quiet=False, use_cookies=False)
if not os.path.exists('./data/GSE112274'):
    url4 = 'https://drive.google.com/drive/folders/1o7MjpiQ08Kc0DgTt5Obhazl2YnsAQ-k8?usp=sharing'
    gdown.download_folder(url4, output='./data/GSE112274', quiet=False, use_cookies=False)
if not os.path.exists('./data/GSE149383'):
    url5 = 'https://drive.google.com/drive/folders/1XHqaKJt-xQXRADgwCRLHhoGtApHmFfTK?usp=sharing'
    gdown.download_folder(url5, output='./data/GSE149383', quiet=False, use_cookies=False)
if not os.path.exists('./data/GSE140440'):
    url6 = 'https://drive.google.com/drive/folders/1G-LnaDTpJm-jy0IDtcbUb7hU4_7uxDHY?usp=sharing'
    gdown.download_folder(url6, output='./data/GSE140440', quiet=False, use_cookies=False)


if not os.path.exists('./saved/adata/GSE110894_I-BET-762.h5ad'):
    gdown.download("https://drive.google.com/file/d/1jZj25AEaeE2kCwwbhqkPgm_Qr4miHztN/view?usp=sharing", output='./saved/adata/GSE110894_I-BET-762.h5ad', quiet=False, use_cookies=False, fuzzy=True)
    st.warning("Result Download finished")
    
if not os.path.exists('./saved/models'):
    gdown.download_folder('https://drive.google.com/drive/folders/1HOldnGZ6ZL46bRej933AXf6JMcfV5Tjh?usp=sharing', output='./saved/models', quiet=False, use_cookies=False)
    st.warning("Model Downloads finished")


app = MultiPage()

# Add all your application here
app.add_page("INTRO", intro.app)
app.add_page("Single Cells", single_cell.app)
app.add_page("Tissues", bulk.app)

# The main app
app.run()