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

def app():
    col1, mid, col2 = st.columns([40,3,20])
    col1.write("""
    # Cell Drug Sensitivity based on RNAseq data.
    ## The bulk RNAseq data model
    ***
    """)
