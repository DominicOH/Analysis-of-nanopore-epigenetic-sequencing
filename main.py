import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import sqrt
from common import *

# consider making the file path an arg for easy switching out
nanopore_ternary_data = readModbam2bed("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/data/cbm1_prom_modbases_mapq60.bed")

##### Load oxBS dataset ##### 
oxbs_df = readBismarkZeroCov("/mnt/data1/doh28/data/CRR008808_oxBS/extraction_output/mapq_filter/CRD018548.gz_val_1_bismark_bt2_pe.deduplicated_mapq10_sorted.bedGraph.gz.bismark.zero.cov",
                             "5mC", True)

##### Load TAB dataset #####
tab_df = readBismarkZeroCov("/mnt/data1/doh28/data/CRR008807_TAB/mapq_filtered/modified_bases/CRR008807_TAB_merged_resorted_q10.bedGraph.gz.bismark.zero.cov",
                            "5hmC", True)
