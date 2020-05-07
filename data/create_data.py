import numpy as np
import pandas as pd

data_node = pd.read_csv('nodes_new.txt', sep=' ')
data_triangle = pd.read_csv('triangles_new.txt', sep=' ')

dmerge = data_triangle.merge(data_node, on='NODE')
dmerge = dmerge.sort_values(by='TRIANGLE')

dmerge.to_csv('data_triangles.txt', sep=' ', index=False)

