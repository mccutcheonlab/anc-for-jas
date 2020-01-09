# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 09:00:27 2020

@author: admin
"""

import pandas as pd

data = [{'a': 20, 'b': 40}]

df1 = pd.DataFrame(data)

data = [{'c': 30, 'd': 60}]

df2 = pd.DataFrame(data)

df = pd.concat([df1, df2], axis=1)


data = {'c': 30, 'd': 60}

for key in data.keys():
    print(key, data[key])