import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sankey

traffic_data = [ [ "amazonia", "bight_of_benin", 1947],
                 [ "amazonia", "bight_of_biafra", 2533],
                 [ "amazonia", "se_africa_ and_indian_ocean_ isls", 1553],
                 [ "amazonia", "senegambia", 83559],
                 [ "amazonia", "sierra_leone", 42],
                 [ "amazonia", "wc_africa_st_helena", 55515],
                 [ "amazonia", "other", 6225],
                 [ "bahia", "bight_of_benin", 722319],
                 [ "bahia", "bight_of_biafra", 67715],
                 [ "bahia", "gold_coast", 10976],
                 [ "bahia", "se_africa_ and_indian_ocean_ isls", 18019],
                 [ "bahia", "senegambia", 5260],
                 [ "bahia", "sierra_leone", 1829],
                 [ "bahia", "wc_africa_st_helena", 455896],
                 [ "bahia", "windward_coast", 2458],
                 [ "bahia", "other", 62842],
                 [ "pernambuco", "asia_e_africa", 490],
                 [ "pernambuco", "bight_of_benin", 66232],
                 [ "pernambuco",        "bight_of_biafra", 23346],
                 [ "pernambuco",        "gold_coast", 28962],
                 [ "pernambuco",        "se_africa_ and_indian_ocean_ isls", 14316],
                 [ "pernambuco",        "senegambia", 4074],
                 [ "pernambuco",        "wc_africa_st_helena", 327033],
                 [ "pernambuco",        "other", 13198],
                 [ "se_brazil_riodjan", "bight_of_benin", 19894],
                 [ "se_brazil_riodjan", "bight_of_biafra", 17576],
                 [ "se_brazil_riodjan", "gold_coast", 3766],
                 [ "se_brazil_riodjan", "se_africa_ and_indian_ocean_ isls", 232821],
                 [ "se_brazil_riodjan", "senegambia", 2264],
                 [ "se_brazil_riodjan", "sierra_leone", 1396],
                 [ "se_brazil_riodjan", "wc_africa_st_helena", 1156304],
                 [ "se_brazil_riodjan", "other", 109097],
                 [ "other", "wc_africa_st_helena", 452]]


df = pd.DataFrame(traffic_data, columns = ['disembark', 'embark', 'count'])

sankey.sankey( left        = df['embark'], 
               right       = df['disembark'], 
               leftWeight = df['count'],
               figureName = "/master/nplatt/sch_man_nwinvasion/results/sankey.svg")

