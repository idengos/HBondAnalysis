import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt


altered = pd.read_fwf(sys.argv[1])

native = pd.read_fwf(sys.argv[2])


# Watson Crick Edge  (don't worry about the H)
looking_for_WC = ['C23', 'N10', 'C19', 'O18', 'O17', 'N16', 'O20', 'O21', 'C22', 'O6', 'N2', 'N1', 'O6']
# ['O2', 'N3', 'O4']
# Hoogsteen Edge
#looking_for_H = []

# based on which atom
ID = 58 # will be read in with the file name later
for k in range(len(looking_for_WC)):
    looking_for_WC[k] = '37@' + looking_for_WC[k]
    #print(looking_for_WC[k])


# define empty DataFrame for loop
altered_pulled_final = pd.DataFrame()
native_pulled_final = pd.DataFrame()



for i in range(len(looking_for_WC)):
    altered_pulled = altered[altered['#Acceptor'].str.contains(looking_for_WC[i])==True]
    native_pulled = native[native['#Acceptor'].str.contains(looking_for_WC[i])==True]

    altered_pulled_D = altered[altered['Donor'].str.contains(looking_for_WC[i])==True]
    native_pulled_D = native[native['Donor'].str.contains(looking_for_WC[i])==True]





    Acceptor_DF = altered_pulled.drop(['AvgDist','AvgAng'], axis=1)
    Acceptor_DF = Acceptor_DF.rename(columns={'Frames': 'alteredFrames', 'Frac': 'alteredFrac'})
    #print(Acceptor_DF)






    Donor_DF = altered_pulled_D.drop(['AvgDist','AvgAng'], axis=1)
    Donor_DF = Donor_DF.rename(columns={'Frames': 'alteredFrames', 'Frac': 'alteredFrac'})
    Donor_DF["nativeFrames"] = ""
    Donor_DF["nativeFrac"] = ""

    for n in Donor_DF.index:

        aa = Donor_DF['#Acceptor'][n]
        ab = Donor_DF['DonorH'][n]

        ab = ab.partition('@')[2]

        dfplay1 = native_pulled_D[((native_pulled_D['#Acceptor'] == Donor_DF['#Acceptor'][n]) & (native_pulled_D['DonorH'].str.contains(ab)))] # & (native_pulled_D['DonorH'] == Donor_DF['DonorH'][n]))]

        if not dfplay1.empty:
            Donor_DF.at[n, 'nativeFrames'] = dfplay1['Frames'].values[0]
            Donor_DF.at[n, 'nativeFrac'] = dfplay1['Frac'].values[0]
            native_pulled_D = native_pulled_D[((native_pulled_D['#Acceptor'] != dfplay1['#Acceptor'].values[0]) | (native_pulled_D['Frames'] != dfplay1['Frames'].values[0]))]
            #print(native_pulled_D)

        else:
            Donor_DF.at[n, 'nativeFrames'] = 0
            Donor_DF.at[n, 'nativeFrac'] = 0

    native_pulled_D = native_pulled_D.drop(columns=['AvgAng','AvgDist'])
    native_pulled_D = native_pulled_D.rename(columns={'Frames': 'nativeFrames', 'Frac': 'nativeFrac'})
    native_pulled_D['alteredFrac'] = 0
    native_pulled_D['alteredFrames'] = 0
    native_AA = pd.concat([Donor_DF, native_pulled_D]) # $$$$$$$$$$$$$$$$

    if not native_AA.empty:
        native_AA.plot(x='#Acceptor', y=['alteredFrac', 'nativeFrac'], kind='bar')
        plt.title(looking_for_WC[i] + " Donors")


plt.show()








# Series.str.startswith(‘Nov’)
