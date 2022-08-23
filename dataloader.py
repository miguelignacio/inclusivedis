import os
import uproot3
import pandas as pd
import numpy as np

def get_Dataframe(path, name='Data', tag=None, verbose=False):
    Files = os.listdir(path) 
    #print('you have entered get_Dataframe')
    #print ('FILES ' , Files)
    df = None
    for i, f in enumerate(Files):
        #if( df is not None):
        #    if(df.shape[0]>5000000): continue
        #if(i>10):continue
        if name not in f: continue
        filename = path+f
        if not(tag is None) and (tag not in f): continue
        if (verbose):
            print ('filename is' , filename)
        
        temp_file = uproot3.open(filename)
        
        hasTree = False 

        if (verbose):
            print (temp_file.keys()) 
                    
        if(len(temp_file.keys())<1):
            if (verbose):
                print('could not find %s, skipping'%name)
            continue
        
        if( not(name in str(temp_file.keys()[0]))):
            if (verbose):
                print('could not find %s, skipping'%name)
            continue
        
        for key in temp_file[name].keys():
            #print (key)
            if('minitree' in str(key)):
                hasTree=True
        if (not hasTree):
            if (verbose):
                print('file has not minitree, skipping')
            continue

        temp_tree = temp_file[name+'/minitree']

        #print('Variables in tree ' , temp_tree.keys())
        temp_df = None
        
        if 'Data' not in name:
            try:
                temp_df   =  temp_tree.pandas.df(['gen_event*','event_*','gene*','e_*','genHFS*','HFS*','wgt','vertex_z','ptmiss','Empz'], entrystop=3e7,flatten=True)
                df = pd.concat([df,temp_df])
            except ValueError:
                if (verbose):
                    print ('oops, there is a problem in flattening the TTree ')
        else:
            try:
                temp_df   =  temp_tree.pandas.df(['event_*','e_*','HFS*','wgt','vertex_z','ptmiss','Empz'], entrystop=3e7,flatten=True) 
                df = pd.concat([df,temp_df])
            except ValueError:
                if (verbose):
                    print ('oops, there is a problem in flattening the TTree ')
        
        #try:
        #    df.shape[0]
        #except ValueError:
        #    print('no valid dataframe')
    if (verbose):
        print('####################################################################')
        if( not(df is None)):
            print('Dataframe has a total of ', df.shape[0], ' entries')
        else:
            print ('Dataframe has no entry, it is None')
        print('####################################################################')

    return df

def applyCut(dataframe, cut, text=None,verbose=False):
    #dataframe = inputDataframe
    nbeforecut = dataframe.shape[0]
    dataframe = dataframe.query(cut)
    if text and verbose:
        print (text, nbeforecut, ' fraction kept: %2.1f'%(100.0*float(dataframe.shape[0])/nbeforecut))
    return dataframe

def applyCuts(df,isMC=False,verbose=False):
    temp = df

    temp['pass_reco'] = np.where(temp['event_Q2_e']>0, 1, 0)
    if (isMC):
        temp.eval('gene_e = sqrt(gene_px*gene_px + gene_py*gene_py +gene_pz*gene_pz)',inplace=True)

        temp['pass_truth'] = np.where(temp['gen_event_Q2_e']>0, 1, 0)
        temp['pass_fiducial'] = np.where(temp['pass_truth']*(temp['gene_e'] > 12),
                                          1, 0)
        
    #temp = applyCut(temp, 'abs(vertex_z)<25 and vertex_z!=0','abs(vertex_z)<25 and and vertex_z!=0')
    #temp = applyCut(temp, 'tau1b>0 and tau1b<1', '0<tau1b<1')


    temp.eval('e_pt = sqrt(e_px*e_px + e_py*e_py)',inplace=True)
    temp.eval('e_e = sqrt(e_px*e_px + e_py*e_py +e_pz*e_pz)',inplace=True)

    temp.eval('e_phi = arctan(e_py/e_px)', inplace=True)

    temp = applyCut(temp, 'pass_reco==0 | ptmiss < 10', 'ptmiss<10',verbose)

    #temp = applyCut(temp, 'pass_reco==0 | 0.08 < event_y_es < 0.7', '0.08 < event_y_es < 0.7',verbose)
    #temp = applyCut(temp, 'pass_reco==0 | event_Q2_es>150', 'event_Q2_es>150',verbose)
   # temp = applyCut(temp, 'pass_reco==0 | Q2<10000', 'Q2<10000')
    temp = applyCut(temp, 'pass_reco==0 | Empz<65', 'Empz<65',verbose)
    temp = applyCut(temp, 'pass_reco==0 | Empz>45', 'Empz>45',verbose)


    if(isMC):
        temp = applyCut(temp,'pass_truth>0',' pass_truth>0',verbose)


        temp.eval('gene_pt = sqrt(gene_px*gene_px + gene_py*gene_py)',inplace=True)

        temp.eval('gene_phi = arctan(gene_py/gene_px)', inplace=True)

        


    #Save only the features we need.
    #if (isMC):
    #    temp = temp[['gene_px','gene_py','gene_pz','e_px','e_py','e_pz','HFS_px','HFS_py','HFS_pz','HFS_E','HFS_eta',
    #                 'genHFS_px','genHFS_py','genHFS_pz','genHFS_E','genHFS_eta',                
    #                 'wgt','pass_reco','pass_truth', 'pass_fiducial']]
    #else:
    #    temp = temp[['e_px','e_py','e_pz','HFS_px','HFS_py','HFS_pz','HFS_E','HFS_eta','wgt','pass_reco']]
        
    #df = applyCut(df, 'n_total>1', ' n>1')
    return temp
