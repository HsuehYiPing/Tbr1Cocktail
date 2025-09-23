import numpy as np
import pandas as pd
import os
import os.path
import re
import fnmatch
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
mpl.use('Agg')
#%matplotlib inline
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns; sns.set(style="white", context = "paper", font_scale = 0.729, palette='colorblind', color_codes = True, rc = {"font.size":7, "xtick.bottom":True, "ytick.left":True, "xtick.major.size": 2, "ytick.major.size": 2})
from sklearn.metrics import auc
from tqdm import tqdm
from sklearn.feature_selection import f_regression 
from scipy.stats import pearsonr
from itertools import chain
from sklearn.metrics.pairwise import cosine_similarity
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 


dsets = ["NO.1", "NO.2", "NO.3", "NO.8", "NO.13", "NO.15", "NO.16", "NO.20"]
dstatus = ["OE", "RSI"]
dtreatments = ["Water", "Cocktail"]

for dtreatment in dtreatments:
      for dset in dsets:
            for dstat in dstatus:
                  fdir1 = "Calcium_imaging_data/" + dtreatment + "/" + dset + "/" + dstat
                  print(fdir1)
                  outputdir = fdir1 + "/output_files"
                  if not os.path.exists(outputdir): os.mkdir(outputdir)

                  accfile = ""
                  blankfile = ""
                  oefile = ""
                  rsifile = ""
                  
                  for filename in os.listdir(fdir1):
                        if fnmatch.fnmatch(filename, '*OE_behavior_period.csv'): oefile = fdir1 + "/" + filename
                        if fnmatch.fnmatch(filename, '*RSI_behavior_period.csv'): rsifile = fdir1 + "/" + filename
                        if fnmatch.fnmatch(filename, '*acc*.csv'): accfile = fdir1 + "/" + filename 
                        if fnmatch.fnmatch(filename, '*blank*.csv'): blankfile = fdir1 + "/" + filename

                  # z-score transformation of calcium traces

                  #import csv traces file, convert to pd Dataframe and than to numpy array
                  accraw = pd.read_csv(accfile, header = None).values
                  accrawt = np.transpose (accraw)
                  accq = np.quantile (accrawt,0.5, axis =1) # set baseline value to 0.5 quantile for all traces
                  #extract baseline value condition less than 0.5 quantile for all traces
                  condition = []
                  for i in range (len(accrawt)):
                        result = accrawt[i] < accq[i]
                        condition.append (result) 
                  # z-score transformation    
                  accrawt_z = []
                  for x in range (len (accrawt)):
                        mean = np.mean (np.extract (condition[x], accrawt[x]))
                        std = np.std (np.extract (condition[x], accrawt[x]))
                        z = (accrawt[x]-mean)/std
                        accrawt_z.append (z)
                        print(z)
                  #print (np.array(accrawt_z).shape)
                  pd.DataFrame(accrawt_z).to_csv(outputdir + '/accrawt_z.csv')
                  accnpt = np.array(accrawt_z)
                  
                  #add 51 nan elements before and after each calcium traces, np.full(cell number, 51) 
                  accnpt_plus = np.concatenate([np.full([len(accrawt), 51], np.nan), accrawt_z, np.full([len(accrawt), 51], np.nan)], axis = 1)
                  
                  ################            
                  BF = []
                  if (blankfile != ""):
                        # #import blank frame time point csv, and insert nan value into these time point
                        BF = pd.read_csv(blankfile, index_col = 0, header = 0).values
                        # print (BF)
                        accrawt_z2= accrawt_z

                        nans = [np.nan]
                        nans = np.transpose (nans)
                        nan = np.array(nans)

                        for a in BF:
                              print("insert nan in time: ", a)
                              accrawt_z2 = np.insert(accrawt_z2, a-1, nan, 1)
                              print (np.array(accrawt_z2).shape)
                        pd.DataFrame(accrawt_z2).to_csv(outputdir + '/accrawt_z_6000.csv')
                        accnpt2 = np.array(accrawt_z2)
                        
                        #add 51 nan elements before and after each calcium traces, np.full(cell number, 51) 
                        accnpt_plus = np.concatenate([np.full([len(accrawt), 51], np.nan), accrawt_z2, np.full([len(accrawt), 51], np.nan)], axis = 1)
                        
                        
                  ################

                  # Generate shuffled intaction time (from 6000 sample points) from interaction time list function            
                  def random_interaction_interaction_mask(interaction_list, interaction_list_1d, maxi, frame_loss):
      
                        '''
                        interaction_list: 2d interaction frame list
                        interaction_list_1d: 1d interaction frame list
                        maxi: maximum of recording frame
                        
                        '''
                        
                        #add maximum value
                        new_interaction_list= interaction_list+ [[maxi]]
                        #     print(new_interaction_list)

                        #create non-interaction list (2d)
                        non_interaction_list= []
                        init= 0
                        for i in new_interaction_list:
                              if init < min(i):
                                    if init not in i:
                                          f = np.arange(init, min(i), 1)
                                          init = max(i)+1
                                          non_interaction_list.append(f.tolist())
                        #     print(non_interaction_list)

                        #combine two list
                        entire_list = interaction_list+ non_interaction_list
                        
                        #randomize 2d list
                        random_list = np.random.permutation(entire_list)
                        
                        #convert 2d list to 1d list
                        new_random_list1 = list(chain(*random_list))
                        
                        #remove frame loss
                        new_random_list2 = [ele for ele in new_random_list1 if ele not in frame_loss]
                        # print(new_random_list2)
                        
                        #convert random list to binary mask
                        random_mask= []
                        for j in new_random_list2:
                              if j in interaction_list_1d:
                                    random_mask.append(1)
                              else:
                                    random_mask.append(0)
                        return random_mask
                  
                  # Calculate the auROC of activity during interaction and non-interaction period
                  def auROC (activity, interaction_mask, number, fdir):
                        
                        df_activity = pd.DataFrame(list(zip(interaction_mask, activity)), columns = ["mask", "value"])

                        mask1_activity = df_activity[df_activity[["mask"]].isin([1]).any(axis=1)]
                        mask0_activity = df_activity[df_activity[["mask"]].isin([0]).any(axis=1)]
                        
                        cell_tracesr = mask1_activity["value"].to_numpy().reshape(-1,1)
                        baseliner = mask0_activity["value"].to_numpy().reshape(-1,1)
                        criterionr = np.linspace (min(activity), max (activity), 100)
                        
                        TPRr = []
                        FPRr = []
                        for ctprr in range (len(criterionr)):
                              tprr = np.sum (cell_tracesr> criterionr[ctprr])/len (cell_tracesr)
                              TPRr.append (tprr)
                        for cfprr in range (len(criterionr)):
                              fprr = np.sum (baseliner > criterionr[cfprr])/len (baseliner)
                              FPRr.append (fprr)
                        arear = auc (FPRr, TPRr)
                        RSr = (arear-0.5)*2     #test statistic : calcium signal auROC generated RS
                        
                        #export auROC_value of each cell into CSV files    
                        pd.DataFrame(cell_tracesr).to_csv(fdir + '/cell_auROC_csv/cell' + str(number+1) + '_mask1_activity.csv')
                        pd.DataFrame(baseliner).to_csv(fdir + '/cell_auROC_csv/cell' + str(number+1) + '_mask0_activity.csv')
                        pd.DataFrame(TPRr).to_csv(fdir + '/cell_auROC_csv/cell' + str(number+1) + '_TPRr.csv')
                        pd.DataFrame(FPRr).to_csv(fdir + '/cell_auROC_csv/cell' + str(number+1) + '_FPRr.csv')
                        pd.DataFrame([RSr]).to_csv(fdir + '/cell_auROC_csv/cell' + str(number+1) + '_auROC.csv')


                  ################
                  def dataprocess (fname):
                        # print(setname)
                        fdir = outputdir

                        if os.path.exists(fname):
                              #import interaction event time csv, and values to list
                              activeint = pd.read_csv(fname, index_col = None, header = None)
                              # print(activeint)
                              # print(len(activeint.T))
                              if len(activeint.T) > 1: 
                                    
                                    #create interaction frame list
                                    interaction_frame = []
                                    for a in range(len(activeint.T)):
                                          b= list(range(activeint[a][0], activeint[a][1], 1))
                                          interaction_frame.append(b)
                                    
                                    #convert 2d list to 1d list
                                    interaction_frame_1d = list(chain(*interaction_frame))
                                    
                                    #convert interaction peroid to binary mask
                                    record_list = np.arange(6000)
                                    new_record_list = [ele for ele in record_list if ele not in BF]
                                    
                                    interaction_mask = []
                                    for j in new_record_list:
                                          if j in interaction_frame_1d:
                                                interaction_mask.append(1)
                                          else:
                                                interaction_mask.append(0)
                                    
                                    #create cosine similarity of each cells
                                    if not os.path.exists(fdir + '/cell_csv'): os.mkdir(fdir + '/cell_csv')
                                    for r in tqdm(np.arange(len(accnpt))):
                                          row = accnpt[r,:]
                                          similarity = cosine_similarity([row], [interaction_mask])
                                          #export aligned clacium traces of each cell into CSV files    
                                          pd.DataFrame(similarity).to_csv(fdir + '/cell_csv/cell' + str(r+1) + '.csv')
                                          
                                          
                                    #############
                                    #generate cosine similarity of shuffled interaction period and calcium signal
                                    np.random.seed(17)
                                    csall = []
                                    for a in tqdm(range (1000)): # shuffled 1000 times
                                          cs_ind = []
                                          for s in np.arange(len(accnpt)):
                                                rowr = accnpt[s,:]
                                                rand_mask = random_interaction_interaction_mask(interaction_frame, interaction_frame_1d, 6000, BF)
                                                cs = cosine_similarity([rowr], [rand_mask]).tolist()
                                                cs2 = list(chain(*cs))
                                                for i in cs2:
                                                      cs_ind.append(i)
                                          csall.append(cs_ind)
                                    #save as CSV files
                                    pd.DataFrame(np.transpose(csall)).to_csv(fdir + '/shuffled_cosine_similarity.csv')
                                    
                                    
                                    ##############
                                    #save calcium signal of interaction onset and offset
                                    if not os.path.exists(fdir + '/cell_onset'): os.mkdir(fdir + '/cell_onset')
                                    if not os.path.exists(fdir + '/cell_offset'): os.mkdir(fdir + '/cell_offset')
                                    for r in tqdm(np.arange(len(accnpt_plus))):
                                          row = accnpt_plus[r,:]
                                          
                                          cell_onset_traces = []
                                          cell_offset_traces = []
                                          
                                          for i in range(len(activeint.T)):
                                                x = int(activeint[0:1][i])
                                                y = int(activeint[1:2][i])
                                                
                                                onset = row[ (x+51) -50 : (x+51) +50]
                                                offset = row[ (y+51) -50 : (y+51) +50]
                                                cell_onset_traces.append(onset)
                                                cell_offset_traces.append(offset)
                                                
                                          #export aligned clacium traces of each cell into CSV files    
                                          pd.DataFrame(cell_onset_traces).to_csv(fdir + '/cell_onset/cell' + str(r+1) + '.csv')
                                          pd.DataFrame(cell_offset_traces).to_csv(fdir + '/cell_offset/cell' + str(r+1) + '.csv')
                                          
                                          
                                    ##############
                                    #calculate mean of cell activity during interaction
                                    behavior_activity_mean = []
                                    for r in tqdm(np.arange(len(accnpt_plus))):
                                          row = accnpt_plus[r,:]
                                          interaction_activity = []
                                          for s in interaction_frame_1d:
                                                interaction_activity.append(row[ (s+51) : (s+51) +1])
                                          behavior_activity_mean.append(np.nanmean(interaction_activity))
                                    
                                    
                                    if not os.path.exists(fdir + '/save_figures'): os.mkdir(fdir + '/save_figures')
                                    
                                    
                                    ##############
                                    #create auROC of each cells
                                    if not os.path.exists(fdir + '/cell_auROC_csv'): os.mkdir(fdir + '/cell_auROC_csv')
                                    for r in tqdm(np.arange(len(accnpt))):
                                          row = accnpt[r,:]
                                          auROC(row, interaction_mask, r, fdir)

                                    ##############
                                    #hypothesis testing, plot and save figures (Activation)
                                    if not os.path.exists(fdir + '/save_figures/Activation'): os.mkdir(fdir + '/save_figures/Activation')
                                    shuffled_cs = pd.read_csv(fdir + '/shuffled_cosine_similarity.csv', index_col= 0, header= 0 ).values

                                    cs_P_all = []
                                    cs_all = []
                                    pearson_r = [] #r_value of pearson
                                    pearson_p = [] #p_value of pearson
                                    
                                    
                                    for roci in range (len(accnpt)):
                                          #read the auROC csv
                                          TPR = pd.read_csv(fdir + '/cell_auROC_csv/cell' + str(roci+1) + '_TPRr.csv', index_col= 0, header= 0).values
                                          FPR = pd.read_csv(fdir + '/cell_auROC_csv/cell' + str(roci+1) + '_FPRr.csv', index_col= 0, header= 0).values
                                          activity = pd.read_csv(fdir + '/cell_auROC_csv/cell' + str(roci+1) + '_mask1_activity.csv', index_col= 0, header= 0).values
                                          baseline = pd.read_csv(fdir + '/cell_auROC_csv/cell' + str(roci+1) + '_mask0_activity.csv', index_col= 0, header= 0).values
                                          
                                          
                                          #hypothesis testing by using sum as test statistic
                                          cs = pd.read_csv(fdir + '/cell_csv/cell' + str(roci+1) + '.csv', index_col= 0, header= 0).values
                                          cs_p =  np.sum(shuffled_cs [roci]>= cs)/len(shuffled_cs[roci])
                                          cs_P_all.append (cs_p)
                                          cs_all.append (cs[0][0])
                                          
                                          onset_traces = pd.read_csv(fdir + '/cell_onset/cell' + str(roci+1) + '.csv', index_col = 0, header =0 )
                                          offset_traces = pd.read_csv(fdir + '/cell_offset/cell' + str(roci+1) + '.csv', index_col = 0, header =0 )

                                          ###################
                                          #calculate correlation between activity and interaction experience
                                          accnpt_z= accnpt
                                          if (blankfile != ""): accnpt_z= accnpt2

                                          row2 = accnpt_z[roci,:]
                                          
                                          signal_mean = []
                                          for i in interaction_frame:
                                                signal= []
                                                for j in i:
                                                      signal.append(row2[j:j+1])
                                                mean = np.nanmean(signal)
                                                # print(mean)
                                                signal_mean.append(mean)
                                          # print("len", len(signal_mean))
                                          
                                          # signal_mean_array = np.array(signal_mean).reshape(-1,1)
                                          interaction_experience = [i+1 for i in range(len(signal_mean))]
                                          
                                          PC_Results = pearsonr(signal_mean, interaction_experience)
                                          PC_corr = [np.round(c, 2) for c in PC_Results]
                                          pearson_r.append(PC_corr[0])
                                          pearson_p.append(PC_corr[1])
                                          ###################

                                          #set significance level (default: 0.05)   
                                          if cs_p < 0.05 :  #change the comparison operator for ploting interaction-responsive or non interaction-response cells
                                                onset_traces_average = np.mean(onset_traces)
                                                onset_traces_sem = np.std (onset_traces)/ np.sqrt (len (onset_traces))
                                                onset_x = np.arange (len(onset_traces_average))
                                                
                                                offset_traces_average = np.mean(offset_traces)
                                                offset_traces_sem = np.std (offset_traces)/ np.sqrt (len (offset_traces))
                                                offset_x = np.arange (len(offset_traces_average))
                                                
                                                print('Activation')

                                                #set figure size
                                                plt.figure (figsize = (4,5))
                                          
                                                #plot z score transformed entire calcium traces
                                                fig = plt.subplot (4,1,1)
                                                
                                                for xc in interaction_frame_1d:
                                                      plt.axvline(x= xc, color= "#e5e5e5", linestyle= '-', linewidth= 1, alpha= 0.3)
                                                acc_z= accrawt_z
                                                if (blankfile != ""): acc_z= accrawt_z2
                                                plt.plot (np.transpose(acc_z[roci]), color = "red", alpha = 1, linewidth = 0.5)
                                                new_xtickstime = np.linspace(0, 6000, 11)
                                                plt.xticks(new_xtickstime,[0,1,2,3,4,5,6,7,8,9,10], rotation=0)
                                                plt.xlabel ("Time (min)")
                                                plt.ylabel ("s.d.")
                                                plt.title("cell"+str(roci+1))

                                                #plot random shuffled_sum distribution and true distribution
                                                fig = plt.subplot (4,3,4)
                                                plt.hist(shuffled_cs[roci], density= False, bins= int(np.sqrt(len(np.transpose(shuffled_cs)))), color= "#FFCCCC", edgecolor='none')
                                                plt.axvline(x= cs, color='red', linewidth= 1)
                                                plt.axvline(x= np.percentile(shuffled_cs[roci], 95), color='gray', linewidth= 1, linestyle= '--')
                                                plt.axvline(x= np.percentile(shuffled_cs[roci], 5), color='gray', linewidth= 1, linestyle= '--')
                                                plt.xlabel ("Cosine_similarity")
                                                plt.ylabel ("Event number")
                                                plt.title("Similarity_p = "+ str(cs_p), loc= "right")
                                                
                                                #plot heat map of calcium traces    
                                                fig = plt.subplot (4,3,5)
                                                heat= plt.imshow(onset_traces, aspect = "auto", cmap="jet", interpolation = "nearest")
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                # y_ticks = np.arange (len(onset_traces.values), step = 2)
                                                # plt.yticks(y_ticks, y_ticks +1, rotation = 0)
                                                ax = plt.gca()
                                                ax.spines['left'].set_color('none')
                                                plt.yticks([])
                                                plt.xlabel ("Time (s)")
                                                plt.ylabel ("Interaction bouts" + " (" + str(len(onset_traces.values))+ ")")
                                                plt.axvline(x=50, color='w', linestyle=':', linewidth=1)
                                                plt.title ("interaction onset")
                                                cbar_ticks = np.linspace (int(min(accnpt[roci])), int(max (accnpt[roci])),2)
                                                cb = plt.colorbar (heat, shrink = 0.8, drawedges = False, spacing = "proportional")
                                                cb.ax.set_title ("s.d.", loc = "left")
                                                cb.outline.set_visible (False)
                                                cb.ax.tick_params (axis = "y", direction = "out", pad = 0.2, length = 1)
                                                
                                                #plot mean + - sem traces of interaction onset
                                                fig = plt.subplot (4,3,6)
                                                plt.plot (onset_x, onset_traces_average, color = "red")
                                                plt.fill_between (onset_x, onset_traces_average - onset_traces_sem, onset_traces_average + onset_traces_sem, 
                                                                  alpha = 0.2, edgecolor = "none", color = "red")    
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                plt.xlabel ("Time (s)")
                                                plt.axvline(x=50, color='black', linestyle=':', linewidth=1)
                                                plt.title ("interaction onset")
                                                plt.ylabel ("s.d.")
                                                
                                                #plot correlation between interaction mean activity and interaction experience
                                                fig = plt.subplot (4,3,7)
                                                regression = pd.DataFrame({"Average_activity":signal_mean, "Interaction_experience":interaction_experience})
                                                sns.regplot(x= "Interaction_experience", y= "Average_activity", data= regression, color= "black", 
                                                            line_kws= {"color":"red", "alpha":0.7, "lw":2}, scatter_kws={'s':2})
                                                plt.title("pearsonr "+" r=%s, p=%s" % (PC_corr[0], PC_corr[1]))
                                                
                                                #plot heat map of calcium traces    
                                                fig = plt.subplot (4,3,8)
                                                heat= plt.imshow(offset_traces, aspect = "auto", cmap="jet", interpolation = "nearest")
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                # y_ticks = np.arange (len(onset_traces.values), step = 2)
                                                # plt.yticks(y_ticks, y_ticks +1, rotation = 0)
                                                ax = plt.gca()
                                                ax.spines['left'].set_color('none')
                                                plt.yticks([])
                                                plt.xlabel ("Time (s)")
                                                plt.ylabel ("Interaction bouts" + " (" + str(len(offset_traces.values))+ ")")
                                                plt.axvline(x=50, color='w', linestyle=':', linewidth=1)
                                                plt.title ("interaction offset")
                                                cbar_ticks = np.linspace (int(min(accnpt[roci])), int(max (accnpt[roci])),2)
                                                cb = plt.colorbar (heat, shrink = 0.8, drawedges = False, spacing = "proportional")
                                                cb.ax.set_title ("s.d.", loc = "left")
                                                cb.outline.set_visible (False)
                                                cb.ax.tick_params (axis = "y", direction = "out", pad = 0.2, length = 1)
                                                
                                                #plot mean + - sem traces of interaction offset
                                                fig = plt.subplot (4,3,9)
                                                plt.plot (offset_x, offset_traces_average, color = "red")
                                                plt.fill_between (offset_x, offset_traces_average - offset_traces_sem, offset_traces_average + offset_traces_sem, 
                                                                  alpha = 0.2, edgecolor = "none", color = "red")    
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                plt.xlabel ("Time (s)")
                                                plt.axvline(x=50, color='black', linestyle=':', linewidth=1)
                                                plt.title ("interaction offset")
                                                plt.ylabel ("s.d.")
                                                
                                                #plot calcium signal 0 to 3 s value distribution and baseline signal value distribution
                                                fig = plt.subplot (4,3,10)
                                                plt.hist (baseline, density = True, bins = 100, histtype = "stepfilled", edgecolor = "none", alpha=0.5)
                                                plt.hist (activity, density = True, bins = 100, histtype = "stepfilled", edgecolor = "none", alpha=0.5, color = "#fcb001")
                                                plt.xlabel ("s.d.")
                                                plt.ylabel ("Probability")
                                                plt.legend(["baseline","signal"], loc='upper right', handlelength=0.5, labelspacing=0.5, borderpad=0.5, framealpha = 0)

                                                #plot Receiver Operating Characteristic curves
                                                fig = plt.subplot (4,3,11)
                                                plt.plot([0, 1], [0, 1], 'k--', linewidth = 0.5)
                                                plt.plot(FPR, TPR, color = "#fcb001", linewidth = 1)
                                                plt.fill_between (list(chain(*FPR)), list(chain(*TPR)), 0, alpha=0.2, color = "#fcb001")
                                                plt.xlabel('P (baseline > threshold)')
                                                plt.ylabel('P (signal > threshold)')
                                                plt.title('RS = (auROC -0.5)*2')
                                                
                                                
                                                plt.tight_layout()
                                                sns.despine(offset = 2)
                                          
                                                #save figures : change the target folder to store figures with different properties
                                                plt.savefig(fdir + '/save_figures/Activation/cell' + str(roci+1) + '.pdf', transparent= True, bbox_inches="tight")
                                                plt.savefig(fdir + '/save_figures/Activation/cell' + str(roci+1) + '.tif', transparent= True, bbox_inches="tight", dpi= 720)
                                                plt.cla()
                                                #plt.show()

                                    ##############
                                    #hypothesis testing, plot and save figures (Inhibition)
                                    if not os.path.exists(fdir + '/save_figures/Inhibition'): os.mkdir(fdir + '/save_figures/Inhibition')
                                    shuffled_cs = pd.read_csv(fdir + '/shuffled_cosine_similarity.csv', index_col= 0, header= 0 ).values

                                    cs_P_all = []
                                    cs_all = []
                                    pearson_r = [] #r_value of pearson
                                    pearson_p = [] #p_value of pearson
                                    
                                    
                                    for roci in range (len(accnpt)):
                                          #hypothesis testing by using sum as test statistic
                                          cs = pd.read_csv(fdir + '/cell_csv/cell' + str(roci+1) + '.csv', index_col= 0, header= 0).values
                                          cs_p =  np.sum(shuffled_cs [roci]>= cs)/len(shuffled_cs[roci])
                                          cs_P_all.append (cs_p)
                                          cs_all.append (cs[0][0])
                                          
                                          onset_traces = pd.read_csv(fdir + '/cell_onset/cell' + str(roci+1) + '.csv', index_col = 0, header =0 )
                                          offset_traces = pd.read_csv(fdir + '/cell_offset/cell' + str(roci+1) + '.csv', index_col = 0, header =0 )
                              
                                          ###################
                                          #calculate correlation between activity and interaction experience
                                          accnpt_z= accnpt
                                          if (blankfile != ""): accnpt_z= accnpt2

                                          row2 = accnpt_z[roci,:]
                                          
                                          signal_mean = []
                                          for i in interaction_frame:
                                                signal= []
                                                for j in i:
                                                      signal.append(row2[j:j+1])
                                                mean = np.nanmean(signal)
                                                # print(mean)
                                                signal_mean.append(mean)
                                          # print("len", len(signal_mean))
                                          
                                          # signal_mean_array = np.array(signal_mean).reshape(-1,1)
                                          interaction_experience = [i+1 for i in range(len(signal_mean))]
                                          
                                          PC_Results = pearsonr(signal_mean, interaction_experience)
                                          PC_corr = [np.round(c, 2) for c in PC_Results]
                                          pearson_r.append(PC_corr[0])
                                          pearson_p.append(PC_corr[1])
                                          ###################

                                          #set significance level (default: 0.05)   
                                          if cs_p > 0.95 :  #change the comparison operator for ploting interaction-responsive or non interaction-response cells
                                                onset_traces_average = np.mean(onset_traces)
                                                onset_traces_sem = np.std (onset_traces)/ np.sqrt (len (onset_traces))
                                                onset_x = np.arange (len(onset_traces_average))
                                                
                                                offset_traces_average = np.mean(offset_traces)
                                                offset_traces_sem = np.std (offset_traces)/ np.sqrt (len (offset_traces))
                                                offset_x = np.arange (len(offset_traces_average))
                                                
                                                print('Inhibition')

                                                #set figure size
                                                plt.figure (figsize = (4,5))
                                          
                                                #plot z score transformed entire calcium traces
                                                fig = plt.subplot (4,1,1)
                                                
                                                for xc in interaction_frame_1d:
                                                      plt.axvline(x= xc, color= "#e5e5e5", linestyle= '-', linewidth= 1, alpha= 0.3)
                                                      
                                                acc_z= accrawt_z
                                                if (blankfile != ""): acc_z= accrawt_z2
                                                plt.plot (np.transpose(acc_z[roci]), color = "Blue", alpha = 1, linewidth = 0.5)
                                                new_xtickstime = np.linspace(0, 6000, 11)
                                                plt.xticks(new_xtickstime,[0,1,2,3,4,5,6,7,8,9,10], rotation=0)
                                                plt.xlabel ("Time (min)")
                                                plt.ylabel ("s.d.")
                                                plt.title("cell"+str(roci+1))

                                                #plot random shuffled_sum distribution and true distribution
                                                fig = plt.subplot (4,3,4)
                                                plt.hist(shuffled_cs[roci], density= False, bins= int(np.sqrt(len(np.transpose(shuffled_cs)))), color= "#b2b2ff", edgecolor='none')
                                                plt.axvline(x= cs, color='blue', linewidth= 1)
                                                plt.axvline(x= np.percentile(shuffled_cs[roci], 95), color='gray', linewidth= 1, linestyle= '--')
                                                plt.axvline(x= np.percentile(shuffled_cs[roci], 5), color='gray', linewidth= 1, linestyle= '--')
                                                plt.xlabel ("Cosine_similarity")
                                                plt.ylabel ("Event number")
                                                plt.title("Similarity_p = "+ str(cs_p), loc= "right")
                                                
                                                #plot heat map of calcium traces    
                                                fig = plt.subplot (4,3,5)
                                                heat= plt.imshow(onset_traces, aspect = "auto", cmap="jet", interpolation = "nearest")
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                # y_ticks = np.arange (len(onset_traces.values), step = 2)
                                                # plt.yticks(y_ticks, y_ticks +1, rotation = 0)
                                                ax = plt.gca()
                                                ax.spines['left'].set_color('none')
                                                plt.yticks([])
                                                plt.xlabel ("Time (s)")
                                                plt.ylabel ("Interaction bouts" + " (" + str(len(onset_traces.values))+ ")")
                                                plt.axvline(x=50, color='w', linestyle=':', linewidth=1)
                                                plt.title ("interaction onset")
                                                cbar_ticks = np.linspace (int(min(accnpt[roci])), int(max (accnpt[roci])),2)
                                                cb = plt.colorbar (heat, shrink = 0.8, drawedges = False, spacing = "proportional")
                                                cb.ax.set_title ("s.d.", loc = "left")
                                                cb.outline.set_visible (False)
                                                cb.ax.tick_params (axis = "y", direction = "out", pad = 0.2, length = 1)
                                                
                                                #plot mean + - sem traces of interaction onset
                                                fig = plt.subplot (4,3,6)
                                                plt.plot (onset_x, onset_traces_average, color = "blue")
                                                plt.fill_between (onset_x, onset_traces_average - onset_traces_sem, onset_traces_average + onset_traces_sem, 
                                                                  alpha = 0.2, edgecolor = "none", color = "blue")    
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                plt.xlabel ("Time (s)")
                                                plt.axvline(x=50, color='black', linestyle=':', linewidth=1)
                                                plt.title ("interaction onset")
                                                plt.ylabel ("s.d.")
                                                
                                                #plot correlation between interaction mean activity and interaction experience
                                                fig = plt.subplot (4,3,7)
                                                regression = pd.DataFrame({"Average_activity":signal_mean, "Interaction_experience":interaction_experience})
                                                sns.regplot(x= "Interaction_experience", y= "Average_activity", data= regression, color= "black", 
                                                            line_kws= {"color":"red", "alpha":0.7, "lw":2}, scatter_kws={'s':2})
                                                plt.title("pearsonr "+" r=%s, p=%s" % (PC_corr[0], PC_corr[1]))
                                                
                                                #plot heat map of calcium traces    
                                                fig = plt.subplot (4,3,8)
                                                heat= plt.imshow(offset_traces, aspect = "auto", cmap="jet", interpolation = "nearest")
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                # y_ticks = np.arange (len(onset_traces.values), step = 2)
                                                # plt.yticks(y_ticks, y_ticks +1, rotation = 0)
                                                ax = plt.gca()
                                                ax.spines['left'].set_color('none')
                                                plt.yticks([])
                                                plt.xlabel ("Time (s)")
                                                plt.ylabel ("Interaction bouts" + " (" + str(len(offset_traces.values))+ ")")
                                                plt.axvline(x=50, color='w', linestyle=':', linewidth=1)
                                                plt.title ("interaction offset")
                                                cbar_ticks = np.linspace (int(min(accnpt[roci])), int(max (accnpt[roci])),2)
                                                cb = plt.colorbar (heat, shrink = 0.8, drawedges = False, spacing = "proportional")
                                                cb.ax.set_title ("s.d.", loc = "left")
                                                cb.outline.set_visible (False)
                                                cb.ax.tick_params (axis = "y", direction = "out", pad = 0.2, length = 1)
                                                
                                                #plot mean + - sem traces of interaction offset
                                                fig = plt.subplot (4,3,9)
                                                plt.plot (offset_x, offset_traces_average, color = "blue")
                                                plt.fill_between (offset_x, offset_traces_average - offset_traces_sem, offset_traces_average + offset_traces_sem, 
                                                                  alpha = 0.2, edgecolor = "none", color = "blue")    
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                plt.xlabel ("Time (s)")
                                                plt.axvline(x=50, color='black', linestyle=':', linewidth=1)
                                                plt.title ("interaction offset")
                                                plt.ylabel ("s.d.")  
                                                
                                                #plot calcium signal 0 to 3 s value distribution and baseline signal value distribution
                                                fig = plt.subplot (4,3,10)
                                                plt.hist (baseline, density = True, bins = 100, histtype = "stepfilled", edgecolor = "none", alpha=0.5)
                                                plt.hist (activity, density = True, bins = 100, histtype = "stepfilled", edgecolor = "none", alpha=0.5, color = "#fcb001")
                                                plt.xlabel ("s.d.")
                                                plt.ylabel ("Probability")
                                                plt.legend(["baseline","signal"], loc='upper right', handlelength=0.5, labelspacing=0.5, borderpad=0.5, framealpha = 0)

                                                #plot Receiver Operating Characteristic curves
                                                fig = plt.subplot (4,3,11)
                                                plt.plot([0, 1], [0, 1], 'k--', linewidth = 0.5)
                                                plt.plot(FPR, TPR, color = "#fcb001", linewidth = 1)
                                                plt.fill_between (list(chain(*FPR)), list(chain(*TPR)), 0, alpha=0.2, color = "#fcb001")
                                                plt.xlabel('P (baseline > threshold)')
                                                plt.ylabel('P (signal > threshold)')
                                                plt.title('RS = (auROC -0.5)*2')
                                                
                                                plt.tight_layout()
                                                sns.despine(offset = 2)

                                                #save figures : change the target folder to store figures with different properties
                                                plt.savefig(fdir + '/save_figures/Inhibition/cell' + str(roci+1) + '.pdf', transparent= True, bbox_inches="tight")
                                                plt.savefig(fdir + '/save_figures/Inhibition/cell' + str(roci+1) + '.tif', transparent= True, bbox_inches="tight", dpi= 720)
                                                plt.cla()
                                                #plt.show()

                                    ##############
                                    #hypothesis testing, plot and save figures (Irrelevant)
                                    if not os.path.exists(fdir + '/save_figures/Irrelevant'): os.mkdir(fdir + '/save_figures/Irrelevant')
                                    shuffled_cs = pd.read_csv(fdir + '/shuffled_cosine_similarity.csv', index_col= 0, header= 0 ).values

                                    cs_P_all = []
                                    cs_all = []
                                    pearson_r = [] #r_value of pearson
                                    pearson_p = [] #p_value of pearson
                                    
                                    
                                    for roci in range (len(accnpt)):
                                          #hypothesis testing by using sum as test statistic
                                          cs = pd.read_csv(fdir + '/cell_csv/cell' + str(roci+1) + '.csv', index_col= 0, header= 0).values
                                          cs_p =  np.sum(shuffled_cs [roci]>= cs)/len(shuffled_cs[roci])
                                          cs_P_all.append (cs_p)
                                          cs_all.append (cs[0][0])
                                          
                                          onset_traces = pd.read_csv(fdir + '/cell_onset/cell' + str(roci+1) + '.csv', index_col = 0, header =0 )
                                          offset_traces = pd.read_csv(fdir + '/cell_offset/cell' + str(roci+1) + '.csv', index_col = 0, header =0 )
                              
                                          ###################
                                          #calculate correlation between activity and interaction experience
                                          accnpt_z= accnpt
                                          if (blankfile != ""): accnpt_z= accnpt2

                                          row2 = accnpt_z[roci,:]
                                          
                                          signal_mean = []
                                          for i in interaction_frame:
                                                signal= []
                                                for j in i:
                                                      signal.append(row2[j:j+1])
                                                mean = np.nanmean(signal)
                                                # print(mean)
                                                signal_mean.append(mean)
                                          # print("len", len(signal_mean))
                                          
                                          # signal_mean_array = np.array(signal_mean).reshape(-1,1)
                                          interaction_experience = [i+1 for i in range(len(signal_mean))]
                                          
                                          PC_Results = pearsonr(signal_mean, interaction_experience)
                                          PC_corr = [np.round(c, 2) for c in PC_Results]
                                          pearson_r.append(PC_corr[0])
                                          pearson_p.append(PC_corr[1])
                                          ###################
                              

                                          #set significance level (default: 0.05)   
                                          if 0.05< cs_p < 0.95 :  #change the comparison operator for ploting interaction-responsive or non interaction-response cells
                                                onset_traces_average = np.mean(onset_traces)
                                                onset_traces_sem = np.std (onset_traces)/ np.sqrt (len (onset_traces))
                                                onset_x = np.arange (len(onset_traces_average))
                                                
                                                offset_traces_average = np.mean(offset_traces)
                                                offset_traces_sem = np.std (offset_traces)/ np.sqrt (len (offset_traces))
                                                offset_x = np.arange (len(offset_traces_average))

                                                print('Irrelevant')

                                                #set figure size
                                                plt.figure (figsize = (4,5))
                                          
                                                #plot z score transformed entire calcium traces
                                                fig = plt.subplot (4,1,1)
                                                
                                                for xc in interaction_frame_1d:
                                                      plt.axvline(x= xc, color= "#e5e5e5", linestyle= '-', linewidth= 1, alpha= 0.3)
                                                      
                                                acc_z= accrawt_z
                                                if (blankfile != ""): acc_z= accrawt_z2
                                                plt.plot (np.transpose(acc_z[roci]), color = "black", alpha = 1, linewidth = 0.5)
                                                new_xtickstime = np.linspace(0, 6000, 11)
                                                plt.xticks(new_xtickstime,[0,1,2,3,4,5,6,7,8,9,10], rotation=0)
                                                plt.xlabel ("Time (min)")
                                                plt.ylabel ("s.d.")
                                                plt.title("cell"+str(roci+1))

                                                #plot random shuffled_sum distribution and true distribution
                                                fig = plt.subplot (4,3,4)
                                                plt.hist(shuffled_cs[roci], density= False, bins= int(np.sqrt(len(np.transpose(shuffled_cs)))), color= "#bfbfbf", edgecolor='none')
                                                plt.axvline(x= cs, color='black', linewidth= 1)
                                                plt.axvline(x= np.percentile(shuffled_cs[roci], 95), color='gray', linewidth= 1, linestyle= '--')
                                                plt.axvline(x= np.percentile(shuffled_cs[roci], 5), color='gray', linewidth= 1, linestyle= '--')
                                                plt.xlabel ("Cosine_similarity")
                                                plt.ylabel ("Event number")
                                                plt.title("Similarity_p = "+ str(cs_p), loc= "right")
                                                
                                                #plot heat map of calcium traces    
                                                fig = plt.subplot (4,3,5)
                                                heat= plt.imshow(onset_traces, aspect = "auto", cmap="jet", interpolation = "nearest")
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                # y_ticks = np.arange (len(onset_traces.values), step = 2)
                                                # plt.yticks(y_ticks, y_ticks +1, rotation = 0)
                                                ax = plt.gca()
                                                ax.spines['left'].set_color('none')
                                                plt.yticks([])
                                                plt.xlabel ("Time (s)")
                                                plt.ylabel ("Interaction bouts" + " (" + str(len(onset_traces.values))+ ")")
                                                plt.axvline(x=50, color='w', linestyle=':', linewidth=1)
                                                plt.title ("interaction onset")
                                                cbar_ticks = np.linspace (int(min(accnpt[roci])), int(max (accnpt[roci])),2)
                                                cb = plt.colorbar (heat, shrink = 0.8, drawedges = False, spacing = "proportional")
                                                cb.ax.set_title ("s.d.", loc = "left")
                                                cb.outline.set_visible (False)
                                                cb.ax.tick_params (axis = "y", direction = "out", pad = 0.2, length = 1)
                                                
                                                #plot mean + - sem traces of interaction onset
                                                fig = plt.subplot (4,3,6)
                                                plt.plot (onset_x, onset_traces_average, color = "black")
                                                plt.fill_between (onset_x, onset_traces_average - onset_traces_sem, onset_traces_average + onset_traces_sem, 
                                                                  alpha = 0.2, edgecolor = "none", color = "gray")    
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                plt.xlabel ("Time (s)")
                                                plt.axvline(x=50, color='black', linestyle=':', linewidth=1)
                                                plt.title ("interaction onset")
                                                plt.ylabel ("s.d.")
                                                
                                                #plot correlation between interaction mean activity and interaction experience
                                                fig = plt.subplot (4,3,7)
                                                regression = pd.DataFrame({"Average_activity":signal_mean, "Interaction_experience":interaction_experience})
                                                sns.regplot(x= "Interaction_experience", y= "Average_activity", data= regression, color= "black", 
                                                            line_kws= {"color":"red", "alpha":0.7, "lw":2}, scatter_kws={'s':2})
                                                plt.title("pearsonr "+" r=%s, p=%s" % (PC_corr[0], PC_corr[1]))
                                                
                                                #plot heat map of calcium traces    
                                                fig = plt.subplot (4,3,8)
                                                heat= plt.imshow(offset_traces, aspect = "auto", cmap="jet", interpolation = "nearest")
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                # y_ticks = np.arange (len(onset_traces.values), step = 2)
                                                # plt.yticks(y_ticks, y_ticks +1, rotation = 0)
                                                ax = plt.gca()
                                                ax.spines['left'].set_color('none')
                                                plt.yticks([])
                                                plt.xlabel ("Time (s)")
                                                plt.ylabel ("Interaction bouts" + " (" + str(len(offset_traces.values))+ ")")
                                                plt.axvline(x=50, color='w', linestyle=':', linewidth=1)
                                                plt.title ("interaction offset")
                                                cbar_ticks = np.linspace (int(min(accnpt[roci])), int(max (accnpt[roci])),2)
                                                cb = plt.colorbar (heat, shrink = 0.8, drawedges = False, spacing = "proportional")
                                                cb.ax.set_title ("s.d.", loc = "left")
                                                cb.outline.set_visible (False)
                                                cb.ax.tick_params (axis = "y", direction = "out", pad = 0.2, length = 1)
                                                
                                                #plot mean + - sem traces of interaction offset
                                                fig = plt.subplot (4,3,9)
                                                plt.plot (offset_x, offset_traces_average, color = "black")
                                                plt.fill_between (offset_x, offset_traces_average - offset_traces_sem, offset_traces_average + offset_traces_sem, 
                                                                  alpha = 0.2, edgecolor = "none", color = "gray")    
                                                new_xticks = np.linspace(0, 100, 11)
                                                plt.xticks(new_xticks,[-5,-4,-3,-2,-1,0,1,2,3,4,5], rotation=0)
                                                plt.xlabel ("Time (s)")
                                                plt.axvline(x=50, color='black', linestyle=':', linewidth=1)
                                                plt.title ("interaction offset")
                                                plt.ylabel ("s.d.")  
                                                
                                                #plot calcium signal 0 to 3 s value distribution and baseline signal value distribution
                                                fig = plt.subplot (4,3,10)
                                                plt.hist (baseline, density = True, bins = 100, histtype = "stepfilled", edgecolor = "none", alpha=0.5)
                                                plt.hist (activity, density = True, bins = 100, histtype = "stepfilled", edgecolor = "none", alpha=0.5, color = "#fcb001")
                                                plt.xlabel ("s.d.")
                                                plt.ylabel ("Probability")
                                                plt.legend(["baseline","signal"], loc='upper right', handlelength=0.5, labelspacing=0.5, borderpad=0.5, framealpha = 0)

                                                #plot Receiver Operating Characteristic curves
                                                fig = plt.subplot (4,3,11)
                                                plt.plot([0, 1], [0, 1], 'k--', linewidth = 0.5)
                                                plt.plot(FPR, TPR, color = "#fcb001", linewidth = 1)
                                                plt.fill_between (list(chain(*FPR)), list(chain(*TPR)), 0, alpha=0.2, color = "#fcb001")
                                                plt.xlabel('P (baseline > threshold)')
                                                plt.ylabel('P (signal > threshold)')
                                                plt.title('RS = (auROC -0.5)*2')
                                                
                                                plt.tight_layout()
                                                sns.despine(offset = 2)

                                                #save figures : change the target folder to store figures with different properties
                                                plt.savefig(fdir + '/save_figures/Irrelevant/cell' + str(roci+1) + '.pdf', transparent= True, bbox_inches="tight")
                                                plt.savefig(fdir + '/save_figures/Irrelevant/cell' + str(roci+1) + '.tif', transparent= True, bbox_inches="tight", dpi= 720)
                                                plt.cla()
                                                #plt.show()
                                          #########################################################################

                                    #Generate statistic summary table and save to csv file
                                    cell_ID = []
                                    for idx in range (len(accnpt)):
                                          cell_id = idx +1
                                          cell_ID.append (cell_id)
                              
                                    behavior_responsive = []
                                    for res in range (len (accnpt)):
                                          if cs_P_all [res] < 0.05:
                                                responsive = str("activation")
                                                behavior_responsive.append (responsive)
                                          elif cs_P_all [res] > 0.95:
                                                responsive = str("inhibition")
                                                behavior_responsive.append (responsive)
                                          else:
                                                responsive = str("")
                                                behavior_responsive.append (responsive)

                                    Regression_index = []
                                    for rp in range (len(accnpt)):
                                          if pearson_r [rp] > 0 and pearson_p [rp] < 0.05 and behavior_responsive[rp] != "":
                                                LR = str("positive")
                                                Regression_index.append(LR)
                                          elif pearson_r [rp] < 0 and pearson_p [rp] < 0.05 and behavior_responsive[rp] != "":
                                                LR = str("negative")
                                                Regression_index.append(LR)
                                          else:
                                                LR = str("")
                                                Regression_index.append(LR)

                                    df = pd.DataFrame (list (zip(cell_ID, cs_all, cs_P_all, behavior_responsive, pearson_r, pearson_p, Regression_index, behavior_activity_mean)), 
                                                      columns = ["cell ID", "similarity", "CS P *", "behavior correlation *", "pearson_r", "pearson_p", "Experience_dependent_cell_type", "mean_activty_during_behavior"])
                                    #save to csv file
                                    df.to_csv(fdir + '/summary_P.csv', index= False)

                              else:
                                    cell_ID = []
                                    for idx in range (len(accnpt)):
                                          cell_id = idx +1
                                          cell_ID.append (cell_id)
                                    df = pd.DataFrame(np.nan, index=cell_ID, columns = ["similarity", "CS P *", "behavior correlation *", "pearson_r", "pearson_p", "Experience_dependent_cell_type", "mean_activty_during_behavior"])
                                    df.index.names= ["cell ID"]
                                    #save to csv file
                                    df.to_csv(fdir + '/summary_P.csv')

                  ########################################
                  #For object and social cells
                  dataprocess (oefile)
                  dataprocess (rsifile)