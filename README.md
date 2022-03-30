# PROB_R
R codes for PROB algorithm
main.R  --- main function, containing codes in the "step-by-step details" part of the protocol

All codes below contain executable functions required in main.R, please execute them in your R environment before analysis.
PROB_GEOinstall.R --- function to download the GSE7390 dataset. The downloaded data will be saved as "GSE7390.csv" in your current working directory.
Progression_Inference.R --- function to infer the disease progression. 
                            Input objects:
                                 Your cross-sectional transcriptomic data, genes in row and samples in column. The last row is grade information of each sample
                            Output objects:
                                 Accumulated_Transition_Matrix --- the accumulated transition matrix used in the step of feature extraction with diffusion maps    
                                 Temporal_Progression --- PPD value of each gene
                                 Ordered_Data -- reordered samples against temporal progression
                                 Sampled_Time -- sampled time marks in trajectory smoothing
                                 Order --- The order of each sample in the progression
ODE_Bayesian_Lasso.R -- function to infer the GRN with ODE Bayesian Lasso method. 
                        Input objects:
                            Data_ordered --- reordered samples against temporal progression of TCGs
                            Time_Sampled ---  sampled time marks in trajectory smoothing
                            Prior_Network --- prior knowledge to build GRNs                          
                        Output objects:
                            Ajacent_Matrix --- the matrix of the posterior mean of each parameter
                            Adjusted_Ajacent_Matrix --- 0-1 matrix, 1 entries occur only when edges have 95% or more credible level
                            Presence_Probability --- credible level of each edge
                            Standard_Deviations --- posterior standard deviation of each parameter
trans_cytoscape.R, BL_to_csv.R --- functions to save the network into a .csv file
Locate_Key_Genes.R --- function to cauculate hub scores of each TCG
                         Input objects:
                            BLin --- the object returned by function ODE_Bayesian_Lasso()
                            tcgname --- symbols of TCGs
                         Output objects:
                            eig_scores --- hubscores of each TCG
Time_course.R --- function to draw timecourse curve of top genes
                         Input objects:
                            eig_scores1 --- hubscores of each TCG
                            cut  --- the number of top genes that you want to plot
                            TCG_series1 --- reordered samples against temporal progression of TCGs
                            pseudo_time1 --- sampled time marks in trajectory smoothing
                         Output objects:
                             A plot of timecouse curves of top genes
