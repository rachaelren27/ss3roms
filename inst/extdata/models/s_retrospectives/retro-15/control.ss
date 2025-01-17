#V3.30
#C file created using the SS_writectl function in the R package r4ss
#C file write time: 2022-11-01 16:34:07
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
7 #_Nblock_Patterns
1 5 5 3 3 1 1 #_blocks_per_pattern
#_begin and end years of blocks
1995 2014
1942 1946 1947 1996 1997 2010 2011 2018 2019 2020
1942 1946 1947 1981 1982 2010 2011 2018 2019 2020
1997 2002 2003 2010 2011 2020
1982 2002 2003 2010 2011 2020
1995 2014
1890 1890
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_no additional input for selected M option; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
30 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
3 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
 0.01	    0.11	  0.0725862	-2.93857	0.438	3	  3	0	0	0	0	0	0	0	#_NatM_p_1_Fem_GP_1        
   22	      35	    25.7207	      22	   99	0	  2	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1       
   60	      70	    62.4569	      66	   99	0	  2	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1       
 0.15	    0.55	   0.343282	    0.25	   99	0	  2	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1       
0.001	    0.15	  0.0572535	    0.05	   99	0	  2	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1        
 0.01	     0.3	   0.109531	    0.11	   99	0	  2	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1          
    0	       1	3.31546e-06	       0	   99	6	-50	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1         
    0	       4	    3.27264	     3.3	   99	6	-50	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1         
   53	      59	      55.19	      55	   99	6	-50	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1          
   -3	       3	     -0.421	   -0.25	   99	6	-50	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1       
   -3	       3	          1	       1	   99	6	-50	0	0	0	0	0	0	0	#_Eggs/kg_inter_Fem_GP_1   
   -3	       3	          0	       0	   99	6	-50	0	0	0	0	0	0	0	#_Eggs/kg_slope_wt_Fem_GP_1
 0.01	    0.11	  0.0604722	-2.89857	0.438	3	  3	0	0	0	0	0	0	0	#_NatM_p_1_Mal_GP_1        
   15	      35	     26.926	       0	   99	0	  2	0	0	0	0	0	0	0	#_L_at_Amin_Mal_GP_1       
   50	      60	    56.6228	       0	   99	0	  2	0	0	0	0	0	0	0	#_L_at_Amax_Mal_GP_1       
  0.2	    0.55	   0.371287	       0	   99	0	  2	0	0	0	0	0	0	0	#_VonBert_K_Mal_GP_1       
0.001	    0.15	  0.0749236	       0	   99	0	  2	0	0	0	0	0	0	0	#_CV_young_Mal_GP_1        
 0.01	     0.3	  0.0783725	       0	   99	0	  2	0	0	0	0	0	0	0	#_CV_old_Mal_GP_1          
    0	       1	3.37087e-06	       0	   99	6	-50	0	0	0	0	0	0	0	#_Wtlen_1_Mal_GP_1         
    0	       4	    3.27008	     3.3	   99	6	-50	0	0	0	0	0	0	0	#_Wtlen_2_Mal_GP_1         
  0.1	      10	          1	       1	    1	0	 -1	0	0	0	0	0	0	0	#_CohortGrowDev            
1e-06	0.999999	        0.5	     0.5	  0.5	0	-99	0	0	0	0	0	0	0	#_FracFemale_GP_1          
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
  8	 12	9.70454	9.8	   99	0	  1	0	0	0	0	0	0	0	#_SR_LN(R0)  
0.2	  1	    0.7	0.6	0.223	2	 -7	0	0	0	0	0	0	0	#_SR_BH_steep
0.2	1.5	    1.4	0.6	   99	0	-50	0	0	0	0	0	0	0	#_SR_sigmaR  
 -1	  1	      0	  0	   99	0	-50	0	0	0	0	0	0	0	#_SR_regime  
 -1	  1	      0	  0	   99	0	-50	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1925 # first year of main recr_devs; early devs can preceed this era
2020 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase
1 # (0/1) to read 13 advanced options
1860 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
3 #_recdev_early_phase
3 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1976 #_last_yr_nobias_adj_in_MPD; begin of ramp
1980 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2019 #_last_yr_fullbias_adj_in_MPD
2020 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.93 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-4 #min rec_dev
4 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.02 # F ballpark
-2000 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
3 # max F or harvest rate, depends on F_Method
4 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
    4	1	0	1	0	0	#_AKSHLF    
    5	1	0	1	0	1	#_AKSLP     
    6	1	0	1	0	1	#_NWSLP     
    7	1	0	1	0	1	#_NWCBO     
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  -15	 15	 0.325075	0	99	0	 1	0	0	0	0	0	1	2	#_LnQ_base_AKSHLF(5) 
0.025	1.3	 0.178521	0	99	0	 2	0	0	0	0	0	0	0	#_Q_extraSD_AKSHLF(5)
  -15	  5	-0.417466	0	 1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_AKSLP(6)  
0.001	0.7	0.0361518	0	99	0	 2	0	0	0	0	0	0	0	#_Q_extraSD_AKSLP(6) 
  -15	 15	-0.735051	0	 1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_NWSLP(7)  
0.001	0.8	 0.162325	0	99	0	 2	0	0	0	0	0	0	0	#_Q_extraSD_NWSLP(7) 
  -15	 15	 -0.52908	0	 1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_NWCBO(8)  
0.001	0.4	        0	0	99	0	-2	0	0	0	0	0	0	0	#_Q_extraSD_NWCBO(8) 
# timevary Q parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
-15	15	-0.0716196	0	99	0	2	#_LnQ_base_AKSHLF(5)_BLK1repl_1995
# info on dev vectors created for Q parms are reported with other devs after tag parameter section
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
0	2	0	0	#_1 FIX   
0	0	0	0	#_2 NONE  
0	2	0	0	#_3 TWL   
0	0	0	0	#_4 AKSHLF
0	0	0	0	#_5 AKSLP 
0	0	0	0	#_6 NWSLP 
0	0	0	0	#_7 NWCBO 
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
20	0	1	0	#_1 FIX   
 0	0	0	0	#_2 NONE  
20	0	0	0	#_3 TWL   
20	0	1	0	#_4 AKSHLF
20	0	0	0	#_5 AKSLP 
20	0	0	0	#_6 NWSLP 
20	0	0	0	#_7 NWCBO 
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
   10	 60	     41	     30	99	0	 -5	0	0	0	0	0	2	2	#_SizeSel_PRet_1_FIX(1)
  0.1	 20	6.00519	      1	99	0	 -5	0	0	0	0	0	0	0	#_SizeSel_PRet_2_FIX(1)
  -10	 10	     10	2.18163	99	0	 -5	0	0	0	0	0	2	2	#_SizeSel_PRet_3_FIX(1)
  -10	 10	      0	      0	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PRet_4_FIX(1)
    8	 70	     28	     18	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_1_FIX(1)
0.001	  2	   0.01	      1	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_2_FIX(1)
 0.01	0.8	    0.2	    0.1	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_3_FIX(1)
  -10	 10	      0	      0	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_4_FIX(1)
   15	 55	     41	     32	99	0	 -5	0	0	0	0	0	3	2	#_SizeSel_PRet_1_TWL(3)
  0.1	 20	2.89822	      1	99	0	 -5	0	0	0	0	0	0	0	#_SizeSel_PRet_2_TWL(3)
  -10	 10	     10	     10	99	0	 -5	0	0	0	0	0	3	2	#_SizeSel_PRet_3_TWL(3)
  -10	 10	      0	      0	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PRet_4_TWL(3)
    8	 70	     28	     18	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_1_TWL(3)
0.001	  2	   0.01	      1	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_2_TWL(3)
  0.1	0.8	    0.5	    0.5	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_3_TWL(3)
  -10	 10	      0	      0	99	0	-50	0	0	0	0	0	0	0	#_SizeSel_PDis_4_TWL(3)
#_AgeSelex
   2	20	        5	  1	99	0	-4	0	0	0	0	0	4	2	#_AgeSel_P_1_FIX(1)       
 -20	 5	       -4	0.3	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_2_FIX(1)       
 -15	10	 0.194766	  5	99	0	 4	0	0	0	0	0	4	2	#_AgeSel_P_3_FIX(1)       
 -10	10	  2.83982	  4	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_4_FIX(1)       
  -5	 5	       -5	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_5_FIX(1)       
  -5	 5	     -1.5	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_6_FIX(1)       
 -15	15	        0	  0	 0	0	-4	0	0	0	0	0	0	0	#_AgeSel_PMale_1_FIX(1)   
 -15	15	0.0567747	  0	 0	0	 4	0	0	0	0	0	0	0	#_AgeSel_PMale_2_FIX(1)   
 -15	15	-0.835168	  0	 0	0	 4	0	0	0	0	0	0	0	#_AgeSel_PMale_3_FIX(1)   
 -15	15	 -1.31113	  0	 0	0	 4	0	0	0	0	0	0	0	#_AgeSel_PMale_4_FIX(1)   
0.01	20	        1	  1	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_1_TWL(3)       
 -20	 5	       -4	0.3	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_2_TWL(3)       
 -20	10	 -2.40262	  5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_3_TWL(3)       
 -10	10	       -9	  4	99	0	-4	0	0	0	0	0	5	2	#_AgeSel_P_4_TWL(3)       
  -5	 5	 -4.02671	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_5_TWL(3)       
  -5	 5	 -1.60114	 -5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_6_TWL(3)       
   1	12	        1	  1	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_1_AKSHLF(5)    
  -5	 5	       -4	0.3	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_2_AKSHLF(5)    
 -10	10	 -9.74336	  5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_3_AKSHLF(5)    
 -10	10	 -1.01065	  4	99	0	 4	0	0	0	0	0	6	2	#_AgeSel_P_4_AKSHLF(5)    
 -10	 5	     -2.5	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_5_AKSHLF(5)    
 -10	 5	  -3.8559	 -5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_6_AKSHLF(5)    
 -15	15	        0	  0	 0	0	-4	0	0	0	0	0	0	0	#_AgeSel_PMale_1_AKSHLF(5)
 -15	15	-0.544676	  0	 0	0	 4	0	0	0	0	0	0	0	#_AgeSel_PMale_2_AKSHLF(5)
 -15	15	-0.170609	  0	 0	0	 4	0	0	0	0	0	0	0	#_AgeSel_PMale_3_AKSHLF(5)
 -15	15	 -6.15559	  0	 0	0	 4	0	0	0	0	0	0	0	#_AgeSel_PMale_4_AKSHLF(5)
   1	12	  1.46677	  1	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_1_AKSLP(6)     
 -20	 5	       -4	0.3	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_2_AKSLP(6)     
 -10	10	       -4	  5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_3_AKSLP(6)     
 -20	10	 -5.96505	 -4	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_4_AKSLP(6)     
  -5	 5	 -1.33799	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_5_AKSLP(6)     
  -5	 5	-0.530103	 -5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_6_AKSLP(6)     
   1	12	  3.58736	  1	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_1_NWSLP(7)     
  -5	 5	       -4	0.3	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_2_NWSLP(7)     
 -10	10	  1.49463	  5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_3_NWSLP(7)     
 -20	50	 -3.29711	  4	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_4_NWSLP(7)     
  -5	 5	   -4.565	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_5_NWSLP(7)     
  -5	 5	 0.189149	 -5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_6_NWSLP(7)     
0.01	 5	0.0909168	  1	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_1_NWCBO(8)     
 -20	 5	       -4	0.3	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_2_NWCBO(8)     
 -20	10	  -8.4499	  5	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_3_NWCBO(8)     
 -10	10	  3.47597	  4	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_4_NWCBO(8)     
  -5	 5	       -4	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_5_NWCBO(8)     
  -5	 5	-0.320292	 -5	99	0	-4	0	0	0	0	0	0	0	#_AgeSel_P_6_NWCBO(8)     
# timevary selex parameters 
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
 25	60	     25	30	99	0	-5	#_SizeSel_PRet_1_FIX(1)_BLK2repl_1942
 25	60	38.9601	30	99	0	-5	#_SizeSel_PRet_1_FIX(1)_BLK2repl_1947
 25	60	 40.349	30	99	0	 5	#_SizeSel_PRet_1_FIX(1)_BLK2repl_1997
 25	60	41.3726	30	99	0	 5	#_SizeSel_PRet_1_FIX(1)_BLK2repl_2011
 10	60	 35.921	30	99	0	 5	#_SizeSel_PRet_1_FIX(1)_BLK2repl_2019
-10	10	     10	10	99	0	-5	#_SizeSel_PRet_3_FIX(1)_BLK2repl_1942
-10	10	     10	10	99	0	-5	#_SizeSel_PRet_3_FIX(1)_BLK2repl_1947
-10	10	2.53757	10	99	0	 5	#_SizeSel_PRet_3_FIX(1)_BLK2repl_1997
-10	10	4.00765	10	99	0	-5	#_SizeSel_PRet_3_FIX(1)_BLK2repl_2011
-10	10	2.13517	10	99	0	 5	#_SizeSel_PRet_3_FIX(1)_BLK2repl_2019
 15	55	     25	32	99	0	-5	#_SizeSel_PRet_1_TWL(3)_BLK3repl_1942
 15	55	45.9292	32	99	0	-5	#_SizeSel_PRet_1_TWL(3)_BLK3repl_1947
 15	55	47.7536	32	99	0	 5	#_SizeSel_PRet_1_TWL(3)_BLK3repl_1982
 15	55	33.7502	32	99	0	 5	#_SizeSel_PRet_1_TWL(3)_BLK3repl_2011
 15	55	42.2739	32	99	0	 5	#_SizeSel_PRet_1_TWL(3)_BLK3repl_2019
-10	10	     10	10	99	0	-5	#_SizeSel_PRet_3_TWL(3)_BLK3repl_1942
-10	10	     10	10	99	0	-5	#_SizeSel_PRet_3_TWL(3)_BLK3repl_1947
-10	10	 3.7424	10	99	0	 5	#_SizeSel_PRet_3_TWL(3)_BLK3repl_1982
-10	10	     10	10	99	0	-5	#_SizeSel_PRet_3_TWL(3)_BLK3repl_2011
-10	10	5.32636	10	99	0	 5	#_SizeSel_PRet_3_TWL(3)_BLK3repl_2019
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
  2	20	 3.15072	1	99	0	 4	#_AgeSel_P_1_FIX(1)_BLK4repl_1997   
  2	20	 5.03792	1	99	0	 4	#_AgeSel_P_1_FIX(1)_BLK4repl_2003   
  2	20	 3.05793	1	99	0	 4	#_AgeSel_P_1_FIX(1)_BLK4repl_2011   
-10	20	-1.23983	4	99	0	-4	#_AgeSel_P_3_FIX(1)_BLK4repl_1997   
-10	20	   1.852	4	99	0	 4	#_AgeSel_P_3_FIX(1)_BLK4repl_2003   
-10	20	-8.20053	4	99	0	 4	#_AgeSel_P_3_FIX(1)_BLK4repl_2011   
-10	10	  2.0564	4	99	0	 4	#_AgeSel_P_4_TWL(3)_BLK5repl_1982   
-10	10	   6.599	4	99	0	 4	#_AgeSel_P_4_TWL(3)_BLK5repl_2003   
-10	10	 9.17751	4	99	0	 4	#_AgeSel_P_4_TWL(3)_BLK5repl_2011   
-10	10	 3.17048	4	99	0	 4	#_AgeSel_P_4_AKSHLF(5)_BLK6repl_1995
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Factor	Fleet	Value
    2	1	       0	#_Variance_adjustment_list1 
    2	3	       0	#_Variance_adjustment_list2 
    3	1	0.302786	#_Variance_adjustment_list3 
    3	3	       0	#_Variance_adjustment_list4 
    4	1	0.095328	#_Variance_adjustment_list5 
    4	3	0.044144	#_Variance_adjustment_list6 
    4	7	0.032931	#_Variance_adjustment_list7 
    5	1	0.101402	#_Variance_adjustment_list8 
    5	3	0.193659	#_Variance_adjustment_list9 
    5	4	       1	#_Variance_adjustment_list10
    5	5	0.109196	#_Variance_adjustment_list11
    5	6	 0.12705	#_Variance_adjustment_list12
    5	7	0.286539	#_Variance_adjustment_list13
-9999	0	       0	#_terminator                
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
-9999 0 0 0 0 # terminator
#
0 # 0/1 read specs for more stddev reporting
#
999
