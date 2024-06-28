# Sim Generator input file


size=6
pathin='/Users/reneehlozek/Dropbox/CIFAR_Sims/240415/true_nonoise/t_100_ar_0.6_rep_1/'
pathout='/Users/reneehlozek/Dropbox/CIFAR_Sims/240415/pythonsims_t_100_ar_0.6_rep_1/'
num_participants = XXX  
num_iterations = 1
save=True

# Covariance choices
covContempName='randn'
covLaggedName='randn'
measurecovName='diag'
ampContemp = X
ampLagged = X
ampMeasure = X

#Masking
maskZero=True

# Starting value
start = [0,0,0,0,0,0]

# Size of simulation
steps = XXX

# Are you clipping/truncating the time series
clip_samples=False
clip_indices=[0,1]
clip_mins=[0.5, 0.7]
clip_maxs=[0.8,1.3]

clip_outliers=True
clip_sigma=20

debug=False