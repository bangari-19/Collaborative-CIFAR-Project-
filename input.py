# Sim Generator input file


size=6
pathin='/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_syncrng_replaceNF/true_noisefree/t_200_ar_0.6_rep_1/'
pathout='/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_syncrng_replaceNF/pythondata/t_200_ar_0.6_rep_1/'
num_participants=150  
num_iterations = 1
save=True

# Covariance choices
covContempName='randn'
covLaggedName='randn'
measurecovName='diag'
ampContemp = 0.01
ampLagged = 0.01
ampMeasure = 1.0

#Masking
maskZero = True

# Starting value
start = [0,0,0,0,0,0]

# Size of simulation
steps=200

# Are you clipping/truncating the time series
clip_samples=False
clip_indices=[0,1]
clip_mins=[0.5, 0.7]
clip_maxs=[0.8,1.3]

clip_outliers=True
clip_sigma=20

debug=False