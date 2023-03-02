# Sim Generator input file


size=6
pathin='/Users/reneehlozek/Code/CIFAR_Network/GIMME/220429.FinalMatrices/'
pathout='/Users/reneehlozek/Dropbox/CIFAR_Sims/'
num_participants=100
num_iterations=10
save=False

# Covariance choices
covContempName='randn'
covLaggedName='randn'
measurecovName='diag'
ampContemp = 0.1
ampLagged = 0.1
ampMeasure = 1.0

#Masking
maskZero=True

# Starting value
start = [0,0,0,0,0,0]

# Size of simulation
steps=100

# Are you clipping/truncating the time series
clip_samples=False
clip_indices=[0,1]
clip_mins=[0.5, 0.7]
clip_maxs=[0.8,1.3]

