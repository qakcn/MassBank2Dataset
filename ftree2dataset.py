from functions import *

ftreepath = Path('../data/sirius/ftrees')
compounds_file = Path('../data/sirius/ms/all.MS2.pkl')

outputpath = Path('./outputs')

# parse_ftrees(ftreepath, compounds_file, outputpath)

gen_dataset(outputpath/'sliced', True)
