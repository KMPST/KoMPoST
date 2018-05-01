# This script:
# 1) takes input parameters for the initial condition ansatz Trento
# 2) writes an config file for Trento
# 3) runs Trento
# 4) converts the Trento output file in a format that can be read as initial conditions by KoMPoST
# 5) write an input file for KoMPoST that can then be run

import configparser
import os.path
import subprocess
import sys
import shutil

##################################################################
################### 1) Trento input parameters ###################
##################################################################

# Grid used to save Trento and also used later to run KoMPoST
# Total number of grid cell
grid_size=512
# Grid step in fermi
dx=0.08

# Type of ion (Au, Pb, ...)
trento_nucleus='Pb'

# Python dictionary containing Trento parameters
trento_parameters = {
'normalization':10,
'reduced-thickness':0.1,
'fluctuation':0.9,
'nucleon-min-dist':0.4,
'cross-section':6.4, # in fermi
'nucleon-width':0.8
}


##################################################################
################## 2) Writing Trento config file ##################
##################################################################

trento_config_file="trento_config"

if (os.path.isfile(trento_config_file)):
        print("I'd rather not overwrite an existing file (\""+trento_config_file+"\")...  Aborting")
        exit(1)

trento_tmp_dir=trento_nucleus+trento_nucleus

with open(trento_config_file, 'w') as trento_config:
    trento_config.write('projectile = '+trento_nucleus+'\n')
    trento_config.write('projectile = '+trento_nucleus+'\n')
    trento_config.write('number-events = 1\n')
    trento_config.write('output = '+trento_tmp_dir+'\n')
    trento_config.write('grid-max = '+str(grid_size*dx/2.)+'\n')
    trento_config.write('grid-step = '+str(dx)+'\n')
    for param, value in trento_parameters.items():
            trento_config.write(param+' = '+str(value)+'\n')


##################################################################
####################### 3) Running Trento ########################
##################################################################

# Location of Trento executable
trento_bin='trento'

# Run Trento
# subprocess.run(args, *, stdin=None, input=None, stdout=None, stderr=None, shell=False, cwd=None, timeout=None, check=False, encoding=None, errors=None)
#print(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),trento_config_file))
#subprocess.Popen([trento_bin, "-c",os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),trento_config_file)], shell=True)
os.system(trento_bin+" -c "+os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),trento_config_file))


#shutil.move(os.path.join(trento_tmp_dir,"0.dat"),"trento_raw_output.dat")

##################################################################
############### 4) Trento output to KoMPoST intput ###############
##################################################################


kompost_init_name="trento_initial_condition_in_kompost_format.txt"

if (os.path.isfile(kompost_init_name)):
        print("I'd rather not overwrite an existing file (\""+kompost_init_name+"\")...  Aborting")
        exit(1)

# Use Soeren's bash script, which we know works
conversion_bash_command="Ns=\""+str(grid_size)+"\"; tail -n+9 "+os.path.join(trento_tmp_dir,"0.dat") +" | perl -pe 's/\s+/\\n/g' | awk -v Ns=${Ns} 'BEGIN {for(x=0;x<Ns;x++){for(y=0;y<Ns;y++){T00[x+Ns*y]=0.0;}} Index=0;} {T00[Index]=$1; Index++;} END {for(y=0;y<Ns;y++){for(x=0;x<Ns;x++){print x,y,T00[x+Ns*y],0.5*T00[x+Ns*y],0.5*T00[x+Ns*y],0.0,0.0,0.0,0.0,0.0,0.0,0.0;} printf(\"\\n\"); }}' > "+kompost_init_name

#print(conversion_bash_command)

subprocess.Popen(conversion_bash_command, shell=True)


##################################################################
################## 5) Writing KoMPoST input file #################
##################################################################

setup_file_template="setup_trento_template.ini"

setup_file_updated="setup_trento.ini"

# Open setup file
config = configparser.ConfigParser()
config.read(setup_file_template)

# 
config['EventInput']['Ns']=str(grid_size)
config['EventInput']['afm']=str(dx)
config['EventInput']['xSTART']=str(0)
config['EventInput']['xEND']=str(grid_size-1)
config['EventInput']['ySTART']=str(0)
config['EventInput']['yEND']=str(grid_size-1)
config['KoMPoSTInputs']['inputfile']=os.path.realpath(kompost_init_name)

# Write back setup file
with open(setup_file_updated, 'w') as configfile:
    config.write(configfile)
