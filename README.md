# frun

To setup please go into the env_setup dir and run the setup_env.sh script


to run the project activate the alphafold2_env conda environment 

for single runs

python complete_run.py -l /PATH/TO/LIGAND/A3M -r /PATH/TO/RECEPTOR/A3M -c /PATH/TO/COMPLEX/PDB -o OUTPUT_DIR -n 

if you are running a nanobody include the -n flag
otherwise include 2 -l paths for the H and L a3ms


for multiple runs 

python running_multiple_complexes.py -l /PATH/TO/LIGAND/A3M -r /PATH/TO/RECEPTOR/A3M -c /PATH/TO/COMPLEX/DIR -o OUTPUT_DIR -n

the complex dir will contain the complexes you want to run, it will run on all complex.pdbs in that directory
the same rules apply for the multiple runs regarding nanobodies and antibodies
