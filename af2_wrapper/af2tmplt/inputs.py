from af2tmplt import global_var

def process_inputs(inputs_dict):
    sequences = inputs_dict.get("sequences", None)
    if sequences == None:
        print("The input json file does not contain the 'sequences' block.")
        exit()
    seq_names_list = []
    seq_sequences_list = []
    seq_chains_list = []
    for i, sequence in enumerate(sequences):
        seq_name = sequence.get("name", None)
        if seq_name == None:
            print(f"No 'name' field specified for sequence {i}")
            exit()
        seq_names_list.append(seq_name)

        seq_sequence = sequence.get("sequence", None)
        if seq_sequence == None:
            print(f"No 'sequence' field specified for sequence {i}")
            exit()
        seq_sequences_list.append(seq_sequence)

        seq_chains = sequence.get("chains", None)
        if seq_chains == None:
            print(f"No 'chains' field specified for sequence {i}")
            exit()
        if len(seq_chains) < 1:
            print(f"The 'chains' field is empty for sequence {i}")
            exit()
        seq_chains_list += seq_chains

    if len(seq_names_list) != len(set(seq_names_list)):
        print("The input json file contains non-unique sequence names")
        exit()
    if len(seq_sequences_list) != len(set(seq_sequences_list)):
        print("The input json file contains non-unique sequence strings")
        exit()
    if len(seq_chains_list) != len(set(seq_chains_list)):
        print("The input json file contains non-unique chain names")
        exit()

    templates = inputs_dict.get("templates", None)
    if templates != None:
        if len(templates) < 1:
            print("The input json file contains the 'templates' block, but it is empty. Will run without templates")
        for i, template in enumerate(templates):
            pdb_file = template.get("pdb_file", None)
            if pdb_file == None:
                print(f"The 'pdb_file' field is empty for template {i}")
                exit()

            query_to_template = template.get("query_to_template", None)
            if query_to_template == None:
                print(f"The 'query_to_template' field is empty for template {i}")
                #exit()
                continue
            map_query_chains_list = []
            map_template_chains_list = []
            for map_query_chain, map_template_chain in query_to_template.items():
                map_query_chains_list.append(map_query_chain)
                map_template_chains_list.append(map_template_chain)
            if len(map_template_chains_list) != len(set(map_template_chains_list)):
                print("The same template chain is used multiple times in template {i}")
                exit()
            for map_query_chain in map_query_chains_list:
                if map_query_chain not in seq_chains_list:
                    print(f"Unrecongnized query chain {map_query_chain} is used in template {i}")
                    exit()
            for seq_chain in seq_chains_list:
                if seq_chain not in map_query_chains_list:
                    print(f"Warning: query chain {seq_chain} is not mapped to template {i}")
                    exit()

    parameters = inputs_dict.get("parameters", None)
    if parameters == None:
        print("The input json file does not contain the 'parameters' field.")
        exit()

    msa = inputs_dict.get("msa", None)
    if msa == None:
        print("The input json file does not contain the 'msa' field.")
        exit()

    job_name = inputs_dict.get("job_name", None)
    if job_name == None:
        print("The input json file does not contain the 'job_name' field.")
        exit()

    num_models = inputs_dict.get("num_models", None)
    if num_models == None:
        print("The input json file does not contain the 'num_models' field. (Set to 5 for default value)")
        exit()

    num_recycle = inputs_dict.get("num_recycle", None)
    if num_recycle == None:
        print("The input json file does not contain the 'num_recycle' field. (Set to 3 for default value)")
        exit()

    num_predictions = inputs_dict.get("num_predictions", None)
    if num_predictions == None:
        print("The input json file does not contain the 'num_predictions' field. (Set to 1 for default value)")
        exit()

    num_msa_clusters = inputs_dict.get("num_msa_clusters", None)
    if num_msa_clusters == None:
        print("The input json file does not contain the 'num_msa_clusters' field. Set to 'default' for default value. (512/252 are monomer/multimer defaults)")
        exit()

    msa_crop_size = inputs_dict.get("msa_crop_size", None)
    if msa_crop_size == None:
        print("The input json file does not contain the 'msa_crop_size' field. (Set to 2048 for default value)")
        exit()

    seed = inputs_dict.get("seed", None)
    if seed == None:
        print("The input json file does not contain the 'seed' field. (Set to 0 for default value)")
        exit()

    hhblits_path = inputs_dict.get("hhblits_path", None)
    if (hhblits_path == None) and (msa == 'standard'):
        print("The input json file does not contain the 'hhblits_path' field.")
        exit()

    hmmer_path = inputs_dict.get("hmmer_path", None)
    if (hmmer_path == None) and (msa == 'standard'):
        print("The input json file does not contain the 'hmmer_path' field")
        exit()

    global_var.MODEL_SET = parameters

    seq_dict = {}
    for sequence in sequences:
        seq_name = sequence['name']
        seq_str = sequence['sequence']
        seq_dict[seq_name] = seq_str

    return seq_dict, seq_chains_list
