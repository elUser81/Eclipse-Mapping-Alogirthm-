from eclipse_mapping.eclipse_mapping_algorithm import Processes, Binary

binary_data_dir = 'path/to/folder/containing/flux data/and/parameter file'

b = Processes.load_binary(binary_data_dir)
Processes.build_compute_save(binary_data_dir)