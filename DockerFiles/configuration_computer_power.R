nWorkers = max(8,parallel::detectCores()-2)  # Define maximum number of workers for parallel processing
nBytesRAM = 4*1024*1024*1024
nBytesRAM_OS = 8*1024*1024*1024
nBytesRAM = max(nBytesRAM, benchmarkme::get_ram()-nBytesRAM_OS)  # Define maximum RAM for shiny
num_max_samples = 1e6
