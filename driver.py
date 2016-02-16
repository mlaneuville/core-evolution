import subprocess
import time
import logging
import threading
import yaml
import os




class ActivePool(object):
    def __init__(self):
        super(ActivePool, self).__init__()
        self.active = []
        self.lock = threading.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
            logging.debug('Running: %s', self.active)
    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)
            logging.debug('Running: %s', self.active)

def generate_config(i):
    '''Example of config generation. Anything can be done here, the point is just
    to generate a new config.yaml file for each run.'''

    ## # if varying mantle temperature
    ## TM = 4000+100*int(i)
    ## config['run_name'] = 'earth-Tc%d' % TM
    ## config['mantle_temperature'] = TM

    diffusivity = (0.5+int(i)/2.)*1.e-5
    config['run_name'] =  'earth-diff%.1e' %diffusivity
    config["constant_diff"] = "True"
    config["constant_diff_value"] = diffusivity
    config["tbl_conductivity"] = 10. 

    stream = file('config.yaml', 'w')
    yaml.dump(config, stream)

def worker(s, pool):
    name = threading.currentThread().getName()
    logging.debug('Worker %02d waiting to join the pool' % int(name))
    with s:
        generate_config(name)
        pool.makeActive(name)
        subprocess.call(['./a.out'], stdout=FNULL)
        pool.makeInactive(name)

        



# batch run information
PARALLEL_JOBS = 3   # max number of parallel jobs
TOTAL_RUNS = 11     # total number of runs in the series

# null pointer to discard direct code output
FNULL = open(os.devnull, "w")

# loads default config file structure

stream = file("default.yaml", "r")
config = yaml.load(stream)
stream.close()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s (%(threadName)-2s) %(message)s')



pool = ActivePool()
s = threading.Semaphore(PARALLEL_JOBS)
for i in range(TOTAL_RUNS):
    t = threading.Thread(target=worker, name=str(i), args=(s, pool))
    time.sleep(0.2) # to avoid workers overwriting each other's config file
    t.start()
