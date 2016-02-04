import subprocess
import logging
import threading
import yaml
import os

FNULL = open(os.devnull, "w")

stream = file("config.yaml", "r")
config = yaml.load(stream)
stream.close()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s (%(threadName)-2s) %(message)s')

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
    config['run_name'] = 'earth_%s' % i
    config['body'] = 'earth'
    config['constant_diff'] = 'true'
    config['constant_diff_value'] = 1e-5*(1+0.5*int(i))
    stream = file('config.yaml', 'w')
    yaml.dump(config, stream)

def worker(s, pool):
    logging.debug('Waiting to join the pool')
    with s:
        # generate config here
        name = threading.currentThread().getName()
        generate_config(name)
        pool.makeActive(name)
        subprocess.call(['./a.out'], stdout=FNULL)
        pool.makeInactive(name)

pool = ActivePool()
s = threading.Semaphore(2)
for i in range(10):
    t = threading.Thread(target=worker, name=str(i), args=(s, pool))
    t.start()
