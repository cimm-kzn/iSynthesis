from CGRdb import load_schema


EPSILON = 0.05
ROOT = 'root'
VISITS = 'visits'
SCORE = 'score'
VALUE = 'value'
DATA = 'data'
DEPTH = 'depth'
MEAN = 'mean_value'

DB_Name = 'bb'
DEEP = False
THREADS = 20
SEED = 11

db = load_schema('zinc', user='postgres', password='password', host='localhost', database='postgres', port=2202)


FPS = 'iSynthesis/data/zinc/zinc_fps.pickle'
