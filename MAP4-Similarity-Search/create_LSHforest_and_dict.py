import pickle
import random
import sys

import numpy as np
import tmap as tm

fps = []
findCID = {}

with open(sys.argv[1]) as inFile:
    for i, line in enumerate(inFile):
        line = line.strip()
        line = line.split(' ')
        findCID[i] = [line[0], line[1], line[2]]
        fp = tm.VectorUint(np.array(list(map(int, line[3].split(';')))))
        fps.append(fp)

pickle.dump(findCID, open('{}_dictionary'.format(sys.argv[1]), 'wb'))


lf = tm.LSHForest(512, 32)
lf.batch_add(fps)

lf.index()
lf.store('{}_LSHforest'.format(sys.argv[1]))
