import reverse_geocoder as rg
import numpy as np

coordinates = (51.5214588,-0.1729636),(9.936033, 76.259952),(37.38605,-122.08385)

results = rg.search(coordinates,mode=1) # default mode = 2

print results