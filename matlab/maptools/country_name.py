import numpy as np
import reverse_geocoder as rg

data = np.genfromtxt('locations_volcanic_usgs.csv',delimiter=',')
outdata = []

for line in data:
    id = line[0]
    if len(line) > 1:
        lat = line[1]
        lon = line[2]
        
        g = rg.search([lat,lon],mode=1)
        if not g:
            cntry = 'unknown'
        else:
            cntry = g[0]['cc']
        line = line.tolist()
        line.append(cntry)
        print line
        
    outdata.append(line)

with open('locations_volcanic_usgs.csv','w') as f:
    for data in outdata:
        data[0]=int(data[0])
        line = ','.join([str(i) for i in data])
        f.write(line+'\n')

#np.savetxt('locations_w_cntry.csv',outdata,delimiter=',')


