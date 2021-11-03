import re
import sqlite3
import numpy as np
import pandas as pd
from skimage.io import imread, imsave
from skimage.draw import rectangle

DB_PATH_TPL = 'inputs/locations/sqlite/{plate_id}.sqlite'
QUERY = '''
SELECT
    Image_FileName_OrigDNA,
    CAST(Nuclei_Location_Center_X AS INTEGER),
    CAST(Nuclei_Location_Center_Y AS INTEGER)
FROM
    Nuclei
JOIN Image ON
    Image.TableNumber = Nuclei.TableNumber
WHERE
    Image_FileName_OrigDNA like :pattern
'''
PATTERN_TPL = '%_{well_id}_s{site_id}%'


def get_locations(plate_id, well_id, site_id):
    locations = []
    dbpath = DB_PATH_TPL.format(plate_id=plate_id)
    with sqlite3.connect(dbpath) as conn:
        print(f'connected to {dbpath}')
        conn.set_trace_callback(print)
        cur = conn.cursor()
        param = PATTERN_TPL.format(well_id=well_id, site_id=site_id)
        result = cur.execute(QUERY, {'pattern': param})
        for _, xpos, ypos in result:
            locations.append((xpos, ypos))
    return np.asarray(locations)


im_root = 'outputs/compressed/images'
index = pd.read_csv('inputs/metadata/index.csv')
sample = index.sample(10, random_state=123)

fname_rgx = re.compile(r'(\d+)/(.*)_(\w\d+)_s(\d)_.*')

for i, row in sample.iterrows():
    plate_id, _, well_id, site_id = fname_rgx.match(row['DNA']).groups()
    impath = im_root + '/' + row['DNA'].split('.')[0] + '.png'
    im = imread(impath)
    print(f'getting centers from {impath}')
    locations = get_locations(plate_id, well_id, site_id)
    for loc in locations:
        topleft = loc[::-1] - 5
        rr, cc = rectangle(topleft, extent=(10, 10), shape=im.shape)
        im[rr, cc] = 255
    imsave(f'{i}.png', im)
