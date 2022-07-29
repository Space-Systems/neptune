#!python3
import requests
import os
from zipfile import ZipFile

url = 'http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/S1B/'
dates = [
    '2019/06', 
    '2019/06', 
    '2019/06', 
    '2019/06', 
    '2019/06', 
    '2019/06', 
    '2019/06']
filenames = [
    'S1B_OPER_AUX_POEORB_OPOD_20190715T110507_V20190624T225942_20190626T005942.EOF',
    'S1B_OPER_AUX_POEORB_OPOD_20190716T110556_V20190625T225942_20190627T005942.EOF',
    'S1B_OPER_AUX_POEORB_OPOD_20190717T110525_V20190626T225942_20190628T005942.EOF',
    'S1B_OPER_AUX_POEORB_OPOD_20190718T110707_V20190627T225942_20190629T005942.EOF',
    'S1B_OPER_AUX_POEORB_OPOD_20190719T110400_V20190628T225942_20190630T005942.EOF',
    'S1B_OPER_AUX_POEORB_OPOD_20190720T110629_V20190629T225942_20190701T005942.EOF',
    'S1B_OPER_AUX_POEORB_OPOD_20190721T110515_V20190630T225942_20190702T005942.EOF']
dest = 'poeData'

with open(f'{dest}/available_files.txt', 'w') as file:
    for filename in filenames:
        file.write('%s\n' % filename)
i = 0
for filename in filenames:
    if not os.path.exists(f'{dest}/{filename}.zip'):
        print(f'Downloading {url}/{dates[i]}/{filename}.zip ...')
        # Download validation file
        
        S1B_file = requests.get(f'{url}/{dates[i]}/{filename}.zip')

        with open(f'{dest}/{filename}.zip', "wb") as file:
            file.write(S1B_file.content)
    if os.path.exists(f'{dest}/{filename}.zip') \
            and not os.path.exists(f'{dest}/{filename}'):
        with ZipFile(f'{dest}/{filename}.zip', 'r') as zipObject:
            # Extract from tmp/<filename> to destination
            zipInfo = zipObject.getinfo(f'tmp/{filename}')
            # Manipulate the filename to drop the tmp/ path
            zipInfo.filename  = os.path.basename(f'{filename}')
            zipObject.extract(zipInfo, f'{dest}')
    else:
        print('Nothing to download.')
    i=i+1