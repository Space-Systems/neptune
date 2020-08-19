#!python3
import requests
import os

url = 'http://aux.sentinel1.eo.esa.int/POEORB/'
dates = [
    '2019/07/15', 
    '2019/07/16', 
    '2019/07/17', 
    '2019/07/18', 
    '2019/07/19', 
    '2019/07/20', 
    '2019/07/21']
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
    if not os.path.exists(f'{dest}/{filename}'):
        print(f'Downloading {filename} ...')
        # Download validation file
        S1B_file = requests.get(f'{url}/{dates[i]}/{filename}')

        with open(f'{dest}/{filename}', "wb") as file:
            file.write(S1B_file.content)
    else:
        print('Nothing to download.')
    i=i+1