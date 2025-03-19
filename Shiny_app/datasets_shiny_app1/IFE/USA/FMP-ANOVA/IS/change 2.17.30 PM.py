import os
files = list(filter(lambda x: 'py' not in x, os.listdir()))

def mapper(x):
    x = x.split('.')
    return ''.join(x[:-1]) + '_FMP.rds'

new_names = list(map(mapper, files))

for old, new in zip(files, new_names):
    os.system(f'mv {old} {new}')
