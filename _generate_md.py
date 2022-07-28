'''
Script to generate .md files for each allele
'''
import json
import os

with open('_data/alleles.json') as fp:
    alleles = json.load(fp)

template = '''---
layout: default
---
{{% assign variant=site.data.alleles[{page_index}] %}}
{{% include variant.html %}}
'''

for allele in alleles:
    target_dir = f'variants/{allele["variant"]}'
    os.makedirs(target_dir, exist_ok=True)
    with open(f'{target_dir}/index.html', 'w') as f:
        f.write(template.format(page_index=allele['index']))
