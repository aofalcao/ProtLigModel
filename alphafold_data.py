import requests, sys, json



def get_url(url, **kwargs):
    
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response



def get_pdb_files(ids):
    
    for i in ids:
        try:
            r = get_url('https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'.format(i))
            with open('pdb_files/{}.pdb'.format(i), 'w') as fp:
                fp.write(r.text)
                
        except:
            print(i, 'not in AlphaFold')
