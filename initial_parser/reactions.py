import requests
import re
import pandas as pd
from bs4 import BeautifulSoup

out_path = '/home/efarr/Documents/metalinks/Data/Intermediate/HMDB/hmdb_reactions_test.csv'

def get_PD(reaction_page):

    # parse the html using beautiful soup
    soup = BeautifulSoup(reaction_page, 'html.parser')

    # find all the reaction panels
    reaction = soup.find(class_='reaction-panel')

    # extract all the metabolite ids
    metabolite_ids = re.findall(r'/metabolites/(HMDB\d+)', str(reaction))

    status = reaction.text.split('Status')[1].strip()
    status = status.split(' ')[0]
    #status =  'reaction id: ' + str(reaction_id) + ' has status: ' + status


    enzyme_id = re.search(r'/proteins/(HMDBP\d+)', str(reaction)).group(1)

    reaction_str = soup.find(class_='panel-heading').text
    # split reaction string by either + or = and count how many object were before the =
    reactands, products = reaction_str.split('=')
    reactands = reactands.split('+')

    reactand_ids = metabolite_ids[:len(reactands)]
    product_ids = metabolite_ids[len(reactands):]

    df = pd.DataFrame({'HMDBP': enzyme_id, 
                       'Metabolite': reactand_ids + product_ids, 
                       'Type': ['Reactand'] * len(reactand_ids) + ['Product'] * len(product_ids),
                       'Status': status})

    return df



dfs = []
for reaction_id in range(1, 18203):  # number of reactions in HMDB is 18 203 these days
    # specify the HMDB webfile URL to be scraped
    hmdb_url = 'https://hmdb.ca/reactions/' + str(reaction_id)

    # send a GET request to the HMDB webfile URL and store the response
    response = requests.get(hmdb_url)

    
    try: 
        df = get_PD(response.text)
    except:
        print(f'Could not parse reaction {reaction_id}')
        continue

    print(f'Found {df.shape[0]} metabolites for reaction {reaction_id}')
    
    dfs.append(df)

res = pd.concat(dfs, ignore_index=True).drop_duplicates()

res.to_csv(out_path, index=False)



