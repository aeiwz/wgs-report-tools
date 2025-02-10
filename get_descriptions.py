from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, WebDriverException
from tqdm import tqdm
import pandas as pd
import time


ols = ['EFO', 'MONDO', 'Orphanet']
amigo = ['GO']
ontobee = ['OBA', 'NCIT', 'HANCESTRO', 'PATO', 'MP']
hpo = ['HP']


# Load the data
df = pd.read_csv('data/gwas_database.csv.gz', compression='gzip')

# Initialize storage
url_ = []
description = []
save_interval = 30  # Save progress every 30 links
output_file = 'data/gwas_database_with_description.csv.gz'

# Set up Chrome options for headless mode
chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--disable-gpu')
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--disable-dev-shm-usage')

# Get unique links to process
unique_links = df['MAPPED_TRAIT_URI'].unique()

# Resume progress if needed
try:
    df_existing = pd.read_csv(output_file, compression='gzip')
    processed_links = set(df_existing['MAPPED_TRAIT_URI'].dropna())
    print(f"Resuming from previous progress. {len(processed_links)} links already processed.")
except FileNotFoundError:
    processed_links = set()

# Start processing
for index, link in enumerate(tqdm(unique_links)):
    if link in processed_links or link == 'No URI':
        description.append('Description not available')
        continue
    
    try:
        # Create a new Chrome driver instance
        driver = webdriver.Chrome(options=chrome_options)
        driver.get(link)
        driver.implicitly_wait(30)
        
        if link.split('/')[-1].split('_')[0] in ols:
            p_tag = driver.find_element(By.CSS_SELECTOR, "p.pb-3")
            text = p_tag.text
        elif link.split('/')[-1].split('_')[0] in amigo:
            p_tag = driver.find_element(By.XPATH, "//dt[contains(text(), 'Definition')]/following-sibling::dd")
            text = p_tag.text
        
        url_.append(link)
        description.append(text)
    
    except NoSuchElementException:
        description.append('Description not found')
    except WebDriverException:
        description.append('WebDriver error occurred')
    except Exception as e:
        description.append(f'Error: {str(e)}')
    
    finally:
        driver.quit()
    
    # Save progress every `save_interval` links
    if (index + 1) % save_interval == 0:
        temp_dict = dict(zip(url_, description))
        df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(temp_dict).fillna(df.get('MAPPED_TRAIT_DESCRIPTION'))
        df.to_csv(output_file, index=False, compression='gzip')
        print(f"Saved progress at {index + 1} links.")
        time.sleep(5)  # Pause to avoid overloading the server

# Final save
temp_dict = dict(zip(url_, description))
df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(temp_dict).fillna(df.get('MAPPED_TRAIT_DESCRIPTION'))
df.to_csv(output_file, index=False, compression='gzip')
print("Completed fetching descriptions and saved the final file.")