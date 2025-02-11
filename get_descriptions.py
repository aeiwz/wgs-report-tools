from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, WebDriverException
from tqdm import tqdm
import pandas as pd
import time
import random
import os

# Define ontology groups
ols = ['EFO', 'MONDO', 'Orphanet']
amigo = ['GO']
hpo = ['HP']

# Load dataset
df = pd.read_csv('data/gwas_database.csv.gz', compression='gzip', low_memory=False)

# Output file path
output_file = 'data/gwas_database_with_description.csv.gz'

# List of error messages that require reprocessing
list_of_reprocess = [
    'Ontology not covered', 'Description not available', 
    'Description not found', 'WebDriver error occurred', 
    'Error fetching description'
]

# Load existing data if available
if os.path.exists(output_file):
    df_existing = pd.read_csv(output_file, compression='gzip', low_memory=False)
    if 'MAPPED_TRAIT_DESCRIPTION' in df_existing.columns:
        df['MAPPED_TRAIT_DESCRIPTION'] = df_existing['MAPPED_TRAIT_DESCRIPTION']
    else:
        df['MAPPED_TRAIT_DESCRIPTION'] = ''
else:
    df['MAPPED_TRAIT_DESCRIPTION'] = ''

# Identify traits that need reprocessing
to_process = df[df['MAPPED_TRAIT_DESCRIPTION'].isna() | df['MAPPED_TRAIT_DESCRIPTION'].isin(list_of_reprocess)]
unique_links = to_process['MAPPED_TRAIT_URI'].dropna().unique()

print(f"Found {len(unique_links)} traits to process.")

# Set up Chrome options
chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--disable-gpu')
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--disable-dev-shm-usage')

# Dictionary for storing descriptions
trait_descriptions = {}

# Start WebDriver
driver = webdriver.Chrome(options=chrome_options)

# Process each link
for index, link in enumerate(tqdm(unique_links)):
    if link == 'No URI':
        trait_descriptions[link] = 'Description not available'
        continue

    try:
        if len(link.split(',')) > 1:
            link = link.split(',')[0]
        driver.get(link)
        driver.implicitly_wait(60)

        ontology = link.split('/')[-1].split('_')[0]

        if ontology in ols:
            text = driver.find_element(By.CSS_SELECTOR, "p.pb-3").text
        elif ontology in amigo:
            text = driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div[2]/dl/dd[6]").text
        elif ontology in hpo:
            text = driver.find_element(By.XPATH, "/html/body/app-root/mat-sidenav-container/mat-sidenav-content/div/app-term/div/div/div/div[2]/mat-card/mat-card-content/div/div[2]/p").text
        else:
            text = 'Ontology not covered'

        trait_descriptions[link] = text

    except NoSuchElementException:
        trait_descriptions[link] = 'Description not found'
    except WebDriverException:
        trait_descriptions[link] = 'WebDriver error occurred'
    except Exception as e:
        trait_descriptions[link] = f'Error fetching description'

    # Randomized sleep to avoid rate limiting (between 2-5 seconds)
    time.sleep(random.uniform(2, 5))

    # Save progress every 30 links
    if (index + 1) % 30 == 0:
        df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(trait_descriptions).fillna(df['MAPPED_TRAIT_DESCRIPTION'])
        df.to_csv(output_file, index=False, compression='gzip')
        print(f"Saved progress at {index + 1} links.")

# Save final data
df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(trait_descriptions).fillna(df['MAPPED_TRAIT_DESCRIPTION'])
df.to_csv(output_file, index=False, compression='gzip')
print("Completed processing and saved data.")

# Quit WebDriver
driver.quit()