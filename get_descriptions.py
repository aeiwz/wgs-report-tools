from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, WebDriverException
from tqdm import tqdm
import pandas as pd
import time
import random

# Define ontology groups
ols = ['EFO', 'MONDO', 'Orphanet']
amigo = ['GO']
ontobee = ['OBA', 'NCIT', 'HANCESTRO', 'PATO', 'MP']
hpo = ['HP']

# Load the dataset
df = pd.read_csv('data/gwas_database.csv.gz', compression='gzip')

# Output file path
output_file = 'data/gwas_database_with_description.csv.gz'

# Set up Chrome options for headless mode
chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--disable-gpu')
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--disable-dev-shm-usage')

# Get unique links and remove processed ones
unique_links = df['MAPPED_TRAIT_URI'].dropna().unique()
try:
    df_existing = pd.read_csv(output_file, compression='gzip')
    processed_links = set(df_existing['MAPPED_TRAIT_URI'].dropna())
    print(f"Resuming from {len(processed_links)} processed links.")
except FileNotFoundError:
    processed_links = set()

# Dictionary for storing new descriptions
trait_descriptions = {}

# Start WebDriver once
driver = webdriver.Chrome(options=chrome_options)

# Process each link
for index, link in enumerate(tqdm(unique_links)):
    if link in processed_links or link == 'No URI':
        trait_descriptions[link] = 'Description not available'
        continue

    try:
        driver.get(link)
        driver.implicitly_wait(20) # Wait for 5 seconds

        ontology = link.split('/')[-1].split('_')[0]

        if ontology in ols:
            text = driver.find_element(By.CSS_SELECTOR, "p.pb-3").text
        elif ontology in amigo:
            text = driver.find_element(By.XPATH, "//dt[contains(text(), 'Definition')]/following-sibling::dd").text
            text = text.split('\n')[0]
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
        trait_descriptions[link] = f'Error: {str(e)}'

    # Randomized sleep to avoid rate limiting (between 2-5 seconds)
    time.sleep(random.uniform(2, 5))

    # Save progress every 30 links
    if (index + 1) % 30 == 0:
        if 'MAPPED_TRAIT_DESCRIPTION' not in df.columns:
            df['MAPPED_TRAIT_DESCRIPTION'] = ''

        df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(trait_descriptions).fillna(df['MAPPED_TRAIT_DESCRIPTION'])
        df.to_csv(output_file, index=False, compression='gzip')
        print(f"Saved progress at {index + 1} links.")

# Save final data
if 'MAPPED_TRAIT_DESCRIPTION' not in df.columns:
    df['MAPPED_TRAIT_DESCRIPTION'] = ''

df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(trait_descriptions).fillna(df['MAPPED_TRAIT_DESCRIPTION'])
df.to_csv(output_file, index=False, compression='gzip')
print("Completed fetching descriptions.")

# Quit WebDriver
driver.quit()