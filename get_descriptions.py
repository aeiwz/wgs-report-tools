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
df = pd.read_csv('data/gwas_database.csv.gz', compression='gzip', low_memory=False)

# Output file path
output_file = 'data/gwas_database_with_description.csv.gz'

# Set up Chrome options for headless mode
chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--disable-gpu')
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--disable-dev-shm-usage')

# Dictionary for storing descriptions
trait_descriptions = {}

# Start WebDriver once
driver = webdriver.Chrome(options=chrome_options)

# Process all links normally
for index, link in enumerate(tqdm(df['MAPPED_TRAIT_URI'].dropna().unique())):
    if link == 'No URI':
        trait_descriptions[link] = 'Description not available'
        continue

    try:
        if len(link.split(',')) > 1:
            link = link.split(',')[0]
        driver.get(link)
        driver.implicitly_wait(20)  # Wait for 5 seconds

        ontology = link.split('/')[-1].split('_')[0]

        if ontology in ols:
            text = driver.find_element(By.CSS_SELECTOR, "p.pb-3").text
        elif ontology in amigo:
            text = driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div[2]/dl/dd[6]/text()").text
            text = text.split('\n')[0]
        elif ontology in hpo:
            text = driver.find_element(By.XPATH, "/html/body/app-root/mat-sidenav-container/mat-sidenav-content/div/app-term/div/div/div/div[2]/mat-card/mat-card-content/div/div[2]/p").text
            text = text.split('\n')[0]
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

# Save final data after normal processing
if 'MAPPED_TRAIT_DESCRIPTION' not in df.columns:
    df['MAPPED_TRAIT_DESCRIPTION'] = ''

df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(trait_descriptions).fillna(df['MAPPED_TRAIT_DESCRIPTION'])
df.to_csv(output_file, index=False, compression='gzip')
print("Completed normal processing and saved data.")

# Identify failed links (links that encountered "WebDriver error occurred")
failed_links = df[df['MAPPED_TRAIT_DESCRIPTION'] == 'WebDriver error occurred']['MAPPED_TRAIT_URI'].dropna().unique()

# Retry processing only failed links
if failed_links:
    print(f"Re-processing {len(failed_links)} failed links.")

    # Clear previous descriptions for failed links
    failed_trait_descriptions = {}

    # Process each failed link
    for index, link in enumerate(tqdm(failed_links)):
        if link == 'No URI':
            failed_trait_descriptions[link] = 'Description not available'
            continue

        try:
            if len(link.split(',')) > 1:
                link = link.split(',')[0]
            driver.get(link)
            driver.implicitly_wait(20)  # Wait for 5 seconds

            ontology = link.split('/')[-1].split('_')[0]

            if ontology in ols:
                text = driver.find_element(By.CSS_SELECTOR, "p.pb-3").text
            elif ontology in amigo:
                text = driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div[2]/dl/dd[6]/text()").text
                text = text.split('\n')[0]
            elif ontology in hpo:
                text = driver.find_element(By.XPATH, "/html/body/app-root/mat-sidenav-container/mat-sidenav-content/div/app-term/div/div/div/div[2]/mat-card/mat-card-content/div/div[2]/p").text
                text = text.split('\n')[0]
            else:
                text = 'Ontology not covered'

            failed_trait_descriptions[link] = text

        except NoSuchElementException:
            failed_trait_descriptions[link] = 'Description not found'
        except WebDriverException:
            failed_trait_descriptions[link] = 'WebDriver error occurred'
        except Exception as e:
            failed_trait_descriptions[link] = f'Error: {str(e)}'

        # Randomized sleep to avoid rate limiting (between 2-5 seconds)
        time.sleep(random.uniform(2, 5))

    # Save the updated data for failed links
    if 'MAPPED_TRAIT_DESCRIPTION' not in df.columns:
        df['MAPPED_TRAIT_DESCRIPTION'] = ''

    df['MAPPED_TRAIT_DESCRIPTION'] = df['MAPPED_TRAIT_URI'].map(failed_trait_descriptions).fillna(df['MAPPED_TRAIT_DESCRIPTION'])
    df.to_csv(output_file, index=False, compression='gzip')
    print(f"Re-processing complete and saved data.")

# Quit WebDriver
driver.quit()