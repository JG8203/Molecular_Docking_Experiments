import os
import re
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException
import undetected_chromedriver as uc

def sanitize_filename(filename):
    # Replace any character that is not alphanumeric, underscore, hyphen, or period
    return re.sub(r'[^a-zA-Z0-9_\-.]', '', filename)

def get_compound_name(cid):
    # set options to use Brave browser
    options = uc.ChromeOptions()
    options.binary_location = '/Applications/Brave Browser.app/Contents/MacOS/Brave Browser'  # path to Brave browser on MacOS
    driver = uc.Chrome(options=options)

    driver.get(f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")

    try:
        # Wait for the compound name element to be present
        compound_name_element = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.XPATH, '/html/body/div[1]/div/main/div/div/div/h1'))
        )
        compound_name = compound_name_element.text
    except NoSuchElementException:
        compound_name = ''

    driver.quit()
    return compound_name

directory = os.getcwd()  # Get the current working directory

files = os.listdir(directory)
for filename in files:
    if filename.endswith('.tar.gz'):
        parts = filename.split('_')
        pdb_code = "1HXE"  # hardcoded PDB code
        if len(parts) > 3:
            cid = parts[3]  # CID seems to be at the third index
            compound_name = get_compound_name(cid)
            if compound_name:
                sanitized_name = sanitize_filename(compound_name)
                new_filename = f"{pdb_code}_{sanitized_name}.tar.gz"
                old_path = os.path.join(directory, filename)
                new_path = os.path.join(directory, new_filename)
                os.rename(old_path, new_path)
                print(f"Renamed '{filename}' to '{new_filename}'")
            else:
                print(f"Failed to retrieve compound name for '{filename}'")
        else:
            print(f"Invalid filename format: '{filename}'")
