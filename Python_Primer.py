import json
from pathlib import Path
from time import sleep
from typing import Sequence
from rich import pretty
import re 
import csv
from io import StringIO
import pandas as pd
from selenium.webdriver.remote.utils import dump_json
from selenium.webdriver import Firefox
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

pretty.install()

def make_gibson_primers(NEB_json):
    NEB_json = str(NEB_json)
    opts = Options()
    opts.headless = True
    browser = Firefox(options=opts)
    browser.get("https://nebuilder.neb.com/#!/")
    WebDriverWait(browser, 50).until(EC.element_to_be_clickable((By.CLASS_NAME, "plusfrag")))
    browser.find_elements_by_class_name("plusfrag")[1].click()
    browser.find_element_by_id("projectFileInput").send_keys(NEB_json)
    WebDriverWait(browser, 50).until(EC.element_to_be_clickable((By.CLASS_NAME, "btn-link-cta")))
    browser.find_elements_by_class_name("btn-link-cta")[0].click()
    WebDriverWait(browser, 50).until(EC.element_to_be_clickable((By.CLASS_NAME, "uib-tab")))
    sleep(1)
    browser.find_elements_by_class_name("uib-tab")[1].click()

    with StringIO(browser.find_element_by_id("exportedOligos").text) as f:
        df = pd.read_csv(f, delimiter = ' ',names = ["name", "seq", "concentration", "purif"])
    browser.close()
    return df


base_NEB = Path("/home/thomas/ZapA_his.json")
gibson_templates_file = Path("/home/thomas/Gibson_template.xlsx")

with open(base_NEB) as f:
    base_NEB = json.load(f)


vector = base_NEB["fragments"][0]
insert = base_NEB["fragments"][1]

vector["customFwdPrimer"] = True
vector["customRevPrimer"] = True

gibson_templates = pd.read_excel(gibson_templates_file)

insert_name = "ZapA"
insert_seq = """GTGACATCTGAGAAAAAAACCTACAATTTCCTAATTGCGGGTGTCCCTTACAAATTAAAAACCTCTCATGACGACGCGACTGTTGAAGAGCTCGTTGAATTTGTAAATTCCAGAATGAATCAAGCCCTTGGTGTCACCAAGAATGGCTCCTATCAGAATGCAGCCGTGCTGACGGCGATGAATCTTGCTGAAGAACTTATTCTATTGAAACGCAAAGCCCACCGCGAGTTGGAAAAACTTGAAGAGAAGGCCTTGCAACTTTCAATGGATCTGGAAAATTCCAAGAGCAACAAGGTTTTGAACAACTGA"""
insert_seq = insert_seq.replace(" ", "").replace("\n", "").upper()


result = pd.DataFrame()
for i, row in gibson_templates.iterrows():
    vector["seq"] = row["seq"]
    vector["name"]= row["name"]
    vector["fwdPrimerSeq"] = row["for"]
    vector["revPrimerSeq"] = row["rev"]
    insert["name"] = f"{insert_name}_{row['name']}"

    if row["fusion"] == "C":
        insert["seq"] = insert_seq[:-3]
    elif row["fusion"] == "N":
        insert["seq"] = insert_seq[3:]
    elif row["fusion"] == "NC":
        insert["seq"] = insert_seq[3:-3]

    insert["subseq"] = [1, len(insert["seq"])]
    NEB_json = Path("/home/thomas/NEB.json")
    with open(NEB_json, 'w') as outfile:
        json.dump(base_NEB, outfile)
    
    primers = make_gibson_primers(NEB_json)
    result = result.append(primers[2:])



json = str(json)
opts = Options()
opts.headless = False
browser = Firefox(options=opts)
browser.get("https://nebuilder.neb.com/#!/")
WebDriverWait(browser, 50).until(EC.element_to_be_clickable((By.CLASS_NAME, "plusfrag")))
browser.find_elements_by_class_name("plusfrag")[1].click()
browser.find_element_by_id("projectFileInput").send_keys(str(NEB_json))
WebDriverWait(browser, 50).until(EC.element_to_be_clickable((By.CLASS_NAME, "btn-link-cta")))
browser.find_elements_by_class_name("btn-link-cta")[0].click()
WebDriverWait(browser, 50).until(EC.element_to_be_clickable((By.CLASS_NAME, "uib-tab")))

browser.find_elements_by_class_name("uib-tab")[1].click()

sleep(2)