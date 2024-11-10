"""
Code to find amino acid sequence in the beta sheet region

Sarah Hui
"""

import requests

url = "https://webclu.bio.wzw.tum.de/stride/stride.tar.gz"
response = requests.get(url)

with open("stride.tar.gz", "wb") as file:
    file.write(response.content)