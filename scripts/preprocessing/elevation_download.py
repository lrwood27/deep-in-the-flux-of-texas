import requests #request needed downloads
from pathlib import Path #path object creator for files
from urllib.parse import urlparse #breaks up urls

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ELEVATION_DIR = PROJECT_ROOT / 'data' / 'raw' / 'elevation'

urls = []
#text file with all of the download links
with open(ELEVATION_DIR / 'urls.txt', 'r') as file:
    for line in file:
        url = line.strip()

        if url: #check that the line isn't empty
            urls.append(url)


elevation_folder = ELEVATION_DIR / 'tiles' #my folder with elevation data
elevation_folder.mkdir(parents=True, exist_ok=True) #if the directory exists already its ok
c = 1

for url in urls:
    filename = Path(urlparse(url).path).name #creates path object for the file
    destination = elevation_folder / filename #file goes to my elevation folder
    temp = destination.with_suffix('.tif.blah') #so I know if my file is not downloading

    if destination.exists():
        print('already downloaded')
        continue
        
    print(f'[{c}/{len(urls)}] downloading {filename}...')

    with requests.get(url , stream = True ,timeout = 120) as response: #request file at address
        #^^stream=True reads file incrementally
        #connection & read/wait times end at 120 seconds if something is slow
        response.raise_for_status() #raises an error if it can't open the url
        with open(temp , 'wb') as file: #writing my file, wb writes bytes which is what i need
            for chunk in response.iter_content(chunk_size = 1024 * 1024): #1 MiB per chunk
                if chunk: #skips empty chunks
                    file.write(chunk)

    temp.replace(destination) #I rewrite the fake suffix with the actual one if it saved correctly
    print('saved')
    c += 1

print('done')