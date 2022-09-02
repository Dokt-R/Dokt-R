from bs4 import BeautifulSoup
import requests
import pandas as pd

# from fake_useragent import UserAgent


def extract(page):
    url = f'https://www.car.gr/used-cars/suzuki/jimny.html?from-models-feed=1&make=12858&model=14897?pg={page}'
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:103.0) Gecko/20100101 Firefox/103.0'}

    r = requests.get(url, headers)
    soup = BeautifulSoup(r.content, 'html.parser')
    return soup


def trasnform(soup):
    li = soup.find_all('li', class_='')
    for item in li:
        model = item.find('h2').text.strip()
        price = item.find('span', class_='price-no-decimals').text
        # Sort to Selling, Leasing or Buying car
        try:
            badge = item.find('span', class_='bg-blue-700').text.strip()
        except:
            try:
                item.find('span', class_='bg-red-700').text.strip()
                badge = 'Buying'
            except:
                badge = 'Selling'
        # Get details since it is hard to get individual details like Year, Distance etc.
        key_features = item.find(
            class_='mt-auto mobile-bottom-left d-flex flex-column align-items-start pr-1').text.strip()
        # Split the key features into manageable list
        feature_list = key_features.split(',')
        year = ''
        print(model, '-->', feature_list)
        car = {
            'model': model,
            'year': year,
            'price': price,
            'badge': badge,
        }
        return
        # carlist.append(car)
    return


carlist = []
for i in range(1, 2):
    c = extract(i)
    t = trasnform(c)
df = pd.DataFrame(carlist)
# print(df.head)
df.to_csv('carlist.csv')
# posting = soup.find('li', class_='')
# title = posting.find('h2', class_='title title')
# print(title)
