from bs4 import BeautifulSoup
import requests
import pandas as pd
import matplotlib.pyplot as plt

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
        price = int(item.find(
            'span', class_='price-no-decimals').text.replace('.', ''))
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

        # spans = item.find_all('span')
        # print(spans)

        car = {
            'model': model,
            'year': feature_list[0],
            'distance': feature_list[1],
            'power': feature_list[2],
            'region': feature_list[3],
            'price': price,
            'badge': badge,
        }
        carlist.append(car)
        return
    return


carlist = []
for i in range(1, 2):
    c = extract(i)
    t = trasnform(c)
df = pd.DataFrame(carlist)


# print(df.head)
# df.to_csv('carlist.csv')
# posting = soup.find('li', class_='')
# title = posting.find('h2', class_='title title')
# print(title)

# Plot
plt.style.use('fivethirtyeight')
fig, ax = plt.subplots()

# hide axes
fig.patch.set_visible(False)
ax.axis('tight')
ax.axis('off')

table = ax.table(cellText=df.values, colLabels=df.columns, loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.5, 1.5)

fig.tight_layout()
# Set to full screen
plt.get_current_fig_manager().window.state('zoomed')

# plt.show()
