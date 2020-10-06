from selenium import webdriver
from time import sleep
from PIL import Image
from selenium import webdriver


def get_identifying_code(chrome, number, scrollheight=800):
    chrome.execute_script(f"window.scrollTo(0,{scrollheight})")
    chrome.save_screenshot("img_screenshot.png")
    element = chrome.find_elements_by_id('ContentPlaceHolder1_imgcode')[0]

    location = element.location  # 取得圖片 element 的位置
    size = element.size  # 取得圖片 element 的大小

    left = location['x']  # 決定上下左右邊界
    top = location['y'] - scrollheight
    right = location['x'] + size['width']
    bottom = location['y'] - scrollheight + size['height']

    img = Image.open("img_screenshot.png")
    img2 = img.crop((left, top, right, bottom))
    img2.save(f"./xcode/img_source_{number}.png")

def main(url, number):
    chrome = webdriver.Chrome()

    chrome .maximize_window()
    chrome .get(url)

    checkboxs = chrome .find_elements_by_css_selector('input[type=checkbox]')
    for checkbox in checkboxs[-3:]:
        checkbox.click()

    sleep(1)
    chrome.find_element_by_name('ctl00$ContentPlaceHolder1$btnagree').click()
    get_identifying_code(chrome, number)
    chrome.quit()

if __name__ == "__main__":
    url = 'https://npm.cpami.gov.tw/apply_1_2.aspx?unit=c951cdcd-b75a-46b9-8002-8ef952ec95fd'
    # for i in range(100):
    main(url, 0)
