from PIL import Image
import pyocr
import cv2, re
import itertools
import numpy as np
import matplotlib.pyplot as plt


def codeocr(offset):
    global result
    img = cv2.imread("img_source.png")
    dst = cv2.fastNlMeansDenoisingColored(img, None, 40, 40, 7, 21)  # 去雜點
    ret, thresh = cv2.threshold(dst, 150, 255, cv2.THRESH_BINARY_INV)  # 黑白
    imgarr = cv2.cvtColor(thresh, cv2.COLOR_BGR2GRAY)  # 灰階
    _, inv = cv2.threshold(imgarr, 150, 255, cv2.THRESH_BINARY_INV)  # 轉為反相黑白
    for i in range(len(inv)):  # i為每一列
        for j in range(len(inv[i])):  # j為每一行
            if inv[i][j] == 255:  # 顏色為白色
                count = 0
                for k in range(-1, 2):
                    for l in range(-1, 2):
                        try:
                            if inv[i + k][j + l] == 255:  # 若是白點就將count加1
                                count += 1
                        except IndexError:
                            pass
                if count <= 6:  # 週圍少於等於6個白點
                    inv[i][j] = 255  # 將白點去除
    plt.imshow(inv)
    plt.show()
    # dilation = cv2.dilate(inv, (7, 7), iterations=1)  # 圖形加粗
    dilation = cv2.fastNlMeansDenoising(inv, None, 20, 7, 21)  # 去雜點
    plt.imshow(dilation)
    plt.show()
    # dilation = cv2.dilate(dilation, (7, 7), iterations=1)  # 圖形加粗
    # ret, thresh = cv2.threshold(dilation, 150, 255, cv2.THRESH_BINARY_INV)  # 黑白

    plt.imshow(dilation)
    plt.show()

    cv2.imwrite("final.png", dilation)

    # 文字辨識
    from PIL import Image
    import sys
    import pyocr
    import pyocr.builders

    tools = pyocr.get_available_tools()
    if len(tools) == 0:
        print("No OCR tool found")
        sys.exit(1)
    tool = tools[0]  # 取得可用工具

    result = tool.image_to_string(
        Image.open('final.png'),
        builder=pyocr.builders.TextBuilder()
    )


if __name__ == "__main__":
    offset = 1
    results = []
    while offset < 2:
        print("\noffset=", offset, end="   ")
        offset += 1
        codeocr(offset)
        result = result.replace(" ", "").strip()
        result = re.findall('[a-zA-Z0-9]*', result)[0]
        if len(result) == 4:
            results.append(result)
            print("find:", result)
        else:
            print("no fonud!")

    print("Results=", results)
