from PIL import Image
import matplotlib.pyplot as plt



def get_bin_table(threshold=100):
    """
    0 is black. 1 is white.
    :param threshold:
    :return:
    """

    table = []
    for i in range(256):
        if i < threshold:
            table.append(0)
        else:
            table.append(1)
    return table


def sum_9_region_new(img, x, y):
    '''确定噪点 '''
    cur_pixel = img.getpixel((x, y))  # 当前像素点的值
    width = img.width
    height = img.height

    if cur_pixel == 1:  # 如果当前点为白色区域,则不统计邻域值
        return 0

    # 因当前图片的四周都有黑点，所以周围的黑点可以去除
    if y < 3:  # 本例中，前两行的黑点都可以去除
        return 1
    elif y > height - 3:  # 最下面两行
        return 1
    else:  # y不在边界
        if x < 3:  # 前两列
            return 1
        elif x == width - 1:  # 右边非顶点
            return 1
        else:  # 具备9领域条件的
            sum = img.getpixel((x - 1, y - 1)) \
                  + img.getpixel((x - 1, y)) \
                  + img.getpixel((x - 1, y + 1)) \
                  + img.getpixel((x, y - 1)) \
                  + cur_pixel \
                  + img.getpixel((x, y + 1)) \
                  + img.getpixel((x + 1, y - 1)) \
                  + img.getpixel((x + 1, y)) \
                  + img.getpixel((x + 1, y + 1))
            return 9 - sum


def collect_noise_point(img):
    '''收集所有的噪点'''
    noise_point_list = []
    for x in range(img.width):
        for y in range(img.height):
            res_9 = sum_9_region_new(img, x, y)
            if (0 < res_9 < 3) and img.getpixel((x, y)) == 0:  # 找到孤立点
                pos = (x, y)
                noise_point_list.append(pos)
    return noise_point_list


def remove_noise_pixel(img, noise_point_list):
    '''根据噪点的位置信息，消除二值图片的黑点噪声'''
    for item in noise_point_list:
        img.putpixel((item[0], item[1]), 1)


def main():
    image = Image.open('img_source_1.png')
    imgry = image.convert('L')
    table = get_bin_table()
    binary = imgry.point(table, '1')
    # binary.save('binary.png')
    # width, height = binary.size
    # lis = binary.getdata()  # 返回图片所有的像素值，要使用list()才能显示出具体数值
    # lis = list(lis)
    # start = 0
    # step = width
    # for i in range(height):
    #     for p in lis[start: start + step]:
    #         if p == 1:  # 将白色的点变成空格，方便人眼看
    #             p = ' '
    #         print(p, end='')
    #     print('')
    #     start += step
    noise_point_list = collect_noise_point(binary)
    remove_noise_pixel(binary, noise_point_list)
    binary.save('final.png')



if __name__ == '__main__':
    main()
    import cv2
    from PIL import Image
    import sys
    import pyocr
    import pyocr.builders

    img = cv2.imread("final.png")
    dilation = cv2.dilate(img, (9, 9), iterations=1)  # 圖形加粗
    cv2.imwrite("final.png", dilation)


    tools = pyocr.get_available_tools()
    if len(tools) == 0:
        print("No OCR tool found")
        sys.exit(1)
    tool = tools[0]  # 取得可用工具

    result = tool.image_to_string(
        Image.open('final.png'),
        builder=pyocr.builders.TextBuilder()
    )

    print(result)
