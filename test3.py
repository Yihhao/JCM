from PIL import Image
import matplotlib.pyplot as plt

def main():
    image = Image.open('xcode/img_source.png')
    imgry = image.convert('L')
    plt.imshow(imgry)
    plt.show()
    imgry.save('gray.png')


if __name__ == '__main__':
    main()