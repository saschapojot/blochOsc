import cv2
import glob
import re


def atoi(text):
    return int(text) if text.isdigit() else text

def nkeys(text):
    return [atoi(c) for c in re.split(r"(\d+)",text)]






files=glob.glob("./pump0/q*.png")
files.sort(key=nkeys)
img_arr=[]
for fl in files:
    img=cv2.imread(fl)
    h,w,l=img.shape
    size=(w,h)
    img_arr.append(img)

out=cv2.VideoWriter("p0.mp4",cv2.VideoWriter_fourcc(*'DIVX'), 40, size)
for img in img_arr:
    out.write(img)
out.release()

