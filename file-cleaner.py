import re

out = open("out.txt", "w")

with open("./data/TV_News_Channel_Commercial_Detection_Dataset/BBC.txt", "r") as f:
    for line in f.readlines():
        clean = " ".join([l.split(":")[1] for l in line.split()[2:]][:17])
        out.write(clean + "\n")

out.close()
