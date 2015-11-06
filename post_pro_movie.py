import os
from images2gif import writeGif
from PIL import Image
from time import time 



print('\nStart\n')
subpwd = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/postprotests/wavepacket_008_4modes_(-0.931112136502+0j)/images'

file_names = sorted((fn for fn in os.listdir(subpwd) if fn.endswith('.png')))
subpwd = subpwd + '/'
images = [Image.open(str(subpwd) + fn) for fn in file_names]

runningtime = 5.0
fileName = subpwd + 'movie'
writeGif(fileName + str(int(time())) + ".gif", images, duration=runningtime, dither=1, nq = 1)

print('\nDone\n')




convert -delay 20 *.png a.gif
