from PIL import Image
from scipy import fftpack
import numpy
import collections
import sys
from bitstream import BitStream
from numpy import *
from collections import OrderedDict

def hexToBytes(hexStr):
    num = len(hexStr)//2
    ret = numpy.zeros([num],dtype=int)
    for i in range(num):
        ret[i] = int(hexStr[2*i:2*i+2],16)
    ret = ret.tolist()
    ret = bytes(ret)
    return ret

#The DC Hoffman coding table for luminance recommended by JPEG
DCLuminanceSizeToCode = [
    [1,1,0],              #0 EOB
    [1,0,1],            #1
    [0,1,1],            #2
    [0,1,0],            #3
    [0,0,0],            #4
    [0,0,1],            #5
    [1,0,0],            #6
    [1,1,1,0],        #7
    [1,1,1,1,0],      #8
    [1,1,1,1,1,0],    #9
    [1,1,1,1,1,1,0],  #10 0A
    [1,1,1,1,1,1,1,0] #11 0B
]

#The DC Hoffman coding table for chrominance recommended by JPEG
DCChrominanceSizeToCode = [
    [0,1],                 #0 EOB
    [0,0],                 #1
    [1,0,0],               #2
    [1,0,1],               #3
    [1,1,0,0],             #4
    [1,1,0,1],             #5
    [1,1,1,0],             #6
    [1,1,1,1,0],           #7
    [1,1,1,1,1,0],     #8
    [1,1,1,1,1,1,0],   #9
    [1,1,1,1,1,1,1,0], #10 0A
    [1,1,1,1,1,1,1,1,0]#11 0B
]

ACLuminanceSizeToCode = {
'01':[0,0],
'02':[0,1],
'03':[1,0,0],
'11':[1,0,1,0],
'04':[1,0,1,1],
'00':[1,1,0,0],
'05':[1,1,0,1,0],
'21':[1,1,0,1,1],
'12':[1,1,1,0,0],
'31':[1,1,1,0,1,0],
'41':[1,1,1,0,1,1],
'51':[1,1,1,1,0,0,0],
'06':[1,1,1,1,0,0,1],
'13':[1,1,1,1,0,1,0],
'61':[1,1,1,1,0,1,1],
'22':[1,1,1,1,1,0,0,0],
'71':[1,1,1,1,1,0,0,1],
'81':[1,1,1,1,1,0,1,0,0],
'14':[1,1,1,1,1,0,1,0,1],
'32':[1,1,1,1,1,0,1,1,0],
'91':[1,1,1,1,1,0,1,1,1],
'A1':[1,1,1,1,1,1,0,0,0],
'07':[1,1,1,1,1,1,0,0,1],
'15':[1,1,1,1,1,1,0,1,0,0],
'B1':[1,1,1,1,1,1,0,1,0,1],
'42':[1,1,1,1,1,1,0,1,1,0],
'23':[1,1,1,1,1,1,0,1,1,1],
'C1':[1,1,1,1,1,1,1,0,0,0],
'52':[1,1,1,1,1,1,1,0,0,1],
'D1':[1,1,1,1,1,1,1,0,1,0],
'E1':[1,1,1,1,1,1,1,0,1,1,0],
'33':[1,1,1,1,1,1,1,0,1,1,1],
'16':[1,1,1,1,1,1,1,1,0,0,0],
'62':[1,1,1,1,1,1,1,1,0,0,1,0],
'F0':[1,1,1,1,1,1,1,1,0,0,1,1],
'24':[1,1,1,1,1,1,1,1,0,1,0,0],
'72':[1,1,1,1,1,1,1,1,0,1,0,1],
'82':[1,1,1,1,1,1,1,1,0,1,1,0,0],
'F1':[1,1,1,1,1,1,1,1,0,1,1,0,1],
'25':[1,1,1,1,1,1,1,1,0,1,1,1,0,0],
'43':[1,1,1,1,1,1,1,1,0,1,1,1,0,1],
'34':[1,1,1,1,1,1,1,1,0,1,1,1,1,0],
'53':[1,1,1,1,1,1,1,1,0,1,1,1,1,1],
'92':[1,1,1,1,1,1,1,1,1,0,0,0,0,0],
'A2':[1,1,1,1,1,1,1,1,1,0,0,0,0,1],
'B2':[1,1,1,1,1,1,1,1,1,0,0,0,1,0,0],
'63':[1,1,1,1,1,1,1,1,1,0,0,0,1,0,1],
'73':[1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0],
'C2':[1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1],
'35':[1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0],
'44':[1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1],
'27':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0],
'93':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,1],
'A3':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,0],
'B3':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1],
'36':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0],
'17':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1],
'54':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,0],
'64':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1],
'74':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0],
'C3':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1],
'D2':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0],
'E2':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,1],
'08':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0],
'26':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1],
'83':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0],
'09':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1],
'0A':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0],
'18':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1],
'19':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,0],
'84':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1],
'94':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0],
'45':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,1],
'46':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,0],
'A4':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1],
'B4':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0],
'56':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,1],
'D3':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0],
'55':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1],
'28':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0],
'1A':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,1],
'F2':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0],
'E3':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1],
'F3':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0],
'C4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,1],
'D4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0],
'E4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1],
'F4':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,0],
'65':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1],
'75':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0],
'85':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1],
'95':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0],
'A5':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1],
'B5':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,0],
'C5':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1],
'D5':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0],
'E5':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1],
'F5':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0],
'66':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],
'76':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
'86':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1],
'96':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0],
'A6':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1],
'B6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0],
'C6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1],
'D6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0],
'E6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1],
'F6':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0],
'37':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,1],
'47':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0],
'57':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1],
'67':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0],
'77':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1],
'87':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0],
'97':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1],
'A7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0],
'B7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1],
'C7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0],
'D7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1],
'E7':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0],
'F7':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1],
'38':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0],
'48':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1],
'58':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0],
'68':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1],
'78':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0],
'88':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1],
'98':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0],
'A8':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1],
'B8':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0],
'C8':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1],
'D8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
'E8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1],
'F8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0],
'29':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1],
'39':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0],
'49':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1],
'59':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0],
'69':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1],
'79':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0],
'89':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1],
'99':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0],
'A9':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1],
'B9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0],
'C9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1],
'D9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0],
'E9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1],
'F9':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
'2A':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1],
'3A':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0],
'4A':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1],
'5A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0],
'6A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1],
'7A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0],
'8A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1],
'9A':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
'AA':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1],
'BA':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0],
'CA':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1],
'DA':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
'EA':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1],
'FA':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]
}
ACLuminanceSizeToCodeList = [
[0,0],
[0,1],
[1,0,0],
[1,0,1,0],
[1,0,1,1],
[1,1,0,0],
[1,1,0,1,0],
[1,1,0,1,1],
[1,1,1,0,0],
[1,1,1,0,1,0],
[1,1,1,0,1,1],
[1,1,1,1,0,0,0],
[1,1,1,1,0,0,1],
[1,1,1,1,0,1,0],
[1,1,1,1,0,1,1],
[1,1,1,1,1,0,0,0],
[1,1,1,1,1,0,0,1],
[1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,0,0,1],
[1,1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,1,0,0,1],
[1,1,1,1,1,1,1,0,1,0],
[1,1,1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,1,1,0,0,1,0],
[1,1,1,1,1,1,1,1,0,0,1,1],
[1,1,1,1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,1,1,1,0,1,1,0,0],
[1,1,1,1,1,1,1,1,0,1,1,0,1],
[1,1,1,1,1,1,1,1,0,1,1,1,0,0],
[1,1,1,1,1,1,1,1,0,1,1,1,0,1],
[1,1,1,1,1,1,1,1,0,1,1,1,1,0],
[1,1,1,1,1,1,1,1,0,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]
]

ACChrominanceToCode = {
'01':[0,0],
'00':[0,1],
'02':[1,0,0],
'11':[1,0,1],
'03':[1,1,0,0],
'04':[1,1,0,1,0],
'21':[1,1,0,1,1],
'12':[1,1,1,0,0,0],
'31':[1,1,1,0,0,1],
'41':[1,1,1,0,1,0],
'05':[1,1,1,0,1,1,0],
'51':[1,1,1,0,1,1,1],
'13':[1,1,1,1,0,0,0],
'61':[1,1,1,1,0,0,1],
'22':[1,1,1,1,0,1,0],
'06':[1,1,1,1,0,1,1,0],
'71':[1,1,1,1,0,1,1,1],
'81':[1,1,1,1,1,0,0,0],
'91':[1,1,1,1,1,0,0,1],
'32':[1,1,1,1,1,0,1,0],
'A1':[1,1,1,1,1,0,1,1,0],
'B1':[1,1,1,1,1,0,1,1,1],
'F0':[1,1,1,1,1,1,0,0,0],
'14':[1,1,1,1,1,1,0,0,1],
'C1':[1,1,1,1,1,1,0,1,0,0],
'D1':[1,1,1,1,1,1,0,1,0,1],
'E1':[1,1,1,1,1,1,0,1,1,0],
'23':[1,1,1,1,1,1,0,1,1,1],
'42':[1,1,1,1,1,1,1,0,0,0],
'15':[1,1,1,1,1,1,1,0,0,1,0],
'52':[1,1,1,1,1,1,1,0,0,1,1],
'62':[1,1,1,1,1,1,1,0,1,0,0],
'72':[1,1,1,1,1,1,1,0,1,0,1],
'F1':[1,1,1,1,1,1,1,0,1,1,0],
'33':[1,1,1,1,1,1,1,0,1,1,1],
'24':[1,1,1,1,1,1,1,1,0,0,0,0],
'34':[1,1,1,1,1,1,1,1,0,0,0,1],
'43':[1,1,1,1,1,1,1,1,0,0,1,0],
'82':[1,1,1,1,1,1,1,1,0,0,1,1],
'16':[1,1,1,1,1,1,1,1,0,1,0,0,0],
'92':[1,1,1,1,1,1,1,1,0,1,0,0,1],
'53':[1,1,1,1,1,1,1,1,0,1,0,1,0],
'25':[1,1,1,1,1,1,1,1,0,1,0,1,1],
'A2':[1,1,1,1,1,1,1,1,0,1,1,0,0],
'63':[1,1,1,1,1,1,1,1,0,1,1,0,1],
'B2':[1,1,1,1,1,1,1,1,0,1,1,1,0],
'C2':[1,1,1,1,1,1,1,1,0,1,1,1,1],
'07':[1,1,1,1,1,1,1,1,1,0,0,0,0,0],
'73':[1,1,1,1,1,1,1,1,1,0,0,0,0,1],
'D2':[1,1,1,1,1,1,1,1,1,0,0,0,1,0],
'35':[1,1,1,1,1,1,1,1,1,0,0,0,1,1,0],
'E2':[1,1,1,1,1,1,1,1,1,0,0,0,1,1,1],
'44':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,0],
'83':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,0],
'17':[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1],
'54':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0],
'93':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1],
'08':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,0],
'09':[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1],
'0A':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0],
'18':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1],
'19':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0],
'26':[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,1],
'36':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0],
'45':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1],
'1A':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0],
'27':[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1],
'64':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0],
'74':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1],
'55':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,0],
'37':[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1],
'F2':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0],
'A3':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,1],
'B3':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,0],
'C3':[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1],
'28':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0],
'29':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,1],
'D3':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0],
'E3':[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1],
'F3':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0],
'84':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,1],
'94':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0],
'A4':[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1],
'B4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0],
'C4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,1],
'D4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0],
'E4':[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1],
'F4':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,0],
'65':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1],
'75':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0],
'85':[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1],
'95':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0],
'A5':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1],
'B5':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,0],
'C5':[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1],
'D5':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0],
'E5':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1],
'F5':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0],
'46':[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],
'56':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
'66':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1],
'76':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0],
'86':[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1],
'96':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0],
'A6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1],
'B6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0],
'C6':[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1],
'D6':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0],
'E6':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,1],
'F6':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0],
'47':[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1],
'57':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0],
'67':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1],
'77':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0],
'87':[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1],
'97':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0],
'A7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1],
'B7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0],
'C7':[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1],
'D7':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0],
'E7':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1],
'F7':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0],
'38':[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1],
'48':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0],
'58':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1],
'68':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0],
'78':[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1],
'88':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0],
'98':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1],
'A8':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0],
'B8':[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1],
'C8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
'D8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1],
'E8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0],
'F8':[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1],
'39':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0],
'49':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1],
'59':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0],
'69':[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1],
'79':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0],
'89':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1],
'99':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0],
'A9':[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1],
'B9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0],
'C9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1],
'D9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0],
'E9':[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1],
'F9':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
'2A':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1],
'3A':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0],
'4A':[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1],
'5A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0],
'6A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1],
'7A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0],
'8A':[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1],
'9A':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
'AA':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1],
'BA':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0],
'CA':[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1],
'DA':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
'EA':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1],
'FA':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
}

ACChrominanceToCodeList = [
[0,0],
[0,1],
[1,0,0],
[1,0,1],
[1,1,0,0],
[1,1,0,1,0],
[1,1,0,1,1],
[1,1,1,0,0,0],
[1,1,1,0,0,1],
[1,1,1,0,1,0],
[1,1,1,0,1,1,0],
[1,1,1,0,1,1,1],
[1,1,1,1,0,0,0],
[1,1,1,1,0,0,1],
[1,1,1,1,0,1,0],
[1,1,1,1,0,1,1,0],
[1,1,1,1,0,1,1,1],
[1,1,1,1,1,0,0,0],
[1,1,1,1,1,0,0,1],
[1,1,1,1,1,0,1,0],
[1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,0,0,1],
[1,1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,1,0,0,1,0],
[1,1,1,1,1,1,1,0,0,1,1],
[1,1,1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,1,1,0,0,0,0],
[1,1,1,1,1,1,1,1,0,0,0,1],
[1,1,1,1,1,1,1,1,0,0,1,0],
[1,1,1,1,1,1,1,1,0,0,1,1],
[1,1,1,1,1,1,1,1,0,1,0,0,0],
[1,1,1,1,1,1,1,1,0,1,0,0,1],
[1,1,1,1,1,1,1,1,0,1,0,1,0],
[1,1,1,1,1,1,1,1,0,1,0,1,1],
[1,1,1,1,1,1,1,1,0,1,1,0,0],
[1,1,1,1,1,1,1,1,0,1,1,0,1],
[1,1,1,1,1,1,1,1,0,1,1,1,0],
[1,1,1,1,1,1,1,1,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]
]

ACChrominance_stats = {}
ACLuminance_stats = {}

#DC components are differentially coded as (SIZE,Value)
def encodeDCToBoolList(value,isLuminance,debugMode = 0):
    boolList = []
    size = int(value).bit_length() # int(0).bit_length()=0
    if(isLuminance==1):
        boolList = boolList + DCLuminanceSizeToCode[size]
    else:
        boolList = boolList + DCChrominanceSizeToCode[size]
    if(value<=0): # if value==0, codeList = [], (SIZE,VALUE)=(SIZE)=EOB
        codeList = list(bin(value)[3:])
        for i in range(len(codeList)):
            if (codeList[i] == '0'):
                codeList[i] = 1
            else:
                codeList[i] = 0
    else:
        codeList = list(bin(value)[2:])
        for i in range(len(codeList)):
            if (codeList[i] == '0'):
                codeList[i] = 0
            else:
                codeList[i] = 1
    boolList = boolList + codeList
    if(debugMode == 1):
        if(isLuminance==1):
            print('isLuminance=',isLuminance,'(size,value)=',size,value,'code=',DCLuminanceSizeToCode[size],codeList)
        else:
            print('isLuminance=', isLuminance, '(size,value)=', size, value, 'code=', DCChrominanceSizeToCode[size],codeList)
    return boolList

def encodeACBlock(bitStream,ACArray,isLuminance,debugMode = 0):

    i = 0
    maxI = numpy.size(ACArray)
    while 1:
        if(i==maxI):
            break
        run = 0

        # check if rest of ACArray are all zero. If so, just write EOB and return
        j = i
        while 1:
            if(ACArray[j]!=0):
                break
            if(j==maxI - 1):
                if (isLuminance == 1):
                    bitStream.write(ACLuminanceSizeToCode['00'], bool)  # EOB
                    if (debugMode == 1):
                        print('EOB', ACLuminanceSizeToCode['00'])
                else:
                    bitStream.write(ACChrominanceToCode['00'], bool)
                    if (debugMode == 1):
                        print('EOB', ACChrominanceToCode['00'])
                return
            j = j + 1

        while 1:
            if(ACArray[i]!=0 or i==maxI - 1 or run==15):
                break
            else:
                run = run + 1
                i = i + 1

        value = ACArray[i]

        if(value==0 and run!=15):
            break # Rest of the components are zeros therefore we simply put the EOB to signify this fact

        size = int(value).bit_length()

        runSizeStr = str.upper(str(hex(run))[2:]) + str.upper(str(hex(size))[2:])

        ### add stats
        if (isLuminance == 1):
            if str(runSizeStr) not in ACLuminance_stats:
                ACLuminance_stats[str(runSizeStr)] = 0
            else:
                ACLuminance_stats[str(runSizeStr)] += 1
        else:
            if str(runSizeStr) not in ACChrominance_stats:
                ACChrominance_stats[str(runSizeStr)] = 0
            else:
                ACChrominance_stats[str(runSizeStr)] += 1

        if (isLuminance == 1):
            bitStream.write(ACLuminanceSizeToCode[runSizeStr], bool)
        else:
            bitStream.write(ACChrominanceToCode[runSizeStr], bool)


        if(value<=0):# if value==0, codeList = [], (SIZE,VALUE)=(SIZE,[])=EOB
            codeList = list(bin(value)[3:])
            for k in range(len(codeList)):
                if (codeList[k] == '0'):
                    codeList[k] = 1
                else:
                    codeList[k] = 0
        else:
            codeList = list(bin(value)[2:])
            for k in range(len(codeList)):
                if (codeList[k] == '0'):
                    codeList[k] = 0
                else:
                    codeList[k] = 1
        bitStream.write(codeList, bool)
        if(debugMode == 1):
            if(isLuminance==1):
                print('isLuminance=',isLuminance,'(run,size,value)=',run,size,value,'code=',ACLuminanceSizeToCode[runSizeStr],codeList)
            else:
                print('isLuminance=',isLuminance,'(run,size,value)=',run,size,value,'code=',ACChrominanceToCode[runSizeStr],codeList)
        i = i + 1


zigzagOrder = numpy.array([0,1,8,16,9,2,3,10,17,24,32,25,18,11,4,5,12,19,26,33,40,48,41,34,27,20,13,6,7,14,21,28,35,42,
                           49,56,57,50,43,36,29,22,15,23,30,37,44,51,58,59,52,45,38,31,39,46,53,60,61,54,47,55,62,63])
#std_quant_tbl from libjpeg::jcparam.c


std_luminance_quant_tbl = numpy.array(
[ 16,  11,  10,  16,  24,  40,  51,  61,
  12,  12,  14,  19,  26,  58,  60,  55,
  14,  13,  16,  24,  40,  57,  69,  56,
  14,  17,  22,  29,  51,  87,  80,  62,
  18,  22,  37,  56,  68, 109, 103,  77,
  24,  35,  55,  64,  81, 104, 113,  92,
  49,  64,  78,  87, 103, 121, 120, 101,
  72,  92,  95,  98, 112, 100, 103,  99],dtype=int)
std_luminance_quant_tbl = std_luminance_quant_tbl.reshape([8,8])

std_chrominance_quant_tbl = numpy.array(
[ 17,  18,  24,  47,  99,  99,  99,  99,
  18,  21,  26,  66,  99,  99,  99,  99,
  24,  26,  56,  99,  99,  99,  99,  99,
  47,  66,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99,
  99,  99,  99,  99,  99,  99,  99,  99],dtype=int)
std_chrominance_quant_tbl = std_chrominance_quant_tbl.reshape([8,8])

def main(tt):

    # inputBMPFileName outputJPEGFilename quality(from 1 to 100) DEBUGMODE(0 or 1)
    # example:
    # ./lena.bmp ./output.jpg 80 0

    if(len(sys.argv)!=5):
        print('inputBMPFileName outputJPEGFilename quality(from 1 to 100) DEBUGMODE(0 or 1)')
        print('example:')
        print('./lena.bmp ./output.jpg 80 0')
        return

    srcFileName = sys.argv[1]
    outputJPEGFileName = sys.argv[2]
    quality = float(sys.argv[3])
    DEBUG_MODE = int(sys.argv[4])

    numpy.set_printoptions(threshold=numpy.inf)
    srcImage = Image.open(srcFileName)
    srcImageWidth, srcImageHeight = srcImage.size
    print('srcImageWidth = %d srcImageHeight = %d' % (srcImageWidth, srcImageHeight))
    print('srcImage info:\n', srcImage)
    srcImageMatrix = numpy.asarray(srcImage)

    imageWidth = srcImageWidth
    imageHeight = srcImageHeight
    # add width and height to %8==0
    if (srcImageWidth % 8 != 0):
        imageWidth = srcImageWidth // 8 * 8 + 8
    if (srcImageHeight % 8 != 0):
        imageHeight = srcImageHeight // 8 * 8 + 8

    print('added to: ', imageWidth, imageHeight)

    # copy data from srcImageMatrix to addedImageMatrix
    addedImageMatrix = numpy.zeros((imageHeight, imageWidth, 3), dtype=numpy.uint8)
    for y in range(srcImageHeight):
        for x in range(srcImageWidth):
            addedImageMatrix[y][x] = srcImageMatrix[y][x]

    # split y u v
    yImage,uImage,vImage = Image.fromarray(addedImageMatrix).convert('YCbCr').split()

    yImageMatrix = numpy.asarray(yImage).astype(int)
    uImageMatrix = numpy.asarray(uImage).astype(int)
    vImageMatrix = numpy.asarray(vImage).astype(int)
    if(DEBUG_MODE==1):
        print('yImageMatrix:\n', yImageMatrix)
        print('uImageMatrix:\n', uImageMatrix)
        print('vImageMatrix:\n', vImageMatrix)

    yImageMatrix = yImageMatrix - 127
    uImageMatrix = uImageMatrix - 127
    vImageMatrix = vImageMatrix - 127

    if(quality <= 0):
        quality = 1
    if(quality > 100):
        quality = 100
    if(quality < 50):
        qualityScale = 5000 / quality
    else:
        qualityScale = 200 - quality * 2
    luminanceQuantTbl = numpy.array(numpy.floor((std_luminance_quant_tbl * qualityScale + 50) / 100))
    luminanceQuantTbl[luminanceQuantTbl == 0] = 1
    luminanceQuantTbl[luminanceQuantTbl > 255] = 255
    luminanceQuantTbl = luminanceQuantTbl.reshape([8, 8]).astype(int)
    print('luminanceQuantTbl:\n', luminanceQuantTbl)
    chrominanceQuantTbl = numpy.array(numpy.floor((std_chrominance_quant_tbl * qualityScale + 50) / 100))
    chrominanceQuantTbl[chrominanceQuantTbl == 0] = 1
    chrominanceQuantTbl[chrominanceQuantTbl > 255] = 255
    chrominanceQuantTbl = chrominanceQuantTbl.reshape([8, 8]).astype(int)
    print('chrominanceQuantTbl:\n', chrominanceQuantTbl)
    blockSum = imageWidth // 8 * imageHeight // 8

    yDC = numpy.zeros([blockSum], dtype=int)
    uDC = numpy.zeros([blockSum], dtype=int)
    vDC = numpy.zeros([blockSum], dtype=int)
    dyDC = numpy.zeros([blockSum], dtype=int)
    duDC = numpy.zeros([blockSum], dtype=int)
    dvDC = numpy.zeros([blockSum], dtype=int)

    print('blockSum = ', blockSum)

    sosBitStream = BitStream()

    blockNum = 0
    for y in range(0, imageHeight, 8):
        for x in range(0, imageWidth, 8):
            # print('block (y,x): ',y, x, ' -> ', y + 8, x + 8)
            yDctMatrix = fftpack.dct(fftpack.dct(yImageMatrix[y:y + 8, x:x + 8], norm='ortho').T, norm='ortho').T
            uDctMatrix = fftpack.dct(fftpack.dct(uImageMatrix[y:y + 8, x:x + 8], norm='ortho').T, norm='ortho').T
            vDctMatrix = fftpack.dct(fftpack.dct(vImageMatrix[y:y + 8, x:x + 8], norm='ortho').T, norm='ortho').T
            if(blockSum<=8):
                print('yDctMatrix:\n',yDctMatrix)
                print('uDctMatrix:\n',uDctMatrix)
                print('vDctMatrix:\n',vDctMatrix)

            yQuantMatrix = numpy.rint(yDctMatrix / luminanceQuantTbl)
            uQuantMatrix = numpy.rint(uDctMatrix / chrominanceQuantTbl)
            vQuantMatrix = numpy.rint(vDctMatrix / chrominanceQuantTbl)
            if(DEBUG_MODE==1):
                print('yQuantMatrix:\n',yQuantMatrix)
                print('uQuantMatrix:\n',uQuantMatrix)
                print('vQuantMatrix:\n',vQuantMatrix)

            yZCode = yQuantMatrix.reshape([64])[zigzagOrder]
            uZCode = uQuantMatrix.reshape([64])[zigzagOrder]
            vZCode = vQuantMatrix.reshape([64])[zigzagOrder]
            yZCode = yZCode.astype(numpy.int)
            uZCode = uZCode.astype(numpy.int)
            vZCode = vZCode.astype(numpy.int)

            yDC[blockNum] = yZCode[0]
            uDC[blockNum] = uZCode[0]
            vDC[blockNum] = vZCode[0]

            if(blockNum==0):
                dyDC[blockNum] = yDC[blockNum]
                duDC[blockNum] = uDC[blockNum]
                dvDC[blockNum] = vDC[blockNum]
            else:
                dyDC[blockNum] = yDC[blockNum] - yDC[blockNum-1]
                duDC[blockNum] = uDC[blockNum] - uDC[blockNum-1]
                dvDC[blockNum] = vDC[blockNum] - vDC[blockNum-1]



            # huffman encode https://www.impulseadventure.com/photo/jpeg-huffman-coding.html
            # encode yDC
            if(DEBUG_MODE==1):
                print("encode dyDC:",dyDC[blockNum])
            sosBitStream.write(encodeDCToBoolList(dyDC[blockNum],1, DEBUG_MODE),bool)
            # encode yAC
            if (DEBUG_MODE == 1):
                print("encode yAC:", yZCode[1:])
            encodeACBlock(sosBitStream, yZCode[1:], 1, DEBUG_MODE)

            # encode uDC
            if(DEBUG_MODE==1):
                print("encode duDC:",duDC[blockNum])
            sosBitStream.write(encodeDCToBoolList(duDC[blockNum],0, DEBUG_MODE),bool)
            # encode uAC
            if (DEBUG_MODE == 1):
                print("encode uAC:", uZCode[1:])
            encodeACBlock(sosBitStream, uZCode[1:], 0, DEBUG_MODE)

            # encode vDC
            if(DEBUG_MODE==1):
                print("encode dvDC:",dvDC[blockNum])
            sosBitStream.write(encodeDCToBoolList(dvDC[blockNum],0, DEBUG_MODE),bool)
            # encode uAC
            if (DEBUG_MODE == 1):
                print("encode vAC:", vZCode[1:])
            encodeACBlock(sosBitStream, vZCode[1:], 0, DEBUG_MODE)

            blockNum = blockNum + 1

    jpegFile = open('number' + str(tt) + '_' + outputJPEGFileName, 'wb+')
    # write jpeg header
    jpegFile.write(hexToBytes('FFD8FFE000104A46494600010100000100010000'))
    # write y Quantization Table
    jpegFile.write(hexToBytes('FFDB004300'))
    luminanceQuantTbl = luminanceQuantTbl.reshape([64])
    jpegFile.write(bytes(luminanceQuantTbl.tolist()))
    # write u/v Quantization Table
    jpegFile.write(hexToBytes('FFDB004301'))
    chrominanceQuantTbl = chrominanceQuantTbl.reshape([64])
    jpegFile.write(bytes(chrominanceQuantTbl.tolist()))
    # write height and width
    jpegFile.write(hexToBytes('FFC0001108'))
    hHex = hex(srcImageHeight)[2:]
    while len(hHex) != 4:
        hHex = '0' + hHex

    jpegFile.write(hexToBytes(hHex))

    wHex = hex(srcImageWidth)[2:]
    while len(wHex) != 4:
        wHex = '0' + wHex

    jpegFile.write(hexToBytes(wHex))

    # 03    01 11 00    02 11 01    03 11 01
    # 1：1	01 11 00	02 11 01	03 11 01
    # 1：2	01 21 00	02 11 01	03 11 01
    # 1：4	01 22 00	02 11 01	03 11 01
    # write Subsamp
    jpegFile.write(hexToBytes('03011100021101031101'))

    ### Huffman table bitstream
    ACLuminanceStr = ""
    for l in ACLuminanceSizeToCodeList:
        for s in ACLuminanceSizeToCode:
            if ACLuminanceSizeToCode[s] == l:
                ACLuminanceStr += (s)
                break

    ACChrominanceStr = ""
    for l in ACChrominanceToCodeList:
        for s in ACChrominanceToCode:
            if ACChrominanceToCode[s] == l:
                ACChrominanceStr += (s)
                break

    #write huffman table
    jpegFile.write(hexToBytes('FFC401A20000000701010101010000000000000000040503020601000708090A0B0100020203010101010100000000000000010002030405060708090A0B1000020103030204020607030402060273'+ACLuminanceStr+'110002020102030505040506040803036D'+ACChrominanceStr))

    # SOS Start of Scan
    # yDC yAC uDC uAC vDC vAC
    sosLength = sosBitStream.__len__()
    filledNum = 8 - sosLength % 8
    if(filledNum!=0):
        sosBitStream.write(numpy.ones([filledNum]).tolist(),bool)

    jpegFile.write(bytes([255, 218, 0, 12, 3, 1, 0, 2, 17, 3, 17, 0, 63, 0])) # FF DA 00 0C 03 01 00 02 11 03 11 00 3F 00

    # write encoded data
    sosBytes = sosBitStream.read(bytes)
    for i in range(len(sosBytes)):
        jpegFile.write(bytes([sosBytes[i]]))
        if(sosBytes[i]==255):
            jpegFile.write(bytes([0])) # FF to FF 00

    # write end symbol
    jpegFile.write(bytes([255,217])) # FF D9
    jpegFile.close()

    ### change ACChrominance
    if tt == 1:
        sorted_chrominance = sorted(ACChrominance_stats.items(), key=lambda item: item[1], reverse = True)
        print(sorted_chrominance)
        cnt = 0
        remaining = []
        for s in ACChrominanceToCode:
            if s == '00':
                continue
            if s not in ACChrominance_stats:
                remaining.append(s)
        for s in sorted_chrominance:
            if cnt == 1:
                ACChrominanceToCode['00'] = ACChrominanceToCodeList[cnt]
                cnt += 1
            ACChrominanceToCode[s[0]] = ACChrominanceToCodeList[cnt]
            cnt += 1
        for s in remaining:
            ACChrominanceToCode[s] = ACChrominanceToCodeList[cnt]
            cnt += 1

    ### change ACLuminance
    if tt == 0:
        sorted_luminance = sorted(ACLuminance_stats.items(), key=lambda item: item[1], reverse = True)    
        cnt = 0
        remaining = []
        for s in ACLuminanceSizeToCode:
            if s == '00':
                continue
            if s not in ACLuminance_stats:
                remaining.append(s)
        for s in sorted_luminance:
            if cnt == 4:
                ACLuminanceSizeToCode['00'] = ACLuminanceSizeToCodeList[cnt]
                cnt += 1
            ACLuminanceSizeToCode[s[0]] = ACLuminanceSizeToCodeList[cnt]
            cnt += 1
        for s in remaining:
            ACLuminanceSizeToCode[s] = ACLuminanceSizeToCodeList[cnt]
            cnt += 1

if __name__ == '__main__':
    main(0)
    main(1)
    main(2)

