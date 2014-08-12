import os
import json
import sys 
from datetime import datetime
#from IPython.core.debugger import Tracer
#Tracer()()
#This function makes a new metadata file (metadata_parsed.txt) from metadata.txt for images collected in micromanager in order to extract the times

#imdirbase = 'C:\\Users\\Ben\\Documents\\Data\\PKA_project\\20140709_37_38_39_sorb_p5M\\Pre\\'

imdirbase = sys.argv[1]

pos_names = os.walk(imdirbase).next()[1]

print 'creating simplified metadata file for files in ' + imdirbase

#use Directory structure to cycle through position names.
for pos in pos_names: 
    imdir = pos + '\\'
    
    parsed_data = open(imdirbase+imdir+'metadata_parsed.txt','w')
    parsed_data.write('Position\tChannel\tZstack\tFrame\tTime\n')
    
    json_data=open(imdirbase+imdir+'metadata.txt')
    data = json.load(json_data)
    json_data.close()
    #Get Initial Time and convert to datetime object
    t0_str = data['Summary']['Time']
    fmt = '%Y-%m-%d %H:%M:%S -0700'
    t0 = datetime.strptime(t0_str,fmt)
    
    nFrames = data['Summary']['Frames']
    channels = data['Summary']['ChNames']
    #Depth might be number of z stacks but only BF usually has z-stacks - do if/then 
    #Cycle through all channels, z-levels and frames
    
    for ff in range(0,nFrames):
        for ch_num in range(0,len(channels)): 
            if channels[ch_num] == 'BF':
                for zz in [0,1]:
                    frame_name = 'FrameKey-' + str(ff) + '-' + str(ch_num) +'-' + str(zz)
                    #print frame_name
                    #parsed_data.write('Position Channel Zstack Frame Time')
                    t1_str = data[frame_name]['Time'] 
                    t1 = datetime.strptime(t1_str,fmt)
                    tt = t1 - t0 
                    tt = tt.seconds/60.0
                    parsed_data.write('%s\t%s\t%d\t%d\t%f\n' % (pos,channels[ch_num],zz,ff,tt))                
            else: 
                zz = 1
                frame_name = 'FrameKey-' + str(ff) + '-' + str(ch_num) + '-' + str(zz)    
                #print frame_name
                #parsed_data.write('Position Channel Zstack Frame Time')
                t1_str = data[frame_name]['Time'] 
                t1 = datetime.strptime(t1_str,fmt)
                tt = t1 - t0 
                tt = tt.seconds/60.0
                parsed_data.write('%s\t%s\t%d\t%d\t%f\n' % (pos,channels[ch_num],zz,ff,tt))
                
    
    parsed_data.close()
