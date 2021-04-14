import os
import json
import sys 
from datetime import datetime
from IPython.core.debugger import Tracer
#Tracer()()
#This function makes a new metadata file (metadata_parsed.txt) from metadata.txt for images collected in micromanager in order to extract the times

#The number of frames is set at the beginning of the experiment:  If you have to stop the experiment early, you need to
#use an extra nFrameOverride argument
#
#Example Call from matlab: 
#without override
#system(['python C:/Users/Ben/Documents/GitHub/image_analysis/times_from_umanager_metadata.py C:/Users/Ben/Documents/Data/PKA_project/20160929_SC_NMPP1_Dose_Resp/Exp/')
# 
#with override if only 6 frames were present
#system(['python C:/Users/Ben/Documents/GitHub/image_analysis/times_from_umanager_metadata.py C:/Users/Ben/Documents/Data/PKA_project/20160929_SC_NMPP1_Dose_Resp/Exp/ 6')
#
#change the number of frames in the metadata file (find and replace "Frames": N) There
#is a powershell routine for this named FrameNoShift.ps1

#imdirbase =  os.path.normpath('C:\Users\Ben\Documents\Data\PKA_project\\20160929_SC_NMPP1_Dose_Resp\Exp') + os.sep
imdirbase = os.path.normpath(sys.argv[1]) + os.sep

print(sys.argv)

#nFrameOverride = 6
if len(sys.argv) == 3:
	#nFrameOverride = int(sys.argv[2])
	raise ValueError('len sys.argv is too long'+str(sys.argv))  
else:
    nFrameOverride = 0



#pos_names = next(os.walk(imdirbase))[1]
pos_names = [filename for filename in os.listdir(imdirbase) if os.path.isdir(os.path.join(imdirbase,filename))]


print('creating simplified metadata file for files in ' + imdirbase)

#use Directory structure to cycle through position names.

#first time through take out minimal time in each position, then find the minimum time in all positions and set that as t0

data = {}

#set initial datetime to be greater than any datetimes from the experiment 
t0 = datetime.now()
t0 = t0.replace(year = t0.year + 2)

#dtg format used in micromanager
fmt = '%Y-%m-%d %H:%M:%S'


for pos in pos_names: 
    #Get Initial Time:  out of each position, check to see if it is the earliest time 
    #convert to datetime object also remove timezone from string since that changes and can screw things up. 
    imdir = pos + os.sep
    json_data=open(imdirbase+imdir+'metadata.txt')
    data[pos] = json.load(json_data)
    json_data.close()
    t0_str = data[pos]['Summary']['Time']
    t0_str_sp = t0_str.split()
    t0_str = t0_str_sp[0] + ' ' + t0_str_sp[1]
    t0_new = datetime.strptime(t0_str,fmt)    

    if t0_new < t0: 
        t0 = t0_new
        print( str(t0) + ' ' + str(t0_new))
   

for pos in pos_names: 
    imdir = pos + os.sep
    
    parsed_data = open(imdirbase+imdir+'metadata_parsed.txt','w')
    parsed_data.write('Position\tChannel\tZstack\tFrame\tTime\n')

    if nFrameOverride != 0: 
        nFrames = nFrameOverride    
    else: 
        nFrames = data[pos]['Summary']['Frames']

    channels = data[pos]['Summary']['ChNames']
    #Depth might be number of z stacks but only BF usually has z-stacks - do if/then 
    #Cycle through all channels, z-levels and frames
    #Tracer()()
    BF_present = 'BF' in channels
    for ff in range(0,nFrames):
        for ch_num in range(0,len(channels)): 
            if channels[ch_num] == 'BF':
                for zz in [0,1]:
                    frame_name = 'FrameKey-' + str(ff) + '-' + str(ch_num) +'-' + str(zz)
                    #print frame_name
                    #parsed_data.write('Position Channel Zstack Frame Time')
                    t1_str = data[pos][frame_name]['Time'] 
                    t1_str_sp = t1_str.split()
                    t1_str = t1_str_sp[0] + ' ' + t1_str_sp[1]
                    t1 = datetime.strptime(t1_str,fmt)
                    tt = t1 - t0 
                    tt = tt.seconds/60.0
                    parsed_data.write('%s\t%s\t%d\t%d\t%f\n' % (pos,channels[ch_num],zz,ff,tt))                
            else: 
                if BF_present: 
                    zz = 1
                else:
                    zz = 0
                    
                frame_name = 'FrameKey-' + str(ff) + '-' + str(ch_num) + '-' + str(zz)    
                #print frame_name
                #parsed_data.write('Position Channel Zstack Frame Time')
                t1_str = data[pos][frame_name]['Time'] 
                t1_str_sp = t1_str.split()
                t1_str = t1_str_sp[0] + ' ' + t1_str_sp[1]
                t1 = datetime.strptime(t1_str,fmt)
                tt = t1 - t0 
                tt = tt.seconds/60.0
                parsed_data.write('%s\t%s\t%d\t%d\t%f\n' % (pos,channels[ch_num],zz,ff,tt))
                
    
    parsed_data.close()
