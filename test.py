import sys 

fname = "C:\\Users\\BMH_work\\github\\image_analysis\\test_out.txt"

print(sys.argv)

with open(fname,'w') as f: 
	f.write(str(len(sys.argv)))
