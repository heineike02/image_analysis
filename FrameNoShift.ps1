#! C:\\Windows\System32\WindowsPowerShell\v1.0\Powershell.exe
#Powershell routine to replace number of frames from 75 to 27 - this is for when you stop your code too early.  
foreach ($file in get-ChildItem *\metadata.txt) { sed -i 's/\"Frames\": 75/\"Frames\": 27/g' $file.FullName }