#! usr/local/bin/python

file = raw_input("What file? ")
infile = open(file, "r")

while(True):
    for iterations in range(10):
	print infile.readline().rstrip()
    response = raw_input("Type \"QUIT\" if you want to quit. Otherwise press any key\n")
    if response == "QUIT":
	break
infile.close()
exit
    
