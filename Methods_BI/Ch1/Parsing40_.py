#! usr/local/bin/python


## sys.stdin.readline().rstrip(), raw_input()

file = raw_input("What file ? ")

infile = open(file, "r")
outfile = open("sample.txt", "w")

for iterations in range(11):
    getline = infile.readline()
    print getline.rstrip()
    outfile.write(getline)

for iterations in range(989):
    getline = infile.readline()
    print getline.rstrip()
    outfile.write(getline)

infile.close()
outfile.close()

print "Your sampled text is in file \"sample.txt\"\n"
exit

