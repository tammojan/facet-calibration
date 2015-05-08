import os,sys


listname = str(raw_input('Enter list of .ps images to turn into single file '))

list = open(listname,'r')

outname = str(raw_input('Enter outname (.tex added) '))

outfile = open(outname +'.tex','w')

num_per_page = 3
count = 0

number_files = 0
total_count = 0
for line in list:
    number_files += 1


list.seek(0,0)

outfile.write("\documentclass{report}\n")
outfile.write("\usepackage{graphicx} \n")
outfile.write("\usepackage{rotating} \n")
outfile.write("\usepackage{float} \n")
outfile.write("\\begin{document} \n")

for line in list:
    line = line[:-1]
    line = line.split(' ')
    while '' in line:
        line = line.remove('')

    filename = line[0]
    cluster_name = (line[0].split('/'))[1]
    count += 1

    total_count += 1

    if count == 1:
        outfile.write("\\begin{figure}\n")
        outfile.write("\centerline{\includegraphics[width=6.0cm,height= 6.0cm,clip=,angle=0.]{"+filename+"}")
        caption = "{" + cluster_name
        
    if count != 1 and count != num_per_page and total_count != number_files:
        outfile.write("\qquad\includegraphics[width=6.0cm,height= 6.0cm,clip=,angle=0.]{" + filename+ "}")
        caption = caption + ' ' + cluster_name

    if count == num_per_page or total_count == number_files:
        outfile.write("\qquad\includegraphics[width=6.0cm,height= 6.0cm,clip=,angle=0.]{" + filename+ "}")
        outfile.write(" }")
        outfile.write("\caption" + caption + ' ' + cluster_name + '}\n')
        outfile.write("\end{figure}\n\n\n")
        outfile.write(" " )
        count = 0


    if total_count == 24 or total_count == 48 or total_count == 72:
        outfile.write('\\clearpage \n')

       

outfile.write("\end{document} \n")

outfile.close()
os.system('latex ' + outname +'.tex')
os.system('dvips ' + outname  + '.dvi')
os.system('gv ' + outname + '.ps')
