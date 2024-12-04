#After using ape to get a tree in nexus forma the following code was used

with open("mrbayes.run1.t") as f:
    file_to_read= f.readlines()

with open('a_2_mb.t', 'a') as file:
    number=0
    for lines in file_to_read:
        if 'TREE * UNTITLED = [&R]' in lines:
            number+=1
            newline=lines.replace('TREE * UNTITLED = [&R]', 'tree gen.'+str(number)+'= [&U]')
            file.write(newline)
        else:
            file.write(lines)

