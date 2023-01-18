import os
import time

total = 0
correct = 0
failLess100Cand = 0
errorPair = 0
while (total < 1000):
    NBsieve = 53
    sizePremier = 256
    stream = os.popen(f'./genSeq {NBsieve} {sizePremier}')
    output = stream.read()
    lines = output.split('\n')
    p = lines[1].split(': ')[1]
    q = lines[3].split('\t')[1]
    stream = os.popen(f'./analyseSeq output_sequence_p.txt output_sequence_q.txt module_rsa.txt {sizePremier}')
    output = stream.read()
    lines = output.split('\n')
    try:
        out = lines[len(lines)-2]
    except:
        continue
    total = total + 1
    
    if (out != "échec de l'attaque"):
        correct += 1

    print(f"{correct} / {total} ")
    time.sleep(1)
