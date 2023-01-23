import os
import time


def SPA_func(NBsieve, sizePremier):
    stream = os.popen(f'./genSeq {NBsieve} {sizePremier}')
    output = stream.read()
    lines = output.split('\n')
    p = lines[1].split(': ')[1]
    q = lines[3].split('\t')[1]
    stream = os.popen(f'./analyseSeq output_sequence_p.txt output_sequence_q.txt module_rsa.txt {NBsieve}')
    output = stream.read()
    lines = output.split('\n')
    try:
        out = lines[len(lines) - 2]
    except:
        return SPA_func(NBsieve, sizePremier)
    return out != "échec de l'attaque"


def CPA_func(NBsieve, sizePremier):
    stream = os.popen(f'./genDots {NBsieve} {sizePremier}')
    output = stream.read()
    lines = output.split('\n')
    nbCand = lines[1].split(': ')[1]
    premier = lines[2].split('\t')[1]
    stream = os.popen(f'./analyse output_hamming_with_func_and_noise.txt {NBsieve} {nbCand}')
    output = stream.read()
    lines = output.split('\n')
    try:
        premier2 = lines[len(lines) - 2].split(' ')[2]
    except:
        print(NBsieve," ",nbCand)
        exit()
        return CPA_func(NBsieve, sizePremier)
    return int(premier) == int(premier2)


SPA_RESULT = {(53, 256): 0,
              (59, 256): 0,
              (69, 256): 0}

CPA_RESULT = {(53, 256): 0}

toRun = 1000
total = 1

while (total <= toRun):
    print("\n\n\n\nSPA : ")
    for x in SPA_RESULT:
        SPA_RESULT[x] += SPA_func(x[0], x[1])
        print(f"{x} : {SPA_RESULT[x]}/{total}  ({round((SPA_RESULT[x] / total) * 100, 2)}%)")

    print("SPA : ")
    for x in CPA_RESULT:
        CPA_RESULT[x] += CPA_func(x[0], x[1])
        print(f"{x} : {CPA_RESULT[x]}/{total}  ({round((CPA_RESULT[x] / total) * 100, 2)}%)")
    total += 1
    time.sleep(1)
