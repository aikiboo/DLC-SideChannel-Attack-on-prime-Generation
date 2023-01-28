import os
import time


def SPA_func(NBsieve, sizePremier):
    stream = os.popen(f'./genSeq {NBsieve} {sizePremier}')
    output = stream.read()
    lines = output.split('\n')
    p = lines[0].split('\t')[1]
    q = lines[1].split('\t')[1]
    stream = os.popen(f'./analyseSeq output_sequence_p.txt output_sequence_q.txt module_rsa.txt {NBsieve}')
    output = stream.read()
    lines = output.split('\n')
    try:
        out = lines[len(lines) - 2]
    except:
        if True:
            exit()
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
        premier2 = int(lines[len(lines) - 2].split(' ')[2])
    except:
        return CPA_func(NBsieve, sizePremier)
    # print(premier + f" vs  {premier2}")
    return (int(premier) == premier2, int(nbCand))


SPA_RESULT = {(53, 256): 0,
              (59, 256): 0,
              (69, 256): 0}

CPA_RESULT = {(53, 256): {0: [0,0],53: [0,0], 126: [0,0], 246: [0,0]},
              (50, 196): {0: [0,0], 126: [0,0], 246: [0,0]}}

toRun = 1000
total = 1

while (total <= toRun):
    print("\n\n\n\nSPA : ")
    for x in SPA_RESULT:
        SPA_RESULT[x] += SPA_func(x[0], x[1])
        print(f"{x} : {SPA_RESULT[x]}/{total}  ({round((SPA_RESULT[x] / total) * 100, 2)}%)")

    print("CPA : ")
    for x in CPA_RESULT:
        r = CPA_func(x[0], x[1])
        for y in CPA_RESULT[x]:
            if r[1] >= y:
                CPA_RESULT[x][y][0] += r[0]
                CPA_RESULT[x][y][1] += 1
            if CPA_RESULT[x][y][1] >0:
                ratio = round((CPA_RESULT[x][y][0] / CPA_RESULT[x][y][1]) * 100, 2)
                print(f"{x} (+ de {y} candidats) : {CPA_RESULT[x][y][0]}/{CPA_RESULT[x][y][1]}  ({ratio}%)")
    total += 1
    time.sleep(1)
