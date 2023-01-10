import os
import time

total = 0
correct = 0

while (total < 1000):
    NBsieve = 53
    sizePremier = 256
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
        continue
    total = total + 1
    print(f" {premier} vs {premier2}")
    if (int(premier) == int(premier2)):
        correct += 1
    print(f"{correct} / {total}")
    time.sleep(1)
