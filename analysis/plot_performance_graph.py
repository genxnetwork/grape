import seaborn as sns
import matplotlib.pyplot as plt


if __name__ == '__main__':

    data_ibis = [12/60, 28/60, 1 + 7/60, 5 + 36/60, 22 + 14/60, 22*60+3]
    data_king = [17/60, 39/60, 1 + 37/60, 8 + 59/60, 31 + 25/60, (22*60+3) + (11*60+39)]

    x_labels = [100, 500, 1000, 2500, 10000, 100000]
    xticks = [1, 4, 8, 8*2.5, 8*2.5*4, 8*2.5*4*10]
    # xticks = list(range(len(x_labels) - 1)) + [6]

    plt.figure(figsize=(18, 12))
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(xticks, [d / 60 for d in data_ibis], label='IBIS + ERSA', marker='o', markersize=18, linewidth=2)
    plt.plot(xticks, [d / 60 for d in data_king], label='IBIS + KING + ERSA', marker='o', markersize=18, linewidth=2)
    plt.grid(True)
    plt.xlabel('Number of Samples', fontsize=24)
    plt.ylabel('Time in Hours', fontsize=24)
    plt.xticks(ticks=xticks, labels=x_labels, fontsize=20)
    plt.yticks(ticks=[1/60, 1, 4, 16], labels=['1m', '1h', '4h', '16h'], fontsize=20)
    plt.legend(fontsize=32)
    plt.ylim(bottom=0.002, top=96)

    text_ibis = ['12s', '28s', '1m 7s', '5m 36s', '22m 14s', '22h 3m 15s']
    text_king = ['17s', '39s', '1m 37s', '8m 59s', '31m 25s', '33h 32m 49s']
    for x, y, txt in zip(xticks, data_ibis, text_ibis):
        plt.annotate(txt, (x - 0.1*x, (y / 60 - 0.37*y / 60)), fontsize=18)

    for x, y, txt in zip(xticks, data_king, text_king):
        plt.annotate(txt, (x - 0.1*x, (y / 60 + 0.3*y / 60)), fontsize=18)
    
    plt.savefig('analysis/performance100k.pdf', bbox_inches='tight')