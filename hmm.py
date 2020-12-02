import numpy as np


def viterbi(a, pi, b, seq):
    trace = np.zeros((len(b), len(seq) - 1), dtype=int)  # initialization
    delta = pi * b[:, seq[0]]

    for t in range(1, len(seq)):  # steps Viterbi algorithm
        values = delta * (a * b[:, seq[t]]).T

        trace[:, t - 1] = np.argmax(values, axis=1)
        delta = values[np.arange(len(b)), trace[:, t - 1]]

    hidden = np.zeros(len(seq), dtype=int)  # computing way
    hidden[-1] = 1 if delta[0] < delta[1] else 0
    for i in range(len(seq) - 2, -1, -1):
        hidden[i] = trace[hidden[i + 1], i]

    return hidden


def forward_backward(a, pi, b, seq):
    alpha = np.zeros((len(b), len(seq)))
    beta = np.ones_like(alpha)

    alpha[:, 0] = pi * b[:, seq[0]]

    for t in range(1, len(seq)):  # steps forward_backward algorithm len(seq)
        alpha[:, t] = np.sum(alpha[:, t - 1] * (a * b[:, seq[t]]).T, axis=1)
        beta[:, - t - 1] = np.sum(beta[:, - t] * a * b[:, seq[-t]], axis=1)

    return alpha * beta / sum(alpha[:, -1])


if __name__ == "__main__":
    # O - 0, P - 1
    a = np.array([[0.8, 0.2], [0.2, 0.8]])
    pi = np.array([0.5, 0.5])
    b = np.array([[0.5, 0.5], [0.1, 0.9]])
    seq = list(map(int, list("ОРОРОРООРРРРРРРРРРОООООООО".replace("О", "0").replace("Р", "1"))))

    print(viterbi(a, pi, b, seq))  # 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
    print(forward_backward(a, pi, b, seq))

    # [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0]
    # [[0.85956173 0.76609555 0.8783159  0.77939405 0.88686798 0.79512375
    #   0.91061422 0.85999559 0.45699127 0.27318414 0.18984842 0.15315539
    #   0.13940407 0.13971214 0.15427858 0.19251224 0.27910923 0.47000478
    #   0.88850328 0.96461882 0.97845654 0.98093931 0.98120371 0.98022273
    #   0.97438637 0.94221681]
    #  [0.14043827 0.23390445 0.1216841  0.22060595 0.11313202 0.20487625
    #  0.08938578 0.14000441 0.54300873 0.72681586 0.81015158 0.84684461
    #  0.86059593 0.86028786 0.84572142 0.80748776 0.72089077 0.52999522
    #  0.11149672 0.03538118 0.02154346 0.01906069 0.01879629 0.01977727
    #  0.02561363 0.05778319]]