import numpy as np


def bayes(s1, s2):
    numerator = (s1 * s2)
    tmp = (1 - s1) * (1 - s2)
    denominator = numerator + tmp
    return numerator / denominator  # ((s1 * s2) + ((1 - s1) * (1 - s2)))


def bayes_weighted(s1, s2, w1=1.0, w2=1.0):
    if w1 > 1:
        w1 = 1
    elif w1 < 0:
        w1 = 0
    if w2 > 1:
        w2 = 1
    elif w2 < 0:
        w2 = 0
    s1 = s1 * w1 + (1 - w1) / 2.0
    s2 = s2 * w2 + (1 - w2) / 2.0
    return bayes(s1, s2)


def bayes_evidence(posterior, prior=0.5):
    numerator = posterior * (1 - prior)
    denominator = posterior + prior - 2 * posterior * prior
    return numerator / denominator


def merge_mean(scores):
    return np.array(np.mean(np.array(scores), 0))


def merge_bayes(scores):
    s_a = bayes(scores[0], scores[1])
    s_b = bayes(scores[2], scores[3])
    return (s_a + s_b) / 2.0