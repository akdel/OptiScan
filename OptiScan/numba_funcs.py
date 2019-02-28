import numpy as np
from numba import njit, jit, prange, vectorize


@jit
def get_density_within_window(array, start, stop):
    n = 0
    for i in range(start, stop, 1):
        if array[i]:
            n += 1
        else:
            continue
    return n


@jit
def get_densities_within_windows(array, windows, length):
    res = np.zeros(windows.shape[0], dtype=np.int64)
    for i in range(windows.shape[0]):
        _, idx, _ = windows[i]
        start = idx - length
        stop = idx + length
        res[i] = get_density_within_window(array, start, stop)
    return res


@jit
def fill_chromosomes(chromosome_name, chromosome_index, chromosome_array, vcf_lines):
    vcf_lines = generate_chr_r_q(vcf_lines, chromosome_name)
    for line in vcf_lines:
        chromosome_array[line[1]+np.int(chromosome_index)] = line[3]
    return chromosome_array


@jit
def fill_chromosomes_with_fasta(fasta, chromosome_index, chromosome_array):
    fasta_array = np.array(tuple(fasta), dtype="|S1")
    chromosome_array[int(chromosome_index): int(chromosome_index + len(fasta))] = fasta_array


@vectorize(["int8(int8,int8)"])
def compare_arrays(nuc1, nuc2):
    if nuc1 == 0 and nuc2 == 0:
        return 0
    else:
        prod = nuc1 * nuc2
        if prod != 0:
            return 0
        elif nuc1 >= 0:
            return nuc1
        else:
            return nuc2


@jit
def check_if_exists_forward(seq_array, pattern):
    for i in range(0, seq_array.shape[0] - (20 + pattern.shape[0]), 1):
        if check_if_equal(seq_array[i: i + pattern.shape[0]], pattern):
            return True
        else:
            continue
    return False


@jit
def check_if_exists_forward_v2(seq_array, pattern, pad):
    res = np.zeros((seq_array.shape[0] - (20 + pattern.shape[0] + pad), 20), dtype="int8")
    res[:] = -1
    for i in range(0, seq_array.shape[0] - (20 + pattern.shape[0] + pad), 1):
        if check_if_equal(seq_array[i: i + pattern.shape[0]], pattern):
            res[i] = seq_array[i + pattern.shape[0] + pad: i + pattern.shape[0] + 20 + pad]
        else:
            continue
    return res


@jit
def check_both_directions(seq_array, pattern):
    if check_if_exists_forward(seq_array, pattern):
        return True
    else:
        return check_if_exists_forward(complement(seq_array), pattern)


@jit
def check_all_windows_for_pam(windows, pattern, pad):
    res = np.zeros(windows.shape[0], dtype=np.dtype("(%s,20)i1" % (windows[0].shape[0] - (20 + pattern.shape[0] + pad))))
    for i in range(windows.shape[0]):
        res[i] = check_if_exists_forward_v2(windows[i], pattern, pad)
    return res


def count_lines(vcf_file: str):
    n = 0
    f = open(vcf_file)
    for line in f:
        if not line.startswith("#"):
            n += 1
    return n


def generate_chr_r_q(lines, chromosome=True):
    for line in lines:
        if not line.startswith("#"):
            line = line.strip().split()
            chr = line[0]
            a = np.int(line[1].strip())
            try:
                b = line[3]
                c = line[4]
            except TypeError:
                continue
            if chr == chromosome and len(b) == 1 and len(c) == 1:
                yield chr, a, b, c


@jit
def extender(arr, length, dtype):
    return np.concatenate((arr, np.zeros(length, dtype=dtype)))


def return_kmers_sequence_window(ref, vcf):
    init_length = 1000000
    kmers = np.zeros(init_length, dtype=np.dtype("|S43, i8, |S1"))
    indices = np.where(vcf.view("int8") != 0)
    n = 0
    for i in indices[0]:
        try:
            kmers[n] = (ref[i-21:i+22].view("|S43")[0], i, vcf[i])
        except IndexError:
            kmers = extender(kmers, init_length, np.dtype("|S43, i8, |S1"))
            kmers[n] = (ref[i-21:i+22].view("|S43")[0], i, vcf[i])
        n += 1
    return kmers[:n]


@jit
def return_windows_for_numba(ref, vcf):
    indices = np.where(vcf != 0)[0]
    kmers = np.zeros(indices.shape[0], dtype=np.dtype("(43,)i1"))
    for i in range(indices.shape[0]):
        kmers[i] = ref[indices[i]-21:indices[i]+22]
    return kmers


@jit
def reverse_complement_return_windows_for_numba(ref, vcf):
    indices = np.where(vcf != 0)[0]
    kmers = np.zeros(indices.shape[0], dtype=np.dtype("(43,)i1"))
    for i in range(indices.shape[0]):
        kmers[i] = reverse_complement_v2(ref[indices[i]-21:indices[i]+22])
    return kmers


def kmerize_for_digestion(window_array, max_length):
    window_length = len(window_array[0][0])
    start_i = (window_length//2) - max_length
    end_i = (window_length//2) + max_length
    sequence_window = reduce_array(window_array, start_i, end_i)
    kmers = np.zeros((window_array.shape[0], max_length + 1), dtype=np.dtype("|S%s" % max_length))
    for i in range(sequence_window.shape[0]):
        for r in range(0, max_length + 1, 1):
            kmers[i, r] = sequence_window[i][r: r + max_length]
    return kmers


def check_digestion(kmers, digestion_site):
    assert len(digestion_site) == len(kmers[0][0])
    return np.where(digestion_site == kmers)


def reduce_array(kmer_array, start, end):
    result = np.zeros(kmer_array.shape, dtype=np.dtype("|S%s" % (end-start)))
    for i in range(kmer_array.shape[0]):
        result[i] = kmer_array[i][0][start:end]
    return result

def reverse_complement(fasta_array: np.ndarray):
    reverse_array = np.zeros(len(fasta_array), dtype=np.dtype("|S1"))
    fasta_array = fasta_array[::-1]
    for i in range(fasta_array.shape[0]):
        if fasta_array[i] == "A":
            reverse_array[i] = "T"
        elif fasta_array[i] == "T":
            reverse_array[i] = "A"
        elif fasta_array[i] == "C":
            reverse_array[i] = "G"
        else:
            reverse_array[i] = "C"
    return reverse_array


@jit
def reverse_complement_v2(fasta_array):
    return complement(fasta_array)[::-1]


def to_ord(sequence):
    ord_sequence = np.zeros(len(sequence), dtype=np.int8)
    for i in range(len(sequence)):
        ord_sequence[i] = ord(sequence[i])
    return ord_sequence


@jit(nopython=True)
def complement(sequence):
    res = np.zeros(sequence.shape[0], dtype=np.int8)
    for i in range(sequence.shape[0]):
        if sequence[i] == 65:
            res[i] = 84
        elif sequence[i] == 84:
            res[i] = 65
        elif sequence[i] == 67:
            res[i] = 71
        elif sequence[i] == 71:
            res[i] = 67
        elif sequence[i] == 78:
            res[i] = 78
    return res


@jit
def check_if_equal(arr1, arr2):
    for i in range(arr1.shape[0]):
        if arr1[i] == arr2[i]:
            continue
        else:
            return False
    return True


@njit(nogil=True, parallel=True)
def digest_fasta(fasta_array, digestion_sequence, digestion_array):
    for i in prange(fasta_array.shape[0] - digestion_sequence.shape[0]):
        if check_if_equal(digestion_sequence, fasta_array[i:i+digestion_sequence.shape[0]]):
            digestion_array[i+digestion_sequence.shape[0]+1] = True
        else:
            continue
    return digestion_array


@vectorize(["boolean(boolean,boolean)"])
def return_all(x, y):
    if x or y:
        return True
    else:
        return False


@jit
def digest_with_complementary(fasta_array, digestor):
    res = digest_fasta(fasta_array[0], digestor, np.zeros(fasta_array[0].shape[0], dtype=np.bool))
    res_rev = digest_fasta(fasta_array[1], digestor, np.zeros(fasta_array[0].shape[0], dtype=np.bool))
    return return_all(res, res_rev)


@njit
def sum_bool(x):
    res = 0
    for i in range(x.shape[0]):
        if x[i]:
            res += 1
        else:
            continue
    return res


@jit
def return_densities(steps, digestion_array):
    res = np.zeros((digestion_array.shape[0]//steps) + 1, dtype="int64")
    print(res.shape)
    for i in range(0, digestion_array.shape[0] - steps, steps):
        res[i//steps] = sum_bool(digestion_array[i:i+steps])
    return res


@jit
def return_def_densities(window, digestion_array):
    res = np.zeros(digestion_array.shape[0] - window, dtype="int64")
    res[:window] = sum_bool(digestion_array[:window])
    for i in range(1, digestion_array.shape[0] - window):
        res[i] = res[i-1]
        if digestion_array[window+i]:
            res[i] += 1
        if digestion_array[i - 1]:
            res[i] -= 1
    return res

@jit
def creating_digestion_array_for_resolving(indices, length, windows):
    res = np.zeros(length, dtype=np.bool)
    for i in range(indices.shape[0]):
        ind = indices[i]
        ref_ind = windows[ind][1]
        res[ref_ind] = True
    return res


@jit
def check_snp_density(snp_indices, nick_array, window_size):
    res = np.zeros(snp_indices.shape[0], dtype="int32")
    window_size //= 2
    nick_array = np.pad(nick_array, window_size, False)
    for i in range(snp_indices.shape[0]):
        idx = snp_indices[i]
        res[i] = sum_bool(nick_array[idx - window_size: idx + window_size])
    return res


@jit(nopython=True)
def check_with_mismatches(arr1, arr2, weights, thr):
    n = 0
    mismatches = np.zeros(arr1.shape[0], dtype=np.int8)
    T1 = 1
    for i in range(arr1.shape[0]):
        if n > thr:
            return 0
        elif arr1[i] == arr2[i]:
            pass
        else:
            mismatches[n] = i
            T1 *= 1 - weights[i]
            n += 1
    if n == 0:
        return 1
    elif n == 1:
        d = 19
    else:
        d = (np.max(mismatches[:n]) - np.min(mismatches[:n]))/(n-1)
    T2 = 1/((19 - d)/19 * 4 + 1)
    T3 = 1/(n**2)
    return T1 * T2 * T3
            

@jit(nopython=True, parallel=True)
def return_all_pam_sites(kmers, pattern):
    off_target_weights = np.array([0 ,0 ,0.014 ,0 ,0 ,0.395 ,0.317 ,0 , 0.389, 0.079,
                                   0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583])
    GG = np.array([71,71], dtype=np.int8)
    res = np.zeros(kmers.shape[0], dtype=np.float32)
    for i in prange(res.shape[0]):
        kmer = kmers[i]
        if check_if_equal(kmer[-2:], GG):
            res[i] = check_with_mismatches(kmer[:20], pattern, off_target_weights, 4)
        else:
            continue
    return res

@jit(nopython=True, parallel=True)
def return_all_pam_sites_in_reverse(kmers, pattern):
    off_target_weights = np.array([0 ,0 ,0.014 ,0 ,0 ,0.395 ,0.317 ,0 , 0.389, 0.079,
                                   0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583])[::-1]
    CC = np.array([67,67], dtype=np.int8)
    res = np.zeros(kmers.shape[0], dtype=np.float32)
    pattern_rev = reverse_complement_v2(pattern)
    for i in prange(res.shape[0]):
        kmer = kmers[i]
        if check_if_equal(kmer[:2], CC):
            res[i] = check_with_mismatches(kmer[-20:], pattern_rev, off_target_weights, 4)
        else:
            continue
    return res

@jit(nopython=True, parallel=True)
def get_pam_kmers(kmers):
    GG = np.array([71,71], dtype=np.int8)
    res = np.zeros(kmers.shape[0], dtype=np.int64)
    n = 0
    for i in range(0, kmers.shape[0]):
        kmer = kmers[i]
        if check_if_equal(kmer[-2:], GG):
            res[n] = i
            n += 1
    return res[:n]

@jit(nopython=True, parallel=True)
def pam_kmers_to_kmers(pam_kmers):
    res = np.zeros((pam_kmers.shape[0], 20), dtype=np.int8)
    for i in prange(pam_kmers.shape[0]):
        res[i] = pam_kmers[i][:20]
    return res

@jit(nopython=True, parallel=True)
def zoomed_out_freq(kmers, pattern, times):
    GG = np.array([71,71], dtype=np.int8)
    res = np.zeros(kmers.shape[0]//times, dtype=np.int8)
    for i in prange(res.shape[0]):
        kmer = kmers[i]
        if check_if_equal(kmer[:20], pattern) and check_if_equal(kmer[-2:], GG):
            res[i//times] = 1
        else:
            continue
    return res

@jit
def get_density_to_plot(arr, window_size):
    res = np.zeros(arr.shape[0]//window_size)
    n = 0
    for i in range(0, arr.shape[0], window_size):
        res[n] = np.sum(arr[i:i+window_size])
        n += 1
    return res

@jit(nopython=True)
def get_multiple_pam_sites(kmers, patterns):
    res = np.zeros((patterns.shape[0], kmers.shape[0]), dtype=np.float32)
    for i in range(patterns.shape[0]):
        res[i] = return_all_pam_sites(kmers, patterns[i])
    return res

@jit(nopython=True)
def get_multiple_pam_sites_in_reverse(kmers, patterns):
    res = np.zeros((patterns.shape[0], kmers.shape[0]), dtype=np.float32)
    for i in range(patterns.shape[0]):
        res[i] = return_all_pam_sites_in_reverse(kmers, patterns[i])
    return res

def obtain_top_densities(sequence, rev_seq, top, limit=0.25):
    kmers = np.lib.stride_tricks.as_strided(sequence, shape=(23, sequence.shape[0]-23 + 1), strides=(1, 1)).T
    kmers2 = np.lib.stride_tricks.as_strided(rev_seq, shape=(23, sequence.shape[0]-23 + 1), strides=(1, 1)).T
    #reverse_kmers = np.lib.stride_tricks.as_strided(rev_seq, shape=(23, sequence.shape[0]-23 + 1), strides=(1, 1)).T
    kmers = np.concatenate((kmers,kmers2))
    pam_kmers = get_pam_kmers(kmers.view("int8"))
    pam_kmers = kmers[pam_kmers]
    digestion_kmers = pam_kmers_to_kmers(pam_kmers.view("int8"))
    digestion_kmers = digestion_kmers.view("|S20").flatten()
    counts = np.unique(digestion_kmers, return_counts=True)
    sorted_idx = np.argsort(counts[1])[::-1][:top]
    frequent_kmers = counts[0][sorted_idx].view("|S1").reshape((-1,20)).view("int8")
    limit *= kmers.shape[0]
    limit = int(limit)
    label_sites = get_multiple_pam_sites(kmers.view("int8")[:limit], frequent_kmers)
    label_sites += get_multiple_pam_sites_in_reverse(kmers.view("int8")[:limit], frequent_kmers)
#    label_sites += reverse_label_sites[::-1]
    return frequent_kmers,sorted_idx,label_sites
